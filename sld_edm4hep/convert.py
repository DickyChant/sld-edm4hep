"""Convert SLD mini-DST Parquet files (jazelle_reader export) to EDM4hep.

Runs inside a Key4hep environment (podio, edm4hep, ROOT, pyarrow, awkward); see
``scripts/run_key4hep.sh``.  Collection layout follows the ALEPH/OPAL legacy
conversions: EventHeader, ReconstructedParticles, ChargedTracks, Clusters,
PrimaryVertex, MCParticles, plus link collections and a few ``UserDataCollection``s
that keep raw SLD words no EDM4hep field can hold.  Conventions are documented in
``docs/MAPPING.md``; every unit conversion and sign choice is a named constant here.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
import time

import awkward as ak
import numpy as np
import pyarrow.parquet as pq

from . import ptypes

CONVERTER_VERSION = "0.1.0"
MAPPING_VERSION = "2026-09-26"

# --- units and detector constants -------------------------------------------------
CM_TO_MM = 10.0
B_FIELD_TESLA = 0.6  # SLD solenoid
# omega [1/mm] = charge * OMEGA_PER_INVPT * (1/pT [1/GeV])  (LCIO/EDM4hep sign convention)
OMEGA_PER_INVPT = 2.99792458e-4 * B_FIELD_TESLA
LAC_BARREL_RADIUS_MM = 1768.0  # nominal inner radius of the barrel liquid-argon calorimeter
LAC_BARREL_COS_THETA = 0.84  # barrel/endcap boundary
LAC_ENDCAP_Z_MM = 2740.0  # nominal endcap front face

# --- MCPART origin bits (jazelle_reader MCPART::Origin) -----------------------------
ORG_DECAYED = 1 << 0  # decayed by the generator
ORG_DECAYFLT = 1 << 1  # decayed in flight (GEANT)
ORG_STOPPED = 1 << 4
ORG_INTERACT = 1 << 5  # hadronic interaction (GEANT)
ORG_ISR = 1 << 17
ORG_BEAM = 1 << 18
ORG_PREFRAG = 1 << 19
ORG_SWUM = 1 << 20  # tracked through the detector by GEANT

GENSTAT_FINAL = 1
GENSTAT_DECAYED = 2
GENSTAT_BEAM = 4
GENSTAT_PARTON = 21  # quarks, gluons, strings, clusters (documented choice)

# MCHEAD.origin flag -> PDG of the primary quark
MCHEAD_FLAVOUR = {1 << 0: 2, 1 << 1: 1, 1 << 2: 3, 1 << 3: 4, 1 << 4: 5, 1 << 5: 6}

# packed lower-triangle index for symmetric matrices (both SLD and EDM4hep use it)
def _tri(i: int, j: int) -> int:
    if i < j:
        i, j = j, i
    return i * (i + 1) // 2 + j


# SLD helix covariance parameters (phi, 1/pT, tanLambda, TSI, ETA) -> EDM4hep TrackState
# covariance slots (d0, phi, omega, z0, tanLambda, time).  TSI (along the normal of the
# trajectory, in the bending plane) is taken as d0 and ETA (along the binormal) as z0.
_SLD_TO_EDM_SLOT = {0: 1, 1: 2, 2: 4, 3: 0, 4: 3}
_SLD_UNIT_SCALE = {0: 1.0, 1: None, 2: 1.0, 3: CM_TO_MM, 4: CM_TO_MM}  # None: omega scale


class Bank:
    """Flattened numpy view of a jagged ``list<struct>`` Parquet column."""

    def __init__(self, arr: ak.Array | None):
        self.present = arr is not None
        if not self.present:
            self.n = None
            return
        counts = ak.to_numpy(ak.num(arr))
        self.n = counts.astype(np.int64)
        self.off = np.concatenate([[0], np.cumsum(self.n)])
        flat = ak.flatten(arr)
        self.cols = {f: ak.to_numpy(flat[f]) for f in flat.fields}

    def rows(self, ev: int) -> slice:
        return slice(int(self.off[ev]), int(self.off[ev + 1]))

    def count(self, ev: int) -> int:
        return int(self.n[ev]) if self.present else 0


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _id_map(ids: np.ndarray) -> dict[int, int]:
    return {int(v): k for k, v in enumerate(ids)}


def convert(inp: str, out: str, max_events: int | None = None, first: int = 0, quiet: bool = False) -> dict:
    import edm4hep
    import podio
    import ROOT
    from podio import root_io

    t0 = time.time()
    table = pq.read_table(inp)
    ev_arr = ak.from_arrow(table)
    n_total = len(ev_arr)
    last = n_total if max_events is None else min(n_total, first + max_events)

    def bank(name: str) -> Bank:
        return Bank(ev_arr[name] if name in ev_arr.fields else None)

    ieventh = ev_arr["IEVENTH"]
    phbm, phpsum, phpoint = bank("PHBM"), bank("PHPSUM"), bank("PHPOINT")
    phchrg, phklus = bank("PHCHRG"), bank("PHKLUS")
    phkelid, phwic, phcrid = bank("PHKELID"), bank("PHWIC"), bank("PHCRID")
    mchead, mcpart, mcpnt = bank("MCHEAD"), bank("MCPART"), bank("MCPNT")
    is_mc = mcpart.present

    # link collections (EDM4hep >= 1.0); their headers are not loaded by the python module
    for hdr_name in ("TrackMCParticleLinkCollection", "RecoMCParticleLinkCollection"):
        ROOT.gInterpreter.Declare(f'#include "edm4hep/{hdr_name}.h"')
    TrackMCLink = ROOT.edm4hep.TrackMCParticleLinkCollection
    RecoMCLink = ROOT.edm4hep.RecoMCParticleLinkCollection

    def set_link(link, reco, sim, weight):
        if hasattr(link, "setFrom"):
            link.setFrom(reco)
            link.setTo(sim)
        else:
            link.setRec(reco)
            link.setSim(sim)
        link.setWeight(float(weight))

    UD_I32 = ROOT.podio.UserDataCollection["int32_t"]
    UD_F = ROOT.podio.UserDataCollection["float"]

    primary_bit = getattr(edm4hep.Vertex, "BITPrimaryVertex", None)
    vertex_type_primary = (1 << primary_bit) if primary_bit is not None else 1

    writer = root_io.Writer(out)
    stats = {"events": 0, "tracks": 0, "clusters": 0, "reco": 0, "mc": 0, "unknown_ptype": {}}

    for ev in range(first, last):
        frame = podio.Frame()
        hdr_row = ieventh[ev]

        # ---------------- EventHeader ----------------
        headers = edm4hep.EventHeaderCollection()
        hdr = headers.create()
        hdr.setRunNumber(int(hdr_row["run"]))
        hdr.setEventNumber(int(hdr_row["event"]))
        hdr.setTimeStamp(int(hdr_row["evttime"]))
        hdr.setWeight(float(hdr_row["weight"]))
        frame.put_parameter("trigger", int(hdr_row["trigger"]))
        frame.put_parameter("evttype", int(hdr_row["evttype"]))
        frame.put_parameter("header_word", int(hdr_row["header"]))

        # ---------------- PHBM: beam parameters and primary vertex ----------------
        vertices = edm4hep.VertexCollection()
        if phbm.count(ev) > 0:
            r = phbm.off[ev]
            frame.put_parameter("ecm", float(phbm.cols["ecm"][r]))
            frame.put_parameter("decm", float(phbm.cols["decm"][r]))
            frame.put_parameter("pol", float(phbm.cols["pol"][r]))
            frame.put_parameter("dpol", float(phbm.cols["dpol"][r]))
            pos = phbm.cols["pos"][r]
            dpos = phbm.cols["dpos"][r]
            vtx = vertices.create()
            vtx.setType(vertex_type_primary)
            vtx.setPosition(edm4hep.Vector3f(*(float(x) * CM_TO_MM for x in pos)))
            cov = edm4hep.CovMatrix3f()
            for k in range(6):
                cov.values[k] = float(dpos[k]) * CM_TO_MM * CM_TO_MM
            vtx.setCovMatrix(cov)
            vtx.setAlgorithmType(0)

        # ---------------- MCHEAD ----------------
        mc_ip = (0.0, 0.0, 0.0)
        if mchead.count(ev) > 0:
            r = mchead.off[ev]
            origin_word = int(mchead.cols["origin"][r])
            frame.put_parameter("mc_ntot", int(mchead.cols["ntot"][r]))
            frame.put_parameter("mc_origin", origin_word)
            frame.put_parameter("mc_primary_quark_pdg", MCHEAD_FLAVOUR.get(origin_word, 0))
            mc_ip = tuple(float(mchead.cols[k][r]) * CM_TO_MM for k in ("ipx", "ipy", "ipz"))
            frame.put_parameter("mc_ip_mm", list(mc_ip))

        # ---------------- MCParticles ----------------
        mcs = edm4hep.MCParticleCollection()
        mc_type_ud, mc_origin_ud, mc_id_ud = UD_I32(), UD_I32(), UD_I32()
        mc_index_by_id: dict[int, int] = {}
        if is_mc and mcpart.count(ev) > 0:
            s = mcpart.rows(ev)
            c = mcpart.cols
            ids = c["id"][s]
            mc_index_by_id = _id_map(ids)
            px, py, pz, e = (c[k][s].astype(np.float64) for k in ("px", "py", "pz", "e"))
            xt = np.stack([c["xt_x"][s], c["xt_y"][s], c["xt_z"][s]], axis=1).astype(np.float64) * CM_TO_MM
            code = c["ptype"][s]
            origin = c["origin"][s]
            parent_id = c["parent_id"][s]
            p2 = px * px + py * py + pz * pz
            generator_level = np.zeros(len(ids), dtype=bool)
            mc_objs = []  # mutable handles (collection indexing returns const objects)
            for i in range(len(ids)):
                pdg, mass, charge, _name = ptypes.lookup(int(code[i]))
                acode = abs(int(code[i]))
                if acode in ptypes.MASS_FROM_ENERGY or acode not in ptypes.TABLE:
                    mass = math.sqrt(max(e[i] * e[i] - p2[i], 0.0))
                if acode not in ptypes.TABLE:
                    stats["unknown_ptype"][int(code[i])] = stats["unknown_ptype"].get(int(code[i]), 0) + 1
                org = int(origin[i])
                pidx = mc_index_by_id.get(int(parent_id[i]), -1)
                if pidx < 0:
                    gen = True
                else:
                    porg = int(origin[pidx])
                    gen = bool(generator_level[pidx]) and bool(porg & ORG_DECAYED) and not (porg & (ORG_DECAYFLT | ORG_INTERACT))
                generator_level[i] = gen

                p = mcs.create()
                mc_objs.append(p)
                p.setPDG(int(pdg))
                p.setMass(float(mass))
                p.setCharge(float(charge))
                p.setMomentum(edm4hep.Vector3d(px[i], py[i], pz[i]))
                p.setEndpoint(edm4hep.Vector3d(*xt[i]))
                if pidx >= 0:
                    p.setVertex(edm4hep.Vector3d(*xt[pidx]))
                    p.addToParents(mcs[pidx])
                    mc_objs[pidx].addToDaughters(p)
                else:
                    p.setVertex(edm4hep.Vector3d(*mc_ip))
                if gen:
                    if acode in ptypes.MASS_FROM_ENERGY or org & ORG_PREFRAG:
                        gs = GENSTAT_PARTON
                    elif org & ORG_BEAM:
                        gs = GENSTAT_BEAM
                    elif org & ORG_DECAYED:
                        gs = GENSTAT_DECAYED
                    else:
                        gs = GENSTAT_FINAL
                    p.setGeneratorStatus(gs)
                else:
                    p.setGeneratorStatus(0)
                    p.setCreatedInSimulation(True)
                if org & ORG_STOPPED:
                    p.setStopped(True)
                mc_type_ud.push_back(int(code[i]))
                mc_origin_ud.push_back(org)
                mc_id_ud.push_back(int(ids[i]))
            stats["mc"] += len(ids)

        # ---------------- ChargedTracks ----------------
        tracks = edm4hep.TrackCollection()
        dedx_ud, imc_ud = UD_I32(), UD_I32()
        b2d_ud, b2dsig_ud, b3d_ud, b3dsig_ud = UD_F(), UD_F(), UD_F(), UD_F()
        track_links = TrackMCLink()
        trk_index_by_id: dict[int, int] = {}
        if phchrg.count(ev) > 0:
            s = phchrg.rows(ev)
            c = phchrg.cols
            trk_index_by_id = _id_map(c["id"][s])
            hlx = c["hlxpar"][s].astype(np.float64)
            dhlx = c["dhlxpar"][s].astype(np.float64)
            q = c["charge"][s].astype(np.float64)
            for i in range(len(q)):
                t = tracks.create()
                t.setType(int(c["status"][s][i]))
                t.setChi2(float(c["chi2dt"][s][i]))
                t.setNdf(int(c["ndfdt"][s][i]))
                t.setNholes(int(c["nmisht"][s][i]))
                for key in ("nhit", "nhitv", "nhite", "nhitp", "nmisht", "nwrght"):
                    t.addToSubdetectorHitNumbers(int(c[key][s][i]))
                omega_scale = q[i] * OMEGA_PER_INVPT
                st = edm4hep.TrackState()
                st.location = edm4hep.TrackState.AtIP
                st.phi = float(hlx[i, 0])
                st.omega = float(omega_scale * hlx[i, 1])
                st.tanLambda = float(hlx[i, 2])
                st.D0 = 0.0
                st.Z0 = 0.0
                st.referencePoint = edm4hep.Vector3f(*(float(x) * CM_TO_MM for x in hlx[i, 3:6]))
                cov = edm4hep.CovMatrix6f()
                for a in range(5):
                    sa = omega_scale if _SLD_UNIT_SCALE[a] is None else _SLD_UNIT_SCALE[a]
                    for b in range(a + 1):
                        sb = omega_scale if _SLD_UNIT_SCALE[b] is None else _SLD_UNIT_SCALE[b]
                        cov.values[_tri(_SLD_TO_EDM_SLOT[a], _SLD_TO_EDM_SLOT[b])] = float(sa * sb * dhlx[i, _tri(a, b)])
                st.covMatrix = cov
                t.addToTrackStates(st)
                dedx_ud.push_back(int(c["dedx"][s][i]))
                imc_raw = int(c["imc"][s][i])
                imc_ud.push_back(imc_raw)
                b2d_ud.push_back(float(c["bnorm"][s][i]) * CM_TO_MM)
                b2dsig_ud.push_back(float(c["impact"][s][i]))
                b3d_ud.push_back(float(c["b3norm"][s][i]) * CM_TO_MM)
                b3dsig_ud.push_back(float(c["impact3"][s][i]))
                if is_mc and imc_raw != 0:
                    j = mc_index_by_id.get(abs(imc_raw), -1)
                    if j >= 0:
                        lk = track_links.create()
                        set_link(lk, t, mcs[j], 1.0 if imc_raw > 0 else -1.0)
            stats["tracks"] += len(q)

        # ---------------- Clusters ----------------
        clusters = edm4hep.ClusterCollection()
        klus_index_by_id: dict[int, int] = {}
        if phklus.count(ev) > 0:
            s = phklus.rows(ev)
            c = phklus.cols
            klus_index_by_id = _id_map(c["id"][s])
            for i in range(phklus.count(ev)):
                cl = clusters.create()
                cl.setType(int(c["status"][s][i]))
                cl.setEnergy(float(c["eraw"][s][i]))
                cth = float(c["cth"][s][i])
                phi = float(c["phi"][s][i])
                sth = math.sqrt(max(1.0 - cth * cth, 1e-12))
                if abs(cth) <= LAC_BARREL_COS_THETA:
                    rho = LAC_BARREL_RADIUS_MM
                    z = rho * cth / sth
                else:
                    z = math.copysign(LAC_ENDCAP_Z_MM, cth)
                    rho = abs(z) * sth / abs(cth)
                cl.setPosition(edm4hep.Vector3f(rho * math.cos(phi), rho * math.sin(phi), z))
                cl.setITheta(math.acos(max(-1.0, min(1.0, cth))))
                cl.setIPhi(phi)
                for key in ("cth", "wcth", "phi", "wphi", "nhit2", "cth2", "wcth2", "phi2", "whphi2", "nhit3", "cth3", "wcth3", "phi3", "wphi3"):
                    cl.addToShapeParameters(float(c[key][s][i]))
                for x in c["elayer"][s][i]:
                    cl.addToSubdetectorEnergies(float(x))
            stats["clusters"] += phklus.count(ev)

        # ---------------- ReconstructedParticles (PHPSUM + PHPOINT) ----------------
        recos = edm4hep.ReconstructedParticleCollection()
        e_ids = edm4hep.ParticleIDCollection()
        mu_ids = edm4hep.ParticleIDCollection()
        crid_ids = edm4hep.ParticleIDCollection()
        reco_links = RecoMCLink()
        reco_index_by_id: dict[int, int] = {}
        if phpsum.count(ev) > 0:
            s = phpsum.rows(ev)
            c = phpsum.cols
            reco_index_by_id = _id_map(c["id"][s])
            has_point = phpoint.count(ev) == phpsum.count(ev)
            sp = phpoint.rows(ev) if has_point else None
            pc = phpoint.cols
            pid_banks = []
            for coll, bnk, key, algo in (
                (e_ids, phkelid, "phkelid_id", 1),
                (mu_ids, phwic, "phwic_id", 2),
                (crid_ids, phcrid, "phcrid_id", 3),
            ):
                if bnk.present and key in pc and bnk.count(ev) > 0:
                    bs = bnk.rows(ev)
                    pid_banks.append((coll, bnk, key, algo, bs, _id_map(bnk.cols["id"][bs])))
            for i in range(phpsum.count(ev)):
                rp = recos.create()
                mom = (float(c["px"][s][i]), float(c["py"][s][i]), float(c["pz"][s][i]))
                rp.setMomentum(edm4hep.Vector3f(*mom))
                rp.setEnergy(float(c["ptot"][s][i]))
                rp.setMass(0.0)
                rp.setCharge(float(c["charge"][s][i]))
                rp.setPDG(0)
                rp.setGoodnessOfPID(0.0)
                rp.setReferencePoint(edm4hep.Vector3f(*(float(c[k][s][i]) * CM_TO_MM for k in ("x", "y", "z"))))
                if has_point:
                    j = trk_index_by_id.get(int(pc["phchrg_id"][sp][i]), -1)
                    if j >= 0:
                        rp.addToTracks(tracks[j])
                    j = klus_index_by_id.get(int(pc["phklus_id"][sp][i]), -1)
                    if j >= 0:
                        rp.addToClusters(clusters[j])
                    for coll, bnk, key, algo, bs, idmap in pid_banks:
                        bid = int(pc[key][sp][i])
                        if bid <= 0:
                            continue
                        idx = idmap.get(bid, -1)
                        if idx < 0:
                            continue
                        pid = coll.create()
                        pid.setAlgorithmType(algo)
                        pid.setParticle(rp)
                        _fill_pid(pid, bnk, bs, idx, algo)
            stats["reco"] += phpsum.count(ev)

        # RP <-> MC links from MCPNT (weight = energy contribution)
        if is_mc and mcpnt.count(ev) > 0 and reco_index_by_id:
            s = mcpnt.rows(ev)
            c = mcpnt.cols
            for i in range(mcpnt.count(ev)):
                j = mc_index_by_id.get(int(c["mcpart_id"][s][i]), -1)
                k = reco_index_by_id.get(int(c["phpoint_id"][s][i]), -1)
                if j >= 0 and k >= 0:
                    lk = reco_links.create()
                    set_link(lk, recos[k], mcs[j], float(c["econtrib"][s][i]))

        # ---------------- store ----------------
        frame.put(headers, "EventHeader")
        frame.put(vertices, "PrimaryVertex")
        frame.put(recos, "ReconstructedParticles")
        frame.put(tracks, "ChargedTracks")
        frame.put(clusters, "Clusters")
        frame.put(e_ids, "ElectronID")
        frame.put(mu_ids, "MuonID")
        frame.put(crid_ids, "CRID")
        frame.put(dedx_ud, "ChargedTracksDedxRaw")
        frame.put(imc_ud, "ChargedTracksImcRaw")
        frame.put(b2d_ud, "ChargedTracksImpact2D")
        frame.put(b2dsig_ud, "ChargedTracksImpact2DSignificance")
        frame.put(b3d_ud, "ChargedTracksImpact3D")
        frame.put(b3dsig_ud, "ChargedTracksImpact3DSignificance")
        if is_mc:
            frame.put(mcs, "MCParticles")
            frame.put(mc_type_ud, "MCParticlesSLDType")
            frame.put(mc_origin_ud, "MCParticlesSLDOrigin")
            frame.put(mc_id_ud, "MCParticlesSLDId")
            frame.put(track_links, "ChargedTracksMCLinks")
            frame.put(reco_links, "ReconstructedParticlesMCLinks")
        writer.write_frame(frame, "events")
        stats["events"] += 1
        if not quiet and (stats["events"] % 1000 == 0 or ev == last - 1):
            print(f"  {stats['events']:7d} events  {time.time() - t0:7.1f} s", flush=True)

    # ---------------- file-level metadata ----------------
    meta = podio.Frame()
    meta.put_parameter("sld_edm4hep_version", CONVERTER_VERSION)
    meta.put_parameter("mapping_version", MAPPING_VERSION)
    meta.put_parameter("input_file", os.path.basename(inp))
    meta.put_parameter("input_sha256", _sha256(inp))
    meta.put_parameter("is_simulation", int(is_mc))
    meta.put_parameter("length_unit", "mm")
    meta.put_parameter("energy_unit", "GeV")
    meta.put_parameter("b_field_tesla", B_FIELD_TESLA)
    meta.put_parameter("subdetector_hit_numbers", ["nhit", "nhitv", "nhite", "nhitp", "nmisht", "nwrght"])
    meta.put_parameter("cluster_shape_parameters", ["cth", "wcth", "phi", "wphi", "nhit2", "cth2", "wcth2", "phi2", "whphi2", "nhit3", "cth3", "wcth3", "phi3", "wphi3"])
    meta.put_parameter("particle_id_algorithms", ["1 PHKELID electron", "2 PHWIC muon", "3 PHCRID cherenkov"])
    meta.put_parameter("sld_ptype_codes", [int(k) for k in sorted(ptypes.TABLE)])
    meta.put_parameter("sld_ptype_pdg", [int(ptypes.TABLE[k].pdg) for k in sorted(ptypes.TABLE)])
    writer.write_frame(meta, "metadata")
    if hasattr(writer, "close"):
        writer.close()
    else:
        del writer
    stats["seconds"] = round(time.time() - t0, 1)
    return stats


def _fill_pid(pid, bnk: Bank, rows: slice, idx: int, algo: int) -> None:
    c = bnk.cols
    if algo == 1:  # PHKELID: electron identification
        pid.setType(int(c["idstat"][rows][idx]))
        pid.setLikelihood(float(c["prob"][rows][idx]))
        pid.setPDG(11)
        for key in ("qp", "em1", "em12", "had1", "isolat", "phiwid", "thewid", "em1x1", "em2x2a", "em2x2b", "em3x3a", "em3x3b"):
            pid.addToParameters(float(c[key][rows][idx]))
    elif algo == 2:  # PHWIC: muon identification in the warm iron calorimeter
        pid.setType(int(c["idstat"][rows][idx]))
        pid.setLikelihood(float(c["missprob"][rows][idx]))
        pid.setPDG(13)
        for key in ("nhit", "nhit45", "npat", "nhitpat", "nlayexp", "nlaybey", "chi2", "ndf", "matchChi2", "matchNdf", "itrlen"):
            pid.addToParameters(float(c[key][rows][idx]))
    else:  # PHCRID: Cherenkov log-likelihoods e, mu, pi, K, p
        pid.setType(int(c["ctlword"][rows][idx]))
        lls = [float(c[k][rows][idx]) for k in ("llik_e", "llik_mu", "llik_pi", "llik_k", "llik_p")]
        best = int(np.argmax(lls))
        pid.setPDG([11, 13, 211, 321, 2212][best])
        pid.setLikelihood(lls[best])
        for v in lls:
            pid.addToParameters(v)
        for key in ("liq_besthyp", "gas_besthyp", "liq_nhits", "gas_nhits"):
            pid.addToParameters(float(c[key][rows][idx]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", help="SLD mini-DST Parquet file")
    ap.add_argument("output", help="EDM4hep ROOT file to write")
    ap.add_argument("--max-events", type=int, default=None)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"converting {args.input} -> {args.output}", flush=True)
    stats = convert(args.input, args.output, args.max_events, args.first, args.quiet)
    print("done:", stats, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
