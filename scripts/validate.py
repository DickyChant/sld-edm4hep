#!/usr/bin/env python
"""Read a converted EDM4hep file back and compare it with the source Parquet.

    scripts/run_key4hep.sh python scripts/validate.py OUT.root IN.parquet [--first N]

Checks event alignment, object counts, kinematics round-trips (track pT from omega,
reconstructed-particle momenta, cluster energies), MC tree consistency and the
physics sanity of the MC charge convention (net charge of generator-level final-state
particles).  Exits non-zero on a hard failure.
"""

from __future__ import annotations

import argparse
import sys

import awkward as ak
import numpy as np
import pyarrow.parquet as pq
import uproot

sys.path.insert(0, ".")
from sld_edm4hep.convert import CM_TO_MM, OMEGA_PER_INVPT


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("parquet")
    ap.add_argument("--first", type=int, default=0)
    args = ap.parse_args()

    ev = uproot.open(args.root)["events"]
    n = ev.num_entries
    src = ak.from_arrow(pq.read_table(args.parquet))[args.first : args.first + n]
    fails: list[str] = []

    def check(cond: bool, msg: str) -> None:
        print(("ok   " if cond else "FAIL ") + msg)
        if not cond:
            fails.append(msg)

    # --- event alignment
    run = ev["EventHeader/EventHeader.runNumber"].array()
    evt = ev["EventHeader/EventHeader.eventNumber"].array()
    check(bool(ak.all(ak.flatten(run) == src.IEVENTH.run)), "run numbers match")
    check(bool(ak.all(ak.flatten(evt) == src.IEVENTH.event)), "event numbers match")

    # --- counts
    def n_of(branch: str) -> np.ndarray:
        return ak.to_numpy(ak.num(ev[branch].array()))

    for coll, bank, key in (
        ("ChargedTracks", "PHCHRG", "ChargedTracks/ChargedTracks.chi2"),
        ("Clusters", "PHKLUS", "Clusters/Clusters.energy"),
        ("ReconstructedParticles", "PHPSUM", "ReconstructedParticles/ReconstructedParticles.energy"),
    ):
        check(bool(np.all(n_of(key) == ak.to_numpy(ak.num(src[bank])))), f"{coll} count == {bank} count")
    is_mc = "MCPART" in src.fields
    if is_mc:
        check(bool(np.all(n_of("MCParticles/MCParticles.PDG") == ak.to_numpy(ak.num(src.MCPART)))), "MCParticles count == MCPART count")

    # --- tracks: pT from omega, charge sign, phi, tanLambda
    ts = ev["_ChargedTracks_trackStates"].array()
    omega = ak.to_numpy(ak.flatten(ts["_ChargedTracks_trackStates.omega"]))
    phi = ak.to_numpy(ak.flatten(ts["_ChargedTracks_trackStates.phi"]))
    tanl = ak.to_numpy(ak.flatten(ts["_ChargedTracks_trackStates.tanLambda"]))
    h = ak.to_numpy(ak.flatten(src.PHCHRG.hlxpar))
    q = ak.to_numpy(ak.flatten(src.PHCHRG.charge))
    pt_edm = OMEGA_PER_INVPT / np.abs(omega)
    pt_sld = 1.0 / h[:, 1]
    check(np.allclose(pt_edm, pt_sld, rtol=1e-5), f"track pT round-trip (max rel dev {np.max(np.abs(pt_edm / pt_sld - 1)):.1e})")
    check(bool(np.all(np.sign(omega) == np.sign(q))), "sign(omega) == track charge")
    check(np.allclose(phi, h[:, 0]) and np.allclose(tanl, h[:, 2]), "phi, tanLambda round-trip")
    cov = ak.to_numpy(ak.flatten(ev["_ChargedTracks_trackStates/_ChargedTracks_trackStates.covMatrix.values[21]"].array()))
    d = ak.to_numpy(ak.flatten(src.PHCHRG.dhlxpar))
    check(np.allclose(cov[:, 2], d[:, 0]), "cov(phi,phi) round-trip")
    check(np.allclose(cov[:, 5], d[:, 2] * OMEGA_PER_INVPT**2), "cov(omega,omega) round-trip")
    check(np.allclose(cov[:, 0], d[:, 9] * CM_TO_MM**2) and np.allclose(cov[:, 9], d[:, 14] * CM_TO_MM**2), "cov(d0,d0), cov(z0,z0) round-trip")
    print(f"     track sigma: phi {np.sqrt(np.median(cov[:,2]))*1e3:.2f} mrad, d0 {np.sqrt(np.median(cov[:,0]))*1e3:.0f} um, z0 {np.sqrt(np.median(cov[:,9]))*1e3:.0f} um")

    # --- reconstructed particles
    rp = ev["ReconstructedParticles"].array()
    for c in "xyz":
        a = ak.to_numpy(ak.flatten(rp[f"ReconstructedParticles.momentum.{c}"]))
        b = ak.to_numpy(ak.flatten(src.PHPSUM[f"p{c}"]))
        check(np.allclose(a, b), f"RP momentum {c} round-trip")
    e = ak.to_numpy(ak.flatten(rp["ReconstructedParticles.energy"]))
    check(np.allclose(e, ak.to_numpy(ak.flatten(src.PHPSUM.ptot)), rtol=1e-6), "RP energy == ptot")
    check(bool(np.all(ak.to_numpy(ak.flatten(rp["ReconstructedParticles.charge"])) == ak.to_numpy(ak.flatten(src.PHPSUM.charge)))), "RP charge round-trip")
    # event parameters (podio stores them as key/value branches)
    keys = ev["GPDoubleKeys"].array()
    vals = ev["GPDoubleValues"].array()
    ecm_edm = np.array([float(v[list(k).index("ecm")][0]) for k, v in zip(keys, vals)])
    ecm_src = ak.to_numpy(ak.fill_none(ak.firsts(src.PHBM.ecm), np.nan))
    check(np.allclose(ecm_edm, ecm_src, equal_nan=True), "ecm event parameter round-trip")
    n_trk_rel = ak.to_numpy(ak.sum(rp["ReconstructedParticles.tracks_end"] - rp["ReconstructedParticles.tracks_begin"], axis=1))
    n_pt_trk = ak.to_numpy(ak.sum(src.PHPOINT.phchrg_id > 0, axis=1))
    check(bool(np.all(n_trk_rel == n_pt_trk)), "RP->track relations == PHPOINT track pointers")
    n_cl_rel = ak.to_numpy(ak.sum(rp["ReconstructedParticles.clusters_end"] - rp["ReconstructedParticles.clusters_begin"], axis=1))
    n_pt_cl = ak.to_numpy(ak.sum(src.PHPOINT.phklus_id > 0, axis=1))
    check(bool(np.all(n_cl_rel == n_pt_cl)), "RP->cluster relations == PHPOINT cluster pointers")

    # --- clusters
    ce = ak.to_numpy(ak.flatten(ev["Clusters/Clusters.energy"].array()))
    check(np.allclose(ce, ak.to_numpy(ak.flatten(src.PHKLUS.eraw))), "cluster energy round-trip")
    pos = ev["Clusters"].array()
    r = np.hypot(ak.to_numpy(ak.flatten(pos["Clusters.position.x"])), ak.to_numpy(ak.flatten(pos["Clusters.position.y"])))
    z = ak.to_numpy(ak.flatten(pos["Clusters.position.z"]))
    cth = z / np.hypot(r, z)
    check(np.allclose(cth, ak.to_numpy(ak.flatten(src.PHKLUS.cth)), atol=1e-5), "cluster cos(theta) from position")

    if not is_mc:
        return 1 if fails else 0

    # --- MC particles
    mc = ev["MCParticles"].array()
    pdg = ak.to_numpy(ak.flatten(mc["MCParticles.PDG"]))
    gs = ak.to_numpy(ak.flatten(mc["MCParticles.generatorStatus"]))
    mass = ak.to_numpy(ak.flatten(mc["MCParticles.mass"]))
    sim = ak.to_numpy(ak.flatten(mc["MCParticles.simulatorStatus"]))
    print(f"     PDG assigned for {np.mean(pdg != 0) * 100:.2f}% of MC particles; generatorStatus counts: "
          + ", ".join(f"{k}:{v}" for k, v in zip(*np.unique(gs, return_counts=True))))
    check(np.mean(pdg != 0) > 0.999, "PDG assigned for >99.9% of MC particles")
    final = mc["MCParticles.generatorStatus"] == 1
    qsum = ak.to_numpy(ak.sum(mc["MCParticles.charge"][final], axis=1))
    print(f"     net charge of generator-level final state: mean {qsum.mean():+.3f}, |q|>0.5 in {np.mean(np.abs(qsum) > 0.5) * 100:.1f}% of events")
    print(f"     generator-level final-state particles per event: {ak.to_numpy(ak.sum(final, axis=1)).mean():.1f}")
    check(np.mean(np.abs(qsum) > 0.5) < 0.05, "net generator charge ~0 in >95% of events (charge/sign convention)")
    created = (sim >> 30) & 1  # BITCreatedInSimulation
    print(f"     created-in-simulation fraction: {created.mean() * 100:.1f}%  (SWUM bit set: {np.mean((ak.to_numpy(ak.flatten(src.MCPART.origin)) >> 20) & 1) * 100:.1f}%)")
    pion = np.abs(pdg) == 211
    check(np.allclose(mass[pion], 0.13957), "pion mass recovered where the energy word was omitted")
    # parent/daughter consistency: daughter's vertex == parent's endpoint
    pb = mc["MCParticles.parents_begin"]
    pe = mc["MCParticles.parents_end"]
    par_idx = ev["_MCParticles_parents/_MCParticles_parents.index"].array()
    has_par = (pe - pb) == 1
    pidx = par_idx[pb[has_par]]
    vx = mc["MCParticles.vertex.x"][has_par]
    ex = mc["MCParticles.endpoint.x"][pidx]
    check(bool(ak.all(vx == ex)), "daughter vertex == parent endpoint")
    n_par = ak.to_numpy(ak.sum(has_par, axis=1))
    n_src = ak.to_numpy(ak.sum(src.MCPART.parent_id > 0, axis=1))
    print(f"     parents resolved: {n_par.sum()} of {n_src.sum()} non-zero parent_id ({n_par.sum() / n_src.sum() * 100:.2f}%)")
    # links
    ntl = ak.to_numpy(ak.num(ev["ChargedTracksMCLinks/ChargedTracksMCLinks.weight"].array())).sum()
    print(f"     track-MC links: {ntl} for {ak.to_numpy(ak.sum(src.PHCHRG.imc != 0, axis=1)).sum()} tracks with imc != 0")
    nrl = ak.to_numpy(ak.num(ev["ReconstructedParticlesMCLinks/ReconstructedParticlesMCLinks.weight"].array())).sum()
    print(f"     RP-MC links: {nrl} for {ak.to_numpy(ak.num(src.MCPNT)).sum()} MCPNT rows")

    print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURES"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
