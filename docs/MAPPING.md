# SLD mini-DST banks → EDM4hep

Target layout: the collections used by the ALEPH/OPAL legacy conversions, so that the
same reader works for all three (collection names in `code`).

| SLD bank | fields | EDM4hep | notes |
|---|---|---|---|
| IEVENTH | run, event, evttime, weight, trigger, evttype | `EventHeader` (runNumber, eventNumber, timeStamp, weight); trigger and evttype as event parameters | one per event |
| PHBM | ecm, decm, pol, dpol, pos[3], dpos[6] | event parameters (ecm, decm, pol, dpol); `PrimaryVertex` (position = pos, covMatrix = dpos) | ecm is 0 in simulation (unfilled), −999 in data when unmeasured |
| PHPSUM | px, py, pz, ptot, charge, x, y, z, status | `ReconstructedParticles` (momentum, energy, charge, referencePoint = (x,y,z)) | no energy or mass stored: energy = ptot, mass = 0 (massless convention, documented) |
| PHPOINT | phpsum_id → phchrg_id, phklus_id, phkelid_id, phwic_id, phcrid_id | relations `ReconstructedParticles.tracks`, `.clusters`, `.particleIDs` | master pointer bank, one row per PHPSUM entry |
| PHCHRG | hlxpar[6], dhlxpar[15], tkpar[5], tkpar0, charge, chi2, ndfdt, nhit, nhitv, nhite, nhitp, length, status, smwstat, dedx, imc | `ChargedTracks` (`Track`: one `TrackState` at the IP with 5 parameters + 15 covariance terms; chi2, ndf; subdetectorHitNumbers = [nhit CDC, nhitv VXD]) ; dedx → `ChargedTracksDqDx`; imc → track–MCParticle link | **helix convention to be confirmed from the corpus**: which 6 parameters, order, units (cm/GeV?), sign of curvature |
| PHKLUS | eraw, cth, phi, elayer[8], nhit2, nhit3, widths | `Clusters` (energy = eraw, position on the LAC surface from (cosθ, φ), subdetectorEnergies = elayer) | **LAC radius / z-extent to be confirmed** |
| PHKCHRG | phchrg_id, phklus_id, match_distance, delta_phi, delta_theta | not needed for relations (PHPOINT has them); kept as `Cluster.shapeParameters`? — decide | |
| PHKELID | electron ID per track: idstat, prob, E/p-type quantities | `ParticleID` (algorithmType = electron, likelihood = prob, parameters) attached to the RP | |
| PHWIC | muon ID per track: idstat, hit patterns, fit | `ParticleID` (algorithmType = muon) | |
| PHCRID | Cherenkov likelihoods (simulation files only in this release?) | `ParticleID` (algorithmType = CRID) | |
| MCHEAD | ntot, origin, ipx, ipy, ipz | event parameters; MC interaction point | simulation only |
| MCPART | px, py, pz, e, ptype, charge, origin, parent_id, xt | `MCParticles` (PDG from ptype **via a mapping table to be taken from the corpus**, momentum, mass = √(e²−p²), charge, vertex = xt, generatorStatus from origin bits, parents/daughters from parent_id) | simulation only; ~500 per event, parton and hadron levels mixed |
| MCPNT | pointers | to be inspected | simulation only |

## Conventions to confirm (from the SLD documentation corpus)

1. PHCHRG.hlxpar: parameterization, order, units, reference point; dhlxpar packing order.
2. PHCHRG.dedx: appears as a packed integer (e.g. 980193684) — decode.
3. MCPART.ptype: SLD/LUND particle code table (a b quark is 504; 1000094 is an
   intermediate object) → PDG.
4. MCPART.origin bit meanings (0x100000 marks detector-produced secondaries, seen in v3).
5. LAC geometry for cluster positions.
6. PHKLUS.status, PHCHRG.status/smwstat bit meanings.

## Deliberate choices

- Energy of a reconstructed particle is set to |p| (massless); the SLD banks carry no
  mass hypothesis. Users who want a pion-mass convention can recompute from momentum.
- Only reconstruction-level content is converted; the jazelle raw format is not touched.
- Every output file carries, as podio metadata, the input file name, its SHA-256, the
  converter version and the mapping-document version.
