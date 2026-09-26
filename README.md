# sld-edm4hep

Conversion of the SLD legacy mini-DSTs (SLC, Z pole, 1996–1998 runs and simulation;
Zenodo record [22732257](https://zenodo.org/records/22732257), Parquet form as exported by
[jazelle_reader](https://github.com/AlkaidCheng/jazelle_reader)) to
[EDM4hep](https://github.com/key4hep/EDM4hep), in the same layout as the ALEPH and OPAL
legacy-data conversions (podio 1.7, EDM4hep 1.0), so that one set of tools reads all
three experiments.

Status: converter and read-back validation complete for every bank in the release;
see [docs/MAPPING.md](docs/MAPPING.md) for the bank-by-bank mapping, the decoded
conventions and the remaining ambiguities.

## Layout

| path | content |
|---|---|
| `sld_edm4hep/convert.py` | the converter: reads Parquet with pyarrow/awkward, writes with podio |
| `sld_edm4hep/ptypes.py` | SLD particle-type code → PDG table (JETSET 6.3 numbering + SLD extensions) |
| `scripts/run_key4hep.sh` | runs any command inside the Key4hep container with CVMFS mounted |
| `scripts/convert_all.sh` | converts and validates every file, several containers in parallel |
| `scripts/validate.py` | read-back checks against the source Parquet |
| `scripts/inspect_sld.py` | schema and content survey of a mini-DST, no podio needed |
| `docs/MAPPING.md` | mapping, conventions, provenance, open questions |

## Environment

The Key4hep software stack is taken from CVMFS (`/cvmfs/sw.hsf.org/key4hep`) inside the
official container `ghcr.io/key4hep/key4hep-images/ubuntu24`; nothing is built locally.

```
scripts/run_key4hep.sh python -m sld_edm4hep.convert IN.parquet OUT.edm4hep.root [--max-events N]
scripts/run_key4hep.sh python scripts/validate.py OUT.edm4hep.root IN.parquet
scripts/convert_all.sh            # all data and simulation files
```

`run_key4hep.sh` mounts `/cvmfs` (read-only, shared), the data disk (`$SLD_DATA`,
default `/media/sqian/WD_RED`) at the same path, and this repository at `/work`. A
30,000-event simulation file converts in about 18 minutes; a 10,000-event data file in
about 70 seconds.

## Output

One `events` frame per SLD event with the collections

`EventHeader`, `PrimaryVertex`, `ReconstructedParticles`, `ChargedTracks`, `Clusters`,
`ElectronID`, `MuonID`, `CRID` (ParticleID), and for simulation `MCParticles`,
`ChargedTracksMCLinks`, `ReconstructedParticlesMCLinks`;

event parameters `ecm, decm, pol, dpol, trigger, evttype, header_word` (+ `mc_ntot,
mc_origin, mc_primary_quark_pdg, mc_ip_mm`); raw SLD words that have no EDM4hep field in
`UserDataCollection`s (`ChargedTracksDedxRaw`, `ChargedTracksImcRaw`,
`ChargedTracksImpact{2D,3D}{,Significance}`, `MCParticlesSLD{Type,Origin,Id}`); and one
`metadata` frame with the input file name and SHA-256, converter and mapping versions,
units, and the parameter orderings.

Lengths are in mm, energies in GeV, following EDM4hep. The reconstructed-particle energy
is |p| (the SLD banks carry no mass hypothesis).
