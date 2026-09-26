# sld-edm4hep

Conversion of the SLD legacy mini-DSTs (SLC, Z pole, 1996–1998 runs and simulation;
Zenodo record [22732257](https://zenodo.org/records/22732257), Parquet form) to
[EDM4hep](https://github.com/key4hep/EDM4hep), in the same layout as the ALEPH and OPAL
legacy-data conversions (podio 1.7, EDM4hep 1.0), so that one set of tools reads all
three experiments.

Status: in development. See [docs/MAPPING.md](docs/MAPPING.md) for the bank-by-bank
mapping and the open questions.

## Layout

| path | content |
|---|---|
| `sld_edm4hep/` | the converter (reads Parquet with pyarrow/awkward, writes with podio) |
| `scripts/inspect_sld.py` | schema and content survey of a mini-DST, no podio needed |
| `scripts/convert.py` | command-line conversion of one or many files |
| `scripts/validate.py` | read-back checks: counts, four-momentum sums, cross-references |
| `docs/` | mapping, conventions, provenance |

## Environment

podio and EDM4hep are not on conda-forge; they are built from source (tags `v01-07`
and `v01-00`) against a conda-forge ROOT. `scripts/setup_env.sh` creates the
environment on a data disk and builds both. Everything heavy lives under
`/media/sqian/WD_RED/sld-edm4hep`; this repository holds source and documentation only.
