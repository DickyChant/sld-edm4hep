#!/usr/bin/env bash
# Run a command inside the Key4hep container with CVMFS mounted.
#
#   scripts/run_key4hep.sh python -m sld_edm4hep.convert IN.parquet OUT.root
#   scripts/run_key4hep.sh bash            # interactive shell with the stack sourced
#
# Requires: docker, /cvmfs/sw.hsf.org mounted on the host (cvmfs client).
# Heavy input/output lives under $SLD_DATA (default /media/sqian/WD_RED), mounted at
# the same path inside the container so that paths in logs stay valid on the host.
set -euo pipefail

IMAGE=${KEY4HEP_IMAGE:-ghcr.io/key4hep/key4hep-images/ubuntu24:latest}
SLD_DATA=${SLD_DATA:-/media/sqian/WD_RED}
REPO=$(cd "$(dirname "$0")/.." && pwd)
SETUP=${KEY4HEP_SETUP:-/cvmfs/sw.hsf.org/key4hep/setup.sh}

if [ ! -d /cvmfs/sw.hsf.org ]; then
  echo "error: /cvmfs/sw.hsf.org is not mounted on the host" >&2
  exit 1
fi

TTY_FLAGS=""
if [ -t 0 ] && [ -t 1 ]; then TTY_FLAGS="-it"; fi

# shellcheck disable=SC2086
exec docker run --rm $TTY_FLAGS \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp/k4home -e PYTHONPATH=/work -e PYTHONDONTWRITEBYTECODE=1 \
  -v /cvmfs:/cvmfs:ro,shared \
  -v "$SLD_DATA:$SLD_DATA" \
  -v "$REPO:/work" -w /work \
  "$IMAGE" bash -c 'mkdir -p /tmp/k4home; source "'"$SETUP"'" > /tmp/k4home/setup.log 2>&1 || { cat /tmp/k4home/setup.log; exit 1; }; exec "$@"' -- "$@"
