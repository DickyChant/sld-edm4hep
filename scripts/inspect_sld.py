#!/usr/bin/env python
"""Survey one SLD mini-DST: banks, fields, per-event multiplicities, value ranges.

Pure pyarrow/awkward; runs without podio. Used to fix the mapping in docs/MAPPING.md.

    python scripts/inspect_sld.py FILE.parquet [--events 2000]
"""

from __future__ import annotations

import argparse

import awkward as ak
import numpy as np
import pyarrow.parquet as pq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--events", type=int, default=2000)
    args = ap.parse_args()

    pf = pq.ParquetFile(args.path)
    print(f"{args.path}\n  rows {pf.metadata.num_rows}, row groups {pf.metadata.num_row_groups}")
    a = ak.from_parquet(args.path, row_groups=[0])[: args.events]
    for bank in a.fields:
        col = a[bank]
        n = ak.num(col, axis=1)
        print(f"\n== {bank}: per event min/mean/max = {ak.min(n)}/{ak.mean(n):.2f}/{ak.max(n)}")
        flat = ak.flatten(col)
        for f in flat.fields:
            v = flat[f]
            try:
                if v.ndim > 1:
                    lens = ak.num(v, axis=1)
                    v2 = ak.flatten(v)
                    print(f"   {f:14s} list len {ak.min(lens)}..{ak.max(lens)}   "
                          f"range [{ak.min(v2):.4g}, {ak.max(v2):.4g}]")
                else:
                    vv = ak.to_numpy(v)
                    uniq = np.unique(vv)
                    extra = f"  {len(uniq)} distinct" + (f": {uniq[:12]}" if len(uniq) <= 12 else "")
                    print(f"   {f:14s} range [{vv.min():.6g}, {vv.max():.6g}]{extra}")
            except Exception as exc:  # noqa: BLE001
                print(f"   {f:14s} ({type(v.type).__name__}) {exc}")


if __name__ == "__main__":
    main()
