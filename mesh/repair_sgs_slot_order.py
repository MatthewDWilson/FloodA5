"""
repair_sgs_slot_order.py
------------------------
Repair the SGS edge tables of a mesh built before the WP3 fix (2026-10).

`build_sgs_tables!` used to store each cell's edge sills and R-A edge curves
(`sgs_edge_sills`, `sgs_edge_area_curve`, `sgs_edge_perim_curve`) in slots
1..5 following `A5Grid.grid_disk_neighbours` order — the order of pya5's
`uncompact(grid_disk(cell, 1), res)`, including cells outside the mesh — while
FloodA5 reads slot k as the k-th entry of the cell's stored `neighbours` list.
The two orders differ for ~95% of cells, so most edges used another
neighbour's sill and R-A curves.

This tool recomputes the build-time order with pya5 0.x (it must be the same
grid_disk the tables were built with) from the ORIGINAL mesh (pya5 0.x cell
IDs, e.g. `*_a5v0.parquet` kept by mesh/migrate_mesh_ids.py), reorders the
three edge arrays into neighbour-list order (any valence), adds
`sgs_layout = 2`, and writes them into TARGET — the same mesh, either the
original itself or its pya5 1.x migration (same row order; checked by cell
centre). Nothing else changes, and no DEM is needed. The result equals a
rebuild with the fixed `build_sgs_tables!`.

Requires pya5 0.x. Usage:
    python mesh/repair_sgs_slot_order.py ORIGINAL_a5v0.parquet TARGET.parquet OUT.parquet
    python mesh/repair_sgs_slot_order.py ORIGINAL_a5v0.parquet TARGET.parquet --in-place
Prints one JSON line; exits non-zero on failure.
"""

import json
import os
import sys


def _fail(msg):
    print(json.dumps({"error": msg}), flush=True)
    sys.exit(1)


def repair(orig_path, target_path, out_path):
    import a5
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    import importlib.metadata as md

    v = md.version("pya5")
    if int(v.split(".")[0]) >= 1:
        raise RuntimeError(f"pya5 {v} installed; this repair needs the pya5 0.x grid_disk "
                           "order the tables were built with (e.g. a venv with pya5==0.7.3)")
    orig = pq.read_table(orig_path)
    tgt = pq.read_table(target_path)
    if orig.num_rows != tgt.num_rows:
        raise RuntimeError("ORIGINAL and TARGET have different numbers of cells")
    for col in ("center_lon", "center_lat"):
        a = np.asarray(orig.column(col).to_pylist()); b = np.asarray(tgt.column(col).to_pylist())
        if np.max(np.abs(a - b)) > 1e-9:
            raise RuntimeError(f"ORIGINAL and TARGET rows do not match ({col})")
    names = orig.column_names
    for col in ("sgs_edge_sills", "sgs_edge_area_curve", "sgs_edge_perim_curve", "sgs_elev_bins"):
        if col not in names:
            raise RuntimeError(f"{orig_path} has no {col} column")
    if "sgs_layout" in names:
        raise RuntimeError(f"{orig_path} already has the new SGS layout")

    ids = orig.column("cell_id").to_pylist()
    nbrs = orig.column("neighbours").to_pylist()
    n_bins = len(orig.column("sgs_elev_bins")[0].as_py())
    sills = orig.column("sgs_edge_sills").to_pylist()
    area = orig.column("sgs_edge_area_curve").to_pylist()
    perim = orig.column("sgs_edge_perim_curve").to_pylist()
    old_slots = len(sills[0])
    n_slots = max(1, max(len(nb or []) for nb in nbrs))

    new_sills, new_area, new_perim = [], [], []
    n_moved = n_lost = 0
    for h, nb, s_old, a_old, p_old in zip(ids, nbrs, sills, area, perim):
        c = int(h, 16)
        res = a5.get_resolution(c)
        gd = [f"{int(x) & 0xFFFFFFFFFFFFFFFF:016x}" for x in a5.uncompact(a5.grid_disk(c, 1), res)]
        gd = [x for x in gd if x != f"{c & 0xFFFFFFFFFFFFFFFF:016x}"]
        src = {x: s for s, x in enumerate(gd[:old_slots])}        # build-time slot of each id
        a_old = np.asarray(a_old, dtype=float).reshape(old_slots, n_bins)   # column-major: slot rows
        p_old = np.asarray(p_old, dtype=float).reshape(old_slots, n_bins)
        s_new = [float("nan")] * n_slots
        a_new = np.zeros((n_slots, n_bins)); p_new = np.zeros((n_slots, n_bins))
        for k, nbk in enumerate(nb or []):
            s = src.get(f"{int(nbk, 16):016x}")
            if s is None:
                n_lost += 1
                continue
            n_moved += (s != k)
            s_new[k] = s_old[s] if s_old[s] is not None else float("nan")
            a_new[k] = a_old[s]; p_new[k] = p_old[s]
        new_sills.append(s_new)
        new_area.append(a_new.reshape(-1).tolist())
        new_perim.append(p_new.reshape(-1).tolist())

    lf = pa.list_(pa.float64())
    cols = {"sgs_edge_sills": pa.array(new_sills, type=lf),
            "sgs_edge_area_curve": pa.array(new_area, type=lf),
            "sgs_edge_perim_curve": pa.array(new_perim, type=lf)}
    out = tgt
    for name, arr in cols.items():
        i = out.column_names.index(name)
        out = out.set_column(i, pa.field(name, lf), arr)
    out = out.append_column("sgs_layout", pa.array([2.0] * out.num_rows, type=pa.float64()))
    meta = dict(tgt.schema.metadata or {})
    info = json.loads(meta.get(b"flooda5", b"{}").decode() or "{}")
    info["sgs_slot_repair"] = "grid_disk order -> neighbour-list order (repair_sgs_slot_order.py)"
    meta[b"flooda5"] = json.dumps(info).encode()
    tmp = out_path + ".tmp"
    pq.write_table(out.replace_schema_metadata(meta), tmp)
    os.replace(tmp, out_path)
    return {"status": "ok", "original": orig_path, "output": out_path, "cells": len(ids),
            "n_slots": n_slots, "slots_reordered": n_moved, "neighbours_without_curve": n_lost,
            "pya5_version": v}


def main():
    a = sys.argv[1:]
    if len(a) == 3 and a[2] == "--in-place":
        a = [a[0], a[1], a[1]]
    if len(a) != 3:
        _fail("Usage: repair_sgs_slot_order.py ORIGINAL_a5v0.parquet TARGET.parquet (OUT.parquet | --in-place)")
    try:
        print(json.dumps(repair(*a)), flush=True)
    except Exception as e:
        _fail(str(e))


if __name__ == "__main__":
    main()
