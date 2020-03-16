#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from datetime import datetime

import scanpy as sc
import pandas as pd


DATA_DIR = Path(__file__).parent.parent.parent / "data"

PATHS = {
    "global_with_dub": DATA_DIR / "h5ad" / "full-final.h5ad",
    "global_no_dub": DATA_DIR / "h5ad" / "full-no-dub-final.h5ad",
    "neuronal": DATA_DIR / "h5ad" / "neuronal-final.h5ad",
    "nonneuronal": DATA_DIR / "h5ad" / "nonneuronal-final.h5ad",
    "supramammillary": DATA_DIR / "h5ad" / "glut8-subclustering.h5ad",
}
MARKER_PATHS = {
    "global": DATA_DIR / "markers" / "global_markers.csv",
    "neuronal": DATA_DIR / "markers" / "neuronal_markers.csv",
    "nonneuronal": DATA_DIR / "markers" / "nonneuronal_markers.csv",
    "mammillary": DATA_DIR / "markers" / "mammillary_markers.csv",
}

CURRENT_DATE = datetime.now().strftime("%Y%m%d")
FIG_PATH = DATA_DIR.parent / "figures" / CURRENT_DATE


def load_adata(key):
    path = PATHS.get(key, None)
    assert path, f"Key must be in {tuple(PATHS.keys())}"
    assert path.exists()
    return sc.read(path)


def save_figure(fig, pdir, name, dpi=600, ext="pdf"):
    if not FIG_PATH.exists():
        FIG_PATH.mkdir(parents=True)

    out = FIG_PATH / pdir / f"{name}.{ext}"
    pout = FIG_PATH / pdir
    if not pout.exists():
        pout.mkdir()
    fig.savefig(out, dpi=dpi)#, bbox="tight")
    print(f"Saved to {out.name}")


def load_markers(key):
    path = MARKER_PATHS.get(key, None)
    assert path, f"Key must be in {tuple(MARKER_PATHS.keys())}"
    assert path.exists()
    return pd.read_csv(path, index_col=0, header=0)
