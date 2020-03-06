#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import scanpy as sc


prefix = "/projects/robson-lab/research/hypothalamus/final_analysis/paper"
full_raw = sc.read(f"{prefix}/full-raw.h5ad")
full = sc.read(f"{prefix}/full-celltyped.h5ad")
neuronal = sc.read(f"{prefix}/neuronal.h5ad")
nonneuronal = sc.read(f"{prefix}/nonneuronal.h5ad")


neuronal_final = sc.read(f"{prefix}/neuronal-final.h5ad")
nonneuronal_final = sc.read(f"{prefix}/nonneuronal-final.h5ad")



raw = sc.read(f"{prefix}/raw.h5ad")


# In[23]:


full_no_dub = full[~full.obs.celltype_broad.isin(["doublet", "putative doublet"])].copy()


subclust_data = {
    "cluster_13_17": sc.read(
        "/projects/robson-lab/research/hypothalamus/final_analysis/"
        "subcluster_hdc_like-ver2/aggr_hdc_like_neurons-clustered_20191212.h5ad"
    ),
    #"cluster_1-7": sc.read(
    #    "/projects/robson-lab/research/hypothalamus/final_analysis/"
    #    "subcluster_/"
    #),
    #"cluster_1-6": sc.read(
    #    "/projects/robson-lab/research/hypothalamus/final_analysis/"
    #    "subcluster_"
    #),
    "cluster_7": sc.read(
        "/projects/robson-lab/research/hypothalamus/final_analysis/"
        "subcluster_premammillary/aggr_premammillary_neurons-clustered_20191212.h5ad"
    ),
    "cluster_8": sc.read(
        "/projects/robson-lab/research/hypothalamus/final_analysis/"
        "subcluster_supramammillary/aggr_supramammillary_neurons-clustered_20191212.h5ad"
    ),    
}

