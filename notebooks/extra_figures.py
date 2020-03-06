#!/usr/bin/env python
# coding: utf-8
from scanpy_recipes import sc
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from plotting.plotting_funcs import *
from utils import *


def plot_doublets():
    from matplotlib.patheffects import withStroke
    fig, ax = plt.subplots(figsize=(4, 4))
    palette = sns.palettes.mpl_palette("tab20c", 16) + sns.palettes.mpl_palette("tab20b", 16)
    sc.pl.umap(full, color="cluster_revised", show=False, title="", ax=ax, palette=palette)#, legend_loc="on data")
    ax.set_xlabel(""); ax.set_ylabel("")
    #for child in ax.get_children():
    #    if isinstance(child, matplotlib.text.Text):
    #        child.set_path_effects([withStroke(linewidth=1, foreground='white')])
    legend = ax.legend(ncol=16, bbox_to_anchor=(0.5, 0.), fontsize=6, frameon=False, 
                       title="Unfiltered Cluster", columnspacing=0.5, loc="upper center")
    legend.get_title().set_fontsize(8)
    fix_aspect_scatter_with_legend(fig)
    save_figure(fig, "SI_UMAP-all-clusters")



    doublet_clusters = ["21", "24", "25", "26", "27", "28", "30", "31"]
    doublet_clusters_plus = ["2", "10", "14", "17", "19", "20", "21", "24", "25", "26", "27", "28", "30", "31"]
    var_group_labels = ["neuron", "astro", "OPC", "oligo", "endothelial", "pericyte"]
    var_group_positions = list(zip(range(0, 18, 3), range(2, 18, 3)))
    doublet_genes = [
        "Tubb3", "Syp", "Elavl2",
        "Aqp4", "Agt", "Slc1a3",
        "Pdgfra", "C1ql1", "Gpr17",
        "Mag", "Mog", "Cldn11",
        "Cd34", "Slc2a1", "Itm2a",
        "Rgs5", "Pdgfrb", "Vtn"
    ]


    from matplotlib.patches import Rectangle


    subset = full[full.obs.cluster_revised.isin(doublet_clusters)].copy()
    #subset.obs.cluster_revised.cat.reorder_categories(["24", "25", "27", "28", "31", "26", "29", "30"], inplace=True)
    cg = sc.pl.matrixplot(
        subset, 
        doublet_genes, use_raw=False, groupby="cluster_revised", 
        var_group_labels=var_group_labels, var_group_positions=var_group_positions,
        var_group_rotation=0, 
        cmap="Reds", vmin=0, vmax=3,
        show=False, figsize=(5, 2)
    )
    fig = plt.gcf()
    axs = fig.get_axes()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha="right", size=8, va="top", rotation_mode="anchor");
    axs[0].tick_params(axis="x", pad=1)
    axs[0].set_ylabel("Unfiltered Cluster")
    for y in range(3, 18, 3):
        axs[0].axvline(y, lw=1, color="k")
    axs[2].set_ylabel("Log Normalized UMIs")
    for child in axs[1].get_children():
        if isinstance(child, matplotlib.text.Text):
            child.set_fontsize(8)
        elif isinstance(child, matplotlib.patches.PathPatch):
            child.set_linewidth(0.75)
    fig.subplots_adjust(right=0.9, left=0.1, bottom=0.2, top=0.9)
    fig.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
    save_figure(fig, "SI_doublet_matrix_ver1")


# In[391]:


    subset = full[full.obs.cluster_revised.isin(doublet_clusters)].copy()
    subset.obs.cluster_revised.cat.reorder_categories(["24", "25", "21", "27", "28", "31", "26", "30"], inplace=True)
    cg = sc.pl.matrixplot(
        subset,
        doublet_genes, use_raw=False, groupby="cluster_revised", 
        var_group_labels=var_group_labels, var_group_positions=var_group_positions,
        var_group_rotation=0, 
        cmap="Reds", vmin=0, vmax=3,
        show=False, figsize=(5, 2)
    )
    fig = plt.gcf()
    axs = fig.get_axes()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha="right", size=8, va="top", rotation_mode="anchor");
    axs[0].tick_params(axis="x", pad=1)
    axs[0].set_ylabel("Unfiltered Cluster")
    rect_params = dict(width=3, height=1, facecolor="none", edgecolor="black", linewidth=1)
    for xy in [(0,0), (3,0), (0,1), (3,1), (0,2), (6,2),
               (0, 3), (9, 3), (0, 4), (9, 4), 
               (0, 5), (12, 5), 
               (3, 6), (9, 6),
               (12, 7), (15, 7)
              ]:
        axs[0].add_patch(Rectangle(xy, **rect_params))
    for child in axs[1].get_children():
        if isinstance(child, matplotlib.text.Text):
            child.set_fontsize(8)
        elif isinstance(child, matplotlib.patches.PathPatch):
            child.set_linewidth(0.75)
    axs[2].set_ylabel("Log Normalized UMIs")
    fig.subplots_adjust(right=0.9, left=0.1, bottom=0.2, top=0.9)
    fig.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
    save_figure(fig, "SI_doublet_matrix_ver2")


    # In[392]:


    subset = full[full.obs.cluster_revised.isin(doublet_clusters_plus)].copy()
    subset.obs.cluster_revised.cat.reorder_categories(
        ["2", "10", "17", "19", "20", "14", "24", "25", "21", "27", "28", "31", "26", "30"], inplace=True)
    cg = sc.pl.matrixplot(
        subset,
        doublet_genes, use_raw=False, groupby="cluster_revised", 
        var_group_labels=var_group_labels, var_group_positions=var_group_positions,
        var_group_rotation=0, 
        cmap="Reds", vmin=0, vmax=3,
        show=False, figsize=(5, 2)
    )
    fig = plt.gcf()
    axs = fig.get_axes()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha="right", size=8, va="top", rotation_mode="anchor");
    axs[0].tick_params(axis="x", pad=1)
    axs[0].set_ylabel("Unfiltered Cluster")
    rect_params = dict(width=3, height=1, facecolor="none", edgecolor="black", linewidth=1)
    for child in axs[1].get_children():
        if isinstance(child, matplotlib.text.Text):
            child.set_fontsize(8)
        elif isinstance(child, matplotlib.patches.PathPatch):
            child.set_linewidth(0.75)
    axs[2].set_ylabel("Log Normalized UMIs")
    fig.subplots_adjust(right=0.9, left=0.1, bottom=0.2, top=0.9)
    fig.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
    save_figure(fig, "SI_doublet_matrix_ver3")


    subset = full[full.obs.cluster_revised.isin(doublet_clusters)].copy()
    #subset.obs.cluster_revised.cat.reorder_categories(["24", "25", "27", "28", "31", "26", "29", "30"], inplace=True)
    cg = sc.pl.matrixplot(
        subset, 
        doublet_genes, use_raw=False, groupby="cluster_revised", 
        var_group_labels=var_group_labels, var_group_positions=var_group_positions,
        var_group_rotation=0, 
        cmap="Reds", vmin=0, vmax=3,
        show=False, figsize=(5, 2)
    )
    fig = plt.gcf()
    axs = fig.get_axes()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha="right", size=8, va="top", rotation_mode="anchor");
    axs[0].tick_params(axis="x", pad=1)
    axs[0].set_ylabel("Unfiltered Cluster")
    axs[2].set_ylabel("Log Normalized UMIs")
    for child in axs[1].get_children():
        if isinstance(child, matplotlib.text.Text):
            child.set_fontsize(8)
        elif isinstance(child, matplotlib.patches.PathPatch):
            child.set_linewidth(0.75)
    fig.subplots_adjust(right=0.9, left=0.1, bottom=0.2, top=0.9)
    fig.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
    save_figure(fig, "SI_doublet_matrix_ver4")


    subset = full[full.obs.cluster_revised.isin(["10", "29", "20"])].copy()
    subset.obs.cluster_revised.cat.reorder_categories(["10", "29", "20"], inplace=True)
    cg = sc.pl.matrixplot(
        subset, 
        ["Agt", "Slc1a3", 'Ntrk2', 'Plpp3', 'Atp1b2', 
        'Ttyh1', 
        'Serpine2', 'Atp1a2', 'Slc6a11',
         'Gpr37l1',
         'Nnat', 'Slc4a4',
         "Cd34", "Slc2a1", "Kdr", "Pecam1",
        ], 
        use_raw=False, groupby="cluster_revised", 
        var_group_rotation=0, 
        cmap="Reds", vmin=0, vmax=3,
        show=False, figsize=(5, 2)
    )
    fig = plt.gcf()
    axs = fig.get_axes()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha="right", size=8, va="top", rotation_mode="anchor");
    axs[0].tick_params(axis="x", pad=1)
    axs[0].set_ylabel("Unfiltered Cluster")
    axs[1].set_ylabel("Log Normalized UMIs")
    for child in axs[1].get_children():
        if isinstance(child, matplotlib.text.Text):
            child.set_fontsize(8)
    fig.subplots_adjust(right=0.9, left=0.1, bottom=0.2, top=0.9)
    fig.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
    save_figure(fig, "SI_doublet_matrix_cluster-29")

    marker_violins(
        full[full.obs.cluster_revised.isin(list("123456789") + ["22"])], 
        ["Slc17a6", "Slc32a1", "Gad1", "Gad2", "Nts", "Onecut2", "Tac2", "Synpr", "Ar"], 
        "cluster_revised", palette[:9] + [palette[21]], cluster_name="Unfiltered Cluster",
        filename="SI_doublet_violins_cluster-22"
    )

    fig, ax = plt.subplots(figsize=(4, 2))
    cluster23_violin_data = full[full.obs.cluster_revised.isin(list("123456789") + ["21", "22", "23"])].obs[
        ["n_genes_by_counts", "total_counts", "cluster_revised"]
    ]
    sns.violinplot(
        data=cluster23_violin_data, x="cluster_revised", y="n_genes_by_counts", ax=ax, 
        inner=None, linewidth=1, scale="width", palette=palette[:9] + palette[20:23]
    )
    ax.axhline(2200, ls="--", color="k", lw=1)
    ax.set_ylabel("Genes")
    ax.set_xlabel("Unfiltered Cluster")
    sns.despine(fig)
    fig.tight_layout()
    save_figure(fig, "SI_doublet_violins_cluster-23")

