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


def figure_S1b():
    neuronal_markers = ["Snap25", "Tubb3", "Elavl2", "Syp"]
    full_no_dub.obs["neuronal_average_exp"] = full_no_dub[:, neuronal_markers].X.mean(axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        full_no_dub, 
        color="neuronal_average_exp",  
        palette=neuronal_nonneuronal_palette,
        title="Neuronal marker mean expression",
        show=False,
        color_map=red_colormap,
        ax=ax
    )
    sns.despine(fig)
    fig.get_children()[-1].set_ylabel("Log-normalized expression")
    fix_aspect_scatter_with_cbar(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_01", "fig1_neuronal-expression")


def figure_S1c():
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        full_no_dub, 
        color="classification_cluster",  
        palette=neuronal_nonneuronal_palette,
        title="Neuronal Classification",
        show=False,
        ax=ax
    )
    ax.legend(bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, loc="lower center")
    sns.despine(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_01", "fig1_neuronal-classification")


def figure_S1a():
    full_comb = sc.read("analysis/aggr_output_ver1/posterior_hypothalamus_aggr-aggr_20190710.h5ad")
    full_comb = full_comb[full_no_dub.obs_names, :]
    full_comb = sc.pp.preprocess(full_comb, n_top_genes=2500, scale=False)
    sc.pp.dimensionality_reduction(full_comb, n_comps=25, n_neighbors=15, min_dist=0.5)

    full_comb.obs["chemistry"] = full_comb.obs.sampleid.map({
        "AJ18003": "v2",
        "AJ18004": "v2",
        "AJ19001": "v3",
        "AJ19002": "v3",
    })

    full_comb.obsm["X_umap"] = full_comb.obsm["X_umap"].dot(-np.eye(2))

    fig, axarr = plt.subplots(2, 3, figsize=(6, 4))
    params = dict(show=False, size=2)
    sc.pl.umap(
        full_comb, color="chemistry", ax=axarr[0, 0], palette=sns.xkcd_palette(["orange", "grass green"]),
        legend_loc=None, title="10X Chemistry", **params
    )
    sc.pl.umap(
        full_no_dub, color="chemistry", ax=axarr[1, 0], palette=sns.xkcd_palette(["orange", "grass green"]),
        title="", **params
    )
    sc.pl.umap(
        full_comb, color="sex", ax=axarr[0, 1], palette=sns.xkcd_palette(["electric blue", "bright red"]),
        legend_loc=None, title="Sex", **params
    )
    sc.pl.umap(
        full_no_dub, color="sex", ax=axarr[1, 1], palette=sns.xkcd_palette(["electric blue", "bright red"]),
        title="", **params
    )
    sc.pl.umap(
        full_comb, color="sample_name", ax=axarr[0, 2], #palette=sns.xkcd_palette(["orange", "grass green"]),
        legend_loc=None, title="Sample", **params
    )
    sc.pl.umap(
        full_no_dub, color="sample_name", ax=axarr[1, 2], #palette=sns.xkcd_palette(["orange", "grass green"]),
        title="", **params
    )

    legend_params = dict(
        bbox_to_anchor=(0.5, -0.0), loc="upper center", ncol=2, frameon=False, columnspacing=0.9, fontsize=7,
        markerscale=0.5
    )
    for k in range(3):
        leg = axarr[1, k].legend(**legend_params)

    for ax in axarr.flat:
        ax.set_xlabel("")
        ax.set_ylabel("")
    axarr[0,0].set_ylabel("No batch correction", size=8)
    axarr[1,0].set_ylabel("With batch correction", size=8)
    axarr[0,0].set_title("10X Chemistry", size=9)
    axarr[0,1].set_title("Sex", size=9)
    axarr[0,2].set_title("Sample", size=9)
    fix_aspect_scatter_with_legend(fig)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_SI_01", "SI_batch-correction-comparison")


def figure_S3a():
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        nonneuronal_final, color="cluster_final", show=False, ax=ax, 
        title="Nonneuronal clusters", palette=nonneuronal_palette
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_02", "fig2_nonneuronal-clusters_v2")



def figure_S3b():
    nonneuronal_violin_genes = [
        "Ascl1", "Pdgfra", "Fyn", "Mag", "Mobp", "Mal",
        "Tspan2", "Ptgds", "Hcn2", "Fth1",
        "Agt", "Aqp4", "Gfap", 
        "Rax", "Ccdc153", 
        "Mrc1", "Tmem119", "Cx3cr1", 
        "Rgs5", "Acta2", "Fxyd5", "Slc47a1", "Col1a1", "Dcn", "Igfbp2",
        "Pecam1", "Slc38a5"
    ]

    marker_violins(
        nonneuronal_final, nonneuronal_violin_genes, "cluster_final", nonneuronal_palette,
        filename="fig2_nonneuronal-marker-violins", figdir="figure_02"
    )

def figure_S3c():
    nonneuronal_markers = sc.tl.find_marker_genes(nonneuronal_final, cluster_key="cluster_final", log_fold_change=0.75)

    nonneuronal_heatmap_data, nonneuronal_col_colors, nonneuronal_row_colors = get_heatmap_data(
        nonneuronal_final, nonneuronal_markers, "cluster_final", nonneuronal_palette
    )

    cg = sns.clustermap(
        nonneuronal_heatmap_data, z_score=0, vmin=-3, vmax=3, cmap=heatmap_cmap,
        xticklabels=False, yticklabels=False,
        row_cluster=False, col_cluster=False,
        col_colors=nonneuronal_col_colors, 
        row_colors=nonneuronal_row_colors,
        robust=True, figsize=(4, 4), cbar_kws=dict(use_gridspec=True)
    )
    cg.ax_row_colors.set_ylabel("Genes", size="small")
    cg.ax_col_colors.set_title("Cells", size="small", zorder=100)
    cg.ax_col_colors.xaxis.tick_top()

    fix_heatmap_colorbar(cg, nonneuronal_final, nonneuronal_heatmap_data, heatmap_cmap)

    cg.fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_heatmap_annotations(cg)
    save_figure(cg.fig, "figure_02", "fig2_nonneuronal-heatmap", dpi=600, ext="png")

    fig, ax = plt.subplots(figsize=(1, 1.4))
    for k, (cluster, inds) in enumerate(nonneuronal_final.obs.groupby("cluster_final")):
        x, y = nonneuronal_final[inds.index, :].obsm["X_umap"].T
        ax.scatter(x, y, c=nonneuronal_palette[k], marker="s", label=cluster)
    leg = fig.legend(
        *ax.get_legend_handles_labels(), 
        frameon=False, 
        markerscale=1.1, 
        scatteryoffsets=[0.5], 
        fontsize="xx-small", 
        ncol=2,
        handletextpad=0.5
    )
    for t in leg.get_texts():
        t.set_va("baseline")
    ax.clear()
    ax.set_axis_off()
    fig.tight_layout()
    save_figure(fig, "nonneuronal-heatmap-legend")


def figure_S3d():
    genes_umi_violins(
        nonneuronal_final, "cluster_final", nonneuronal_palette, 
        figdir="figure_02", filename="fig2_nonneuronal-umi-genes-violins"
    )

def figure_S6b():
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    sc.pl.umap(
        subclust_data["cluster_8"], 
        color="cluster_revised", show=False, ax=ax, title="", size=16
    )
    leg = ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False, fontsize="xx-small")
    for handle in leg.legendHandles:
        handle._sizes = [16]
    fix_aspect_scatter_with_legend(fig)
    #add_scatter_borders(ax)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_09", "fig9_cluster8_subclusters")



def figure_S6a():
    #hdc_slc32a1_gpr101
    genes = ["Gad1", "Gad2", "Slc32a1", "Slc17a6"]
    scatter_2by2(
        subclust_data["cluster_8"], genes, 
        filename="fig9_cluster8_key-genes-1", figdir="figure_09",
        add_marker_edges=False, label_va="top", label_ha="left", markersize=16
    )

def figure_S6c():
    genes = ["Sema3c", "Inhba", "Rxfp1", "Cpne7"]
    scatter_2by2(
        subclust_data["cluster_8"], genes, 
        filename="fig9_cluster8_key-genes-2", figdir="figure_09",
        add_marker_edges=False, label_va="top", label_ha="left", markersize=16
    )


def figure_S2a():
    neuronal_cluster_props = get_cluster_proportions(neuronal_final, "cluster_revised", sample_key="sample_name")
    fig1 = plot_celltype_proportions(neuronal_cluster_props, neuronal_palette)
    save_figure(fig1, "figure_S2", "SI_neuronal_cluster_props-v2")

def figure_S2b():
    neuronal_sample_props = get_sample_proportions(neuronal_final, "cluster_revised", "sample_name")
    fig2 = plot_celltype_proportions(neuronal_sample_props, sns.mpl_palette("tab20", 4)[::-1])
    save_figure(fig2, "figure_S2", "SI_neuronal_sample_props-v2")


def figure_S2c():
    nonneuronal_cluster_props = get_cluster_proportions(nonneuronal_final, "cluster_final", sample_key="sample_name")
    fig1 = plot_celltype_proportions(nonneuronal_cluster_props, nonneuronal_palette)
    save_figure(fig1, "figure_S2", "SI_nonneuronal_cluster_props-v2")
def figure_S2d():
    fig2 = plot_celltype_proportions(nonneuronal_sample_props, sns.mpl_palette("tab20", 4)[::-1])
    nonneuronal_sample_props = get_sample_proportions(nonneuronal_final, "cluster_final", "sample_name")
    save_figure(fig2, "figure_S2", "SI_nonneuronal_sample_props-v2")

def figure_S5a():
    genes = ["Tac2", "Pdyn", "Esr1", "Prlr"]
    scatter_2by2(neuronal_final, genes, filename="fig4_arcuate-umap-markers", figdir="figure_04", add_marker_edges=False)

def figure_S5b():
    arcuate_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6",
        "Prlr", "Tac2", "Tmem35a", "Esr1", "Nhlh2", "Pdyn", "Mrap2", "Inhbb",
        "Nr5a2", "Rxfp1", "Ucp2", "Lxn", "Trp53i11", "Fndc9", "Col2a1",
        "Kiss1", "Kcnk2"
    ])

    marker_violins(
        neuronal_final, arcuate_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig4_arcuate-marker-violins", figdir="figure_04"
    )

def figure_S7a():
    genes = ["Cck", "Synpr", "Tac1", "Nos1"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig4_premammillary-umap-markers", figdir="figure_04", add_marker_edges=False
    )
def figure_S7b():
    premammillary_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6",
        "Foxb1", "Cck",
        "Dlk1", "Synpr", "Spock1",
        "Nxph1", "Negr1",
        "Ar", "Nos1", "Calb2",
        "Tac1", "Foxp2", "C1ql3",
    ])

    marker_violins(
        neuronal_final, premammillary_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig4_premammillary-marker-violins", figdir="figure_04"
    )


def figure_S4():
    figdir = "figure_SX"
    gene_list = (
        "Pdgfra, Cspg4, Fyn, Tcf7l2, Ctps, Mal, Tspan2, Opalin, Apod, "
        "Klk6, Agt, Aqp4, Gfap, Slc7a10, Htra1, Itih3, C4b, Slc38a1, Rax, S100a6"
    ).split(", ")
    for gene in gene_list:
        fig, ax = plt.subplots(figsize=(2, 2))
        sc.pl.umap(
            nonneuronal_final, color=gene, 
            size=10, cmap=red_colormap, title="",
            ax=ax, show=False, use_raw=False
        )
        fig.delaxes(fig.get_children()[-1])
        ax.text(0.025, 0.975, gene, ha="left", va="top", size="small", style="italic", transform=ax.transAxes)
        ax.set_xlabel("")
        ax.set_ylabel("")
        fix_aspect_scatter_with_legend(fig)
        
        save_figure(fig, figdir, f"glial-lineage_{gene}")
        del fig
