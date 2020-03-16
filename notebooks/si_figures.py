#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from plotting.settings import *
from plotting.plot_funcs import *
from plotting.palettes import *
from utils.load_data import load_adata, save_figure, load_markers


def figure_S1a():
    full_no_dub = load_adata("global_with_dub")
    full_comb = full_no_dub.copy()
    sc.pp.normalize_per_cell(full_comb)
    sc.pp.log1p(full_comb)
    sc.pp.highly_variable_genes(full_comb, n_top_genes=2500, flavor="cell_ranger")
    sc.pp.pca(full_comb, n_comps=25, svd_solver="arpack", use_highly_variable=True)
    sc.pp.neighbors(full_comb, n_neighbors=15, metric="correlation")
    sc.tl.umap(full_comb, min_dist=0.5)

    full_comb.obs["chemistry"] = full_comb.obs.sampleid.map(
        {"AJ18003": "v2", "AJ18004": "v2", "AJ19001": "v3", "AJ19002": "v3",}
    )

    full_comb.obsm["X_umap"] = full_comb.obsm["X_umap"].dot(-np.eye(2))

    fig, axarr = plt.subplots(2, 3, figsize=(6, 4))
    params = dict(show=False, size=2)
    sc.pl.umap(
        full_comb,
        color="chemistry",
        ax=axarr[0, 0],
        palette=sns.xkcd_palette(["orange", "grass green"]),
        legend_loc=None,
        title="10X Chemistry",
        **params,
    )
    sc.pl.umap(
        full_no_dub,
        color="chemistry",
        ax=axarr[1, 0],
        palette=sns.xkcd_palette(["orange", "grass green"]),
        title="",
        **params,
    )
    sc.pl.umap(
        full_comb,
        color="sex",
        ax=axarr[0, 1],
        palette=sns.xkcd_palette(["electric blue", "bright red"]),
        legend_loc=None,
        title="Sex",
        **params,
    )
    sc.pl.umap(
        full_no_dub,
        color="sex",
        ax=axarr[1, 1],
        palette=sns.xkcd_palette(["electric blue", "bright red"]),
        title="",
        **params,
    )
    sc.pl.umap(
        full_comb,
        color="sample_name",
        ax=axarr[0, 2],  # palette=sns.xkcd_palette(["orange", "grass green"]),
        legend_loc=None,
        title="Sample",
        **params,
    )
    sc.pl.umap(
        full_no_dub,
        color="sample_name",
        ax=axarr[1, 2],  # palette=sns.xkcd_palette(["orange", "grass green"]),
        title="",
        **params,
    )

    legend_params = dict(
        bbox_to_anchor=(0.5, -0.0),
        loc="upper center",
        ncol=2,
        frameon=False,
        columnspacing=0.9,
        fontsize=7,
        markerscale=0.5,
    )
    for k in range(3):
        leg = axarr[1, k].legend(**legend_params)

    for ax in axarr.flat:
        ax.set_xlabel("")
        ax.set_ylabel("")
    axarr[0, 0].set_ylabel("No batch correction", size=8)
    axarr[1, 0].set_ylabel("With batch correction", size=8)
    axarr[0, 0].set_title("10X Chemistry", size=9)
    axarr[0, 1].set_title("Sex", size=9)
    axarr[0, 2].set_title("Sample", size=9)
    fix_aspect_scatter_with_legend(fig)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S01", "figS1a_batch-correction-comparison")


def figure_S1b():
    full_no_dub = load_adata("global_no_dub")
    neuronal_markers = ["Snap25", "Tubb3", "Elavl2", "Syp"]
    full_no_dub.obs["neuronal_average_exp"] = full_no_dub[:, neuronal_markers].X.mean(
        axis=1
    )

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        full_no_dub,
        color="neuronal_average_exp",
        palette=neuronal_nonneuronal_palette,
        title="Neuronal marker mean expression",
        show=False,
        color_map=red_colormap,
        ax=ax,
    )
    sns.despine(fig)
    fig.get_children()[-1].set_ylabel("Log-normalized expression")
    fix_aspect_scatter_with_cbar(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S01", "figS1b_neuronal-expression")


def figure_S1c():
    full_no_dub = load_adata("global_no_dub")

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        full_no_dub,
        color="classification_cluster",
        palette=neuronal_nonneuronal_palette,
        title="Neuronal Classification",
        show=False,
        ax=ax,
    )
    ax.legend(bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, loc="lower center")
    sns.despine(fig)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S01", "figS1c_neuronal-classification")


def figure_S2a():
    neuronal_final = load_adata("neuronal")
    neuronal_cluster_props = get_cluster_proportions(
        neuronal_final, "cluster_revised", sample_key="sample_name"
    )
    fig1 = plot_celltype_proportions(neuronal_cluster_props, neuronal_palette)
    save_figure(fig1, "figure_S02", "figS2a_neuronal_cluster_props")


def figure_S2b():
    neuronal_final = load_adata("neuronal")
    neuronal_sample_props = get_sample_proportions(
        neuronal_final, "cluster_revised", "sample_name"
    )
    fig2 = plot_celltype_proportions(
        neuronal_sample_props, sns.mpl_palette("tab20", 4)[::-1]
    )
    save_figure(fig2, "figure_S02", "figS2b_neuronal_sample_props")


def figure_S2c():
    nonneuronal_final = load_adata("nonneuronal")
    nonneuronal_cluster_props = get_cluster_proportions(
        nonneuronal_final, "cluster_final", sample_key="sample_name"
    )
    fig1 = plot_celltype_proportions(nonneuronal_cluster_props, nonneuronal_palette)
    save_figure(fig1, "figure_S02", "figS2c_nonneuronal_cluster_props")


def figure_S2d():
    nonneuronal_final = load_adata("nonneuronal")
    nonneuronal_sample_props = get_sample_proportions(
        nonneuronal_final, "cluster_final", "sample_name"
    )
    fig2 = plot_celltype_proportions(
        nonneuronal_sample_props, sns.mpl_palette("tab20", 4)[::-1]
    )
    save_figure(fig2, "figure_S02", "figS2d_nonneuronal_sample_props")


def figure_S3a():
    nonneuronal_final = load_adata("nonneuronal")
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        nonneuronal_final,
        color="cluster_final",
        show=False,
        ax=ax,
        title="Nonneuronal clusters",
        palette=nonneuronal_palette,
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S03", "fig3a_nonneuronal-clusters")


def figure_S3b():
    nonneuronal_final = load_adata("nonneuronal")

    nonneuronal_violin_genes = [
        "Ascl1",
        "Pdgfra",
        "Fyn",
        "Mag",
        "Mobp",
        "Mal",
        "Tspan2",
        "Ptgds",
        "Hcn2",
        "Fth1",
        "Agt",
        "Aqp4",
        "Gfap",
        "Rax",
        "Ccdc153",
        "Mrc1",
        "Tmem119",
        "Cx3cr1",
        "Rgs5",
        "Acta2",
        "Fxyd5",
        "Slc47a1",
        "Col1a1",
        "Dcn",
        "Igfbp2",
        "Pecam1",
        "Slc38a5",
    ]

    marker_violins(
        nonneuronal_final,
        nonneuronal_violin_genes,
        "cluster_final",
        nonneuronal_palette,
        filename="figS3b_nonneuronal-marker-violins",
        figdir="figure_S03",
    )


def figure_S3c():
    nonneuronal_final = load_adata("nonneuronal")

    nonneuronal_markers = load_markers("nonneuronal")

    (
        nonneuronal_heatmap_data,
        nonneuronal_col_colors,
        nonneuronal_row_colors,
    ) = get_heatmap_data(
        nonneuronal_final, nonneuronal_markers, "cluster_final", nonneuronal_palette
    )

    cg = sns.clustermap(
        nonneuronal_heatmap_data,
        z_score=0,
        vmin=-3,
        vmax=3,
        cmap=heatmap_cmap,
        xticklabels=False,
        yticklabels=False,
        row_cluster=False,
        col_cluster=False,
        col_colors=nonneuronal_col_colors,
        row_colors=nonneuronal_row_colors,
        robust=True,
        figsize=(4, 4),
        cbar_kws=dict(use_gridspec=True),
    )
    cg.ax_row_colors.set_ylabel("Genes", size="small")
    cg.ax_col_colors.set_title("Cells", size="small", zorder=100)
    cg.ax_col_colors.xaxis.tick_top()

    fix_heatmap_colorbar(cg, nonneuronal_final, nonneuronal_heatmap_data, heatmap_cmap)

    # cg.fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_heatmap_annotations(cg)
    save_figure(cg.fig, "figure_S03", "figS3c_nonneuronal-heatmap", dpi=600, ext="png")


def figure_S3d():
    nonneuronal_final = load_adata("nonneuronal")
    genes_umi_violins(
        nonneuronal_final,
        "cluster_final",
        nonneuronal_palette,
        figdir="figure_S03",
        filename="figS3d_nonneuronal-umi-genes-violins",
    )


def figure_S4a():
    nonneuronal_final = load_adata("nonneuronal")
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        nonneuronal_final[nonneuronal_final.obs.cluster_final.isin(list("123456"))],
        color="cluster_final",
        show=False,
        ax=ax,
        title="",
        palette=nonneuronal_palette,
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S04", "figS4a_oligo-lineage")


def figure_S4b():
    nonneuronal_final = load_adata("nonneuronal")
    gene_list = (
        ("Cspg4", "Mal"),
        ("Pdgfra", "Tspan2"),
        ("Fyn", "Opalin"),
        ("Tcf7l2", "Apod"),
        ("Ctps", "Klk6"),
    )
    for genes in gene_list:
        gene_label = "-".join(genes)
        scatter_1by2(
            nonneuronal_final,
            genes,
            selection=nonneuronal_final.obs.cluster_final.isin(list("123456")),
            label_ha="right",
            filename=f"figS4b_oligo-lineage_{gene_label}",
            figdir="figure_S04",
            add_marker_edges=False,
        )


def figure_S4c():
    nonneuronal_final = load_adata("nonneuronal")
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        nonneuronal_final[
            nonneuronal_final.obs.cluster_final.isin(["7", "8", "9", "10", "11"])
        ],
        color="cluster_final",
        show=False,
        ax=ax,
        title="",
        palette=nonneuronal_palette,
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S04", "figS4c_astro-lineage")


def figure_S4d():
    nonneuronal_final = load_adata("nonneuronal")
    gene_list = (
        ("Agt", "C4b"),
        ("Aqp4", "Gfap"),
        ("Slc7a10", "Slc38a1"),
        ("Htra1", "Rax"),
        ("Itih3", "S100a6"),
    )
    for genes in gene_list:
        gene_label = "-".join(genes)
        scatter_1by2(
            nonneuronal_final,
            genes,
            selection=nonneuronal_final.obs.cluster_final.isin(
                ["7", "8", "9", "10", "11"]
            ),
            label_ha="right",
            filename=f"figS4d_astro-lineage_{gene_label}",
            figdir="figure_S04",
            add_marker_edges=False,
        )


def figure_S5a():
    neuronal_final = load_adata("neuronal")
    genes = ["Tac2", "Pdyn", "Esr1", "Prlr"]
    scatter_2by2(
        neuronal_final,
        genes,
        filename="figS5a_arcuate-umap-markers",
        figdir="figure_S05",
        add_marker_edges=False,
    )


def figure_S5b():
    neuronal_final = load_adata("neuronal")
    arcuate_violin_genes = pd.Index(
        [
            "Slc32a1",
            "Slc17a6",
            "Prlr",
            "Tac2",
            "Tmem35a",
            "Esr1",
            "Nhlh2",
            "Pdyn",
            "Mrap2",
            "Inhbb",
            "Nr5a2",
            "Rxfp1",
            "Ucp2",
            "Lxn",
            "Trp53i11",
            "Fndc9",
            "Col2a1",
            "Kiss1",
            "Kcnk2",
        ]
    )

    marker_violins(
        neuronal_final,
        arcuate_violin_genes,
        "cluster_revised",
        neuronal_palette,
        filename="figS5b_arcuate-marker-violins",
        figdir="figure_S05",
    )


def figure_S6a():
    supramammillary = load_adata("supramammillary")
    genes = ["Gad1", "Gad2", "Slc32a1", "Slc17a6"]
    scatter_2by2(
        supramammillary,
        genes,
        filename="figS6a_cluster8_key-genes-1",
        figdir="figure_S06",
        add_marker_edges=False,
        label_va="top",
        label_ha="left",
        markersize=16,
    )


def figure_S6b():
    supramammillary = load_adata("supramammillary")

    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    sc.pl.umap(
        supramammillary, color="cluster_revised", show=False, ax=ax, title="", size=16
    )
    leg = ax.legend(
        bbox_to_anchor=(1.0, 0.5),
        loc="center left",
        ncol=1,
        frameon=False,
        fontsize="xx-small",
    )
    for handle in leg.legendHandles:
        handle._sizes = [16]
    fix_aspect_scatter_with_legend(fig)
    # add_scatter_borders(ax)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_S06", "fig6b_cluster8_subclusters")


def figure_S6c():
    supramammillary = load_adata("supramammillary")
    genes = ["Sema3c", "Inhba", "Rxfp1", "Cpne7"]
    scatter_2by2(
        supramammillary,
        genes,
        filename="figS6c_cluster8_key-genes-2",
        figdir="figure_S06",
        add_marker_edges=False,
        label_va="top",
        label_ha="left",
        markersize=16,
    )


def figure_S7a():
    neuronal_final = load_adata("neuronal")
    genes = ["Cck", "Synpr", "Tac1", "Nos1"]
    scatter_2by2(
        neuronal_final,
        genes,
        filename="figS7a_premammillary-umap-markers",
        figdir="figure_S07",
        add_marker_edges=False,
    )


def figure_S7b():
    neuronal_final = load_adata("neuronal")
    premammillary_violin_genes = pd.Index(
        [
            "Slc32a1",
            "Slc17a6",
            "Foxb1",
            "Cck",
            "Dlk1",
            "Synpr",
            "Spock1",
            "Nxph1",
            "Negr1",
            "Ar",
            "Nos1",
            "Calb2",
            "Tac1",
            "Foxp2",
            "C1ql3",
        ]
    )

    marker_violins(
        neuronal_final,
        premammillary_violin_genes,
        "cluster_revised",
        neuronal_palette,
        filename="figS7b_premammillary-marker-violins",
        figdir="figure_S07",
    )
