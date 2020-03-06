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


def figure_1c():
    fig, axarr = plt.subplots(1, 2, gridspec_kw=dict(hspace=0, wspace=0), sharex=True, sharey=True, figsize=(7, 4))
        
    params = dict(size=3, show=False, alpha=1., legend_fontsize=8)
    legend_params = dict(
        loc="upper center", bbox_to_anchor=(0.5, 0.), ncol=2, frameon=False,
        fontsize="small"
    )

    # by sex
    colors = sns.xkcd_palette(["electric blue", "bright red"])
    sc.pl.umap(full_no_dub, color="sex", title="", ax=axarr[0], palette=colors, **params)
    axarr[0].legend(**legend_params)

    # by 10x chemistry
    colors = sns.xkcd_palette(["orange", "grass green"])
    sc.pl.umap(full_no_dub, color="chemistry", title="", ax=axarr[1], palette=colors, **params)
    axarr[1].legend(**legend_params)

    for ax in axarr.flat:
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.set_aspect(1)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_01", "fig1_umap-sex-chemistry")


def figure_1d():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    axs = (ax1, ax2)

    params = dict(norm_hist=False, kde=False, hist_kws=dict(color=sns.xkcd_rgb["cerulean"], alpha=0.9))
    sns.distplot(full_no_dub.obs.total_counts, ax=ax1, bins=np.linspace(2000, 35000, 50), **params)
    sns.distplot(full_no_dub.obs.n_genes_by_counts, ax=ax2, **params)

    median1 = full_no_dub.obs.total_counts.median().astype(int)
    median2 = full_no_dub.obs.n_genes_by_counts.median().astype(int)
    ymax = 0.85
    ax1.axvline(median1, ymax=ymax, lw=1, ls="--", color="k", )
    ax2.axvline(median2, ymax=ymax, lw=1, ls="--", color="k")

    ax1.text(median1, ymax*1.02*ax1.get_ylim()[1], str(median1), transform=ax1.transData, ha="center", size="small")
    ax2.text(median2, ymax*1.02*ax2.get_ylim()[1], str(median2), transform=ax2.transData, ha="center", size="small")

    ax1.set_xlim(0, 35000)
    ax2.set_xlim(0, 10000)

    ax1.set_xlabel("No. of transcripts")
    ax2.set_xlabel("No. of genes")
    ax1.set_ylabel("No. of cells")
    ax2.set_ylabel("No. of cells")

    for ax in axs:
        sns.despine(fig, ax)

    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_01", "fig1_qc-histograms")


def figure_1f():
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        full_no_dub, 
        color="cluster_revised",  
        palette=main_palette,
        title="",
        show=False,
        ax=ax
    )
    legend = ax.legend(bbox_to_anchor=(1, 1), ncol=1, frameon=False, loc="upper left", fontsize=7.6)
    for leg in legend.legendHandles:
        leg._sizes = [20]
    sns.despine(fig)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_01", "fig1_all-clusters")


def figure_1e():
    full_markers = sc.tl.find_marker_genes(full_no_dub, "cluster_revised", log_fold_change=0.75)
    full_heatmap_data, full_col_colors, full_row_colors = get_heatmap_data(
        full_no_dub, full_markers, "cluster_revised", main_palette
    )
    cg = sns.clustermap(
        full_heatmap_data, z_score=0, vmin=-3, vmax=3, cmap=heatmap_cmap,
        xticklabels=False, yticklabels=False,
        row_cluster=False, col_cluster=False,
        col_colors=full_col_colors, 
        row_colors=full_row_colors,
        robust=True, figsize=(4, 4), cbar_kws=dict(use_gridspec=True)
    )
    cg.ax_row_colors.set_ylabel("Genes", size="small")
    cg.ax_col_colors.set_title("Cells", size="small", zorder=100)

    fix_heatmap_colorbar(cg, full_no_dub, full_heatmap_data, cmap=heatmap_cmap)

    cg.ax_col_colors.xaxis.tick_top()

    cg.fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_heatmap_annotations(cg)
    save_figure(cg.fig, "figure_01", "fig1_main-heatmap", dpi=600, ext="png")


def figure_2a():
    from sklearn.mixture import GaussianMixture
    mixture = GaussianMixture(n_components=3, random_state=3)
    mixture.fit(neuronal_final[:, ["Slc17a6", "Slc32a1", "Hdc"]].X.toarray())
    gaba_vs_glut_vs_hdc = mixture.predict(neuronal_final[:, ["Slc17a6", "Slc32a1", "Hdc"]].X.toarray())

    neuronal_final.obs["gaba_vs_glut_vs_hdc"] = gaba_vs_glut_vs_hdc
    neuronal_final.obs["gaba_vs_glut_vs_hdc_cluster"] = neuronal_final.obs.groupby("cluster_revised")["gaba_vs_glut_vs_hdc"].transform(np.mean).round().astype(int).map({
        0: "HA", 2: "Glut", 1: "GABA"
    })
    gaba_vs_glut_palette = sns.xkcd_palette(["light red", "medium blue", "electric purple", "light grey"])
    gaba_vs_glut_vs_hdc_palette = sns.xkcd_palette(["yellow orange", "medium blue", "tree green"])

    gs = plt.GridSpec(1, 5)
    fig = plt.figure(figsize=(18, 4))
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])
    ax4 = fig.add_subplot(gs[0,3])
    ax6 = fig.add_subplot(gs[0,4])
    axs = (ax1, ax2, ax3, ax4, ax6)

    genes = ["Slc17a6", "Slc32a1", "Gad1", "Hdc"]
    for ax, gene in zip(axs, genes):
        sc.pl.umap(neuronal_final, color=gene, use_raw=False, show=False, ax=ax, title="", color_map=red_colormap,
                   size=10)
        ax.text(0.05, 0.95, gene, ha="left", va="top", size="large", transform=ax.transAxes)

    sc.pl.umap(
        neuronal_final, color="gaba_vs_glut_vs_hdc_cluster", 
        ax=ax6, show=False, title="", legend_loc="right margin",
        palette=gaba_vs_glut_vs_hdc_palette,
        size=10
    )
    ax6.legend(frameon=False, loc="upper right", labelspacing=0, borderaxespad=0.1, markerscale=0.7)
        
    filtered_children = list(filter(lambda c: isinstance(c, matplotlib.axes.Axes), fig.get_children()))
    for k, child in enumerate(filtered_children[len(axs):]):
        fig.delaxes(child)
        
    for ax in axs:
        ax.set_xlabel("")
        ax.set_ylabel("")
        fix_aspect_scatter_with_legend(fig)

    ax1.set_xlabel("UMAP1")
    ax2.set_xlabel("UMAP1")
    ax3.set_xlabel("UMAP1")
    ax4.set_xlabel("UMAP1")
    ax6.set_xlabel("UMAP1")
    ax6.set_ylabel("UMAP2")
    ax6.yaxis.set_label_position("right")

    fig.subplots_adjust(wspace=0.0)
    save_figure(fig, "figure_03", "fig3_gaba-vs-glut")


def figure_2b():
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    sc.pl.umap(
        neuronal_final, color="cluster_revised", show=False, ax=ax, size=40, 
        title="Neuronal clusters", palette=neuronal_palette
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_03", "fig3_neuronal-clusters")


def figure_2c():
    neuronal_marker_gene_list = sc.tl.find_marker_genes(neuronal_final, cluster_key="cluster_revised", log_fold_change=0.7)

    neuronal_heatmap_data, neuronal_col_colors, neuronal_row_colors = get_heatmap_data(
        neuronal_final, neuronal_marker_gene_list, "cluster_revised", neuronal_palette
    )

    cg = sns.clustermap(
        neuronal_heatmap_data, z_score=0, vmin=-3, vmax=3, cmap=heatmap_cmap,
        xticklabels=False, yticklabels=False,
        row_cluster=False, col_cluster=False,
        col_colors=neuronal_col_colors, 
        row_colors=neuronal_row_colors,
        robust=True, figsize=(4, 4), cbar_kws=dict(use_gridspec=True)
    )
    cg.ax_row_colors.set_ylabel("Genes", size="small")
    cg.ax_col_colors.set_title("Cells", size="small", zorder=100)
    cg.ax_col_colors.xaxis.tick_top()

    fix_heatmap_colorbar(cg, neuronal_final, neuronal_heatmap_data, heatmap_cmap)

    cg.fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_heatmap_annotations(cg)
    save_figure(cg.fig, "figure_03", "fig3_neuronal-heatmap", dpi=600, ext="png")

    fig, ax = plt.subplots(figsize=(1, 1.5))
    for k, (cluster, inds) in enumerate(neuronal_final.obs.groupby("cluster_revised")):
        x, y = neuronal_final[inds.index, :].obsm["X_umap"].T
        ax.scatter(x, y, c=neuronal_palette[k], marker="s", label=cluster)
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
    save_figure(fig, "neuronal-heatmap-legend")


def figure_2d():
    neuronal_violin_genes = pd.Index([
        "Slc17a6", "Slc32a1", "Gad1", "Gad2",
        "Nts", "Gpr83", "Onecut2", "Cox6a2", "Tac2", "Synpr", "Ar", "Nr4a2",
        "Calcr", "C1ql3", "Gabre", "Asb4", "Nrxn3", "Sox14", "Six3", "Rxfp1",
        "Hdc", "Npy2r", "Sst", "Otp", "Fgf1",
    ])

    marker_violins(
        neuronal_final, neuronal_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig3_neuronal-marker-violins_v2", figdir="figure_03"
    )


def figure_2e():
    genes_umi_violins(
        neuronal_final, "cluster_revised", neuronal_palette, 
        filename="fig3_neuronal-umi-genes-violins_v2", figdir="figure_03",
        ylims=[80000, 10000]
    )


def figure_3a():
    genes = ["Tac1", "Nos1", "Calb2", "Foxp2"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig4_premammillary-cluster7-umap-markers", figdir="figure_04", add_marker_edges=False
    )


def figure_3b():
    premammillary_cluster7_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6", "Calb2", "Foxp2", "Htr2c", "Nos1", "Glra3",
        "Tac1", "C1ql3", "Irs4",
    ])

    marker_violins(
        neuronal_final, premammillary_cluster7_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig4_premammillary-cluster7-marker-violins", figdir="figure_04",
    )


def figure_3e():
    cluster_7_markers = "Tac1, Slc17a6, Slc6a3, Ddc, Slc18a2, Th".split(", ")

    scatter_2by3(
        neuronal_final, cluster_7_markers, 
        filename="fig4_cluster7-umap-markers", figdir="figure_04", 
        add_marker_edges=False, selection=neuronal_final.obs.cluster_revised.isin(["7"]),
        label_va="top", label_ha="right"
    )


def figure_4a():
    genes = ["Hdc", "Slc18a2", "Wif1", "Maob"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig5_tuberomammillary-umap-markers", figdir="figure_05", add_marker_edges=False
    )


def figure_4c():
    tuberomammillary_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6", "Gad1", 
        "Slc18a2", "Hdc", "Wif1",
        "Maoa", "Maob", "Msrb2", "Itm2a",
        "Tspan12", "Bsx", "Prph", "Sncg",
    ])

    marker_violins(
        neuronal_final, tuberomammillary_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig5_tuberomammillary-marker-violins", figdir="figure_05"
    )


def figure_5a():
    genes = ["Tac2", "Tcf4", "Cplx1", "Pvalb"]
    scatter_2by2(
        neuronal_final, genes, 
        "fig6_lateral-mammillary-umap-markers", figdir="figure_06", add_marker_edges=False
    )

def figure_5b():
    lateral_mammillary_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6", "Tac2", 
        "Tcf4", "Pvalb",
        "Fgf1",
        "Cnr1", "Cplx1", "Cbln2",
        "Cabp7", "Parm1", 
        "Inf2", "Nefm", "Myo1a", "Syt2"
    ])

    marker_violins(
        neuronal_final, lateral_mammillary_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig6_lateral-mammillary-marker-violins", figdir="figure_06"
    )


def figure_6a():
    genes = ["Cartpt", "Foxb1", "Cck", "Adcyap1"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig7_mammillary-umap-global-markers", figdir="figure_07", add_marker_edges=False
    )


def figure_6b():
    mammillary_violin_genes = pd.Index([
        "Slc32a1", "Slc17a6",
        "Ctxn3", "Gpr83", "Rprm",  
        "Cck", "Adcyap1", "Fam19a1", "Slc24a2", "Foxb1",
        "Cpne9", "Cox6a2", "Cnih3", "Pvalb", 
        "Nts", "Onecut2", "Tac2",  "Cxcl14",
    ])

    marker_violins(
        neuronal_final, mammillary_violin_genes, "cluster_revised", neuronal_palette,
        filename="fig7_mammillary-marker-violins_v2", figdir="figure_07"
    )

def figure_7a():
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sc.pl.umap(
        neuronal_final[neuronal_final.obs.cluster_revised.isin(list("12345"))], 
        color="cluster_revised", show=False, ax=ax, title="Mammillary clusters",
        palette=neuronal_palette,
    )
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", ncol=1, frameon=False)
    fix_aspect_scatter_with_legend(fig)
    #add_scatter_borders(ax)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    save_figure(fig, "figure_07", "fig7_mammillary-umap-clusters")


def figure_7b():
    def get_heatmap_data_tmp(adata, markers, cluster_key, palette):
        grpd = markers.groupby(cluster_key).head(100)
        genes = grpd.gene_name#.unique()
        duplicated = genes.duplicated()
        clusters = grpd[cluster_key].astype(int).values - 1
        genes = genes[~duplicated]
        clusters = clusters[~duplicated]

        #heatmap_data = nonneuronal_final[nonneuronal_final.obs.sort_values("cluster_final").index, :]
        data = adata[:, genes].X.toarray()
        sort_order = adata.obs[cluster_key].reset_index(drop=True).sort_values().index
        data = data[sort_order, :].T
        col_colors = np.array(palette)[adata.obs[cluster_key].cat.codes[sort_order]]
        row_colors = np.array(palette)[clusters]

        weights = np.array([[0, 0, 1, 0, 0],
                            [0, 2, 4, 2, 0],
                            [1, 4, 8, 4, 1],
                            [0, 2, 4, 2, 0],
                            [0, 0, 1, 0, 0]],
                           dtype=np.float)
        weights = weights / np.sum(weights[:])

        smoothed = convolve(data, weights, mode='constant')
        return smoothed, col_colors, row_colors

    mammillary = neuronal_final[neuronal_final.obs.cluster_revised.isin(list("12345")), :].copy()
    mammillary_heatmap_genes = pd.Index([
        "Nts", "Alcam", "Col25a1", "Ctxn3", 
        "Spock3", "Serpini1", "Tshz2", "Trp53i11", "Grin3a", "Hpcal1", "Zbtb20",
        "Onecut2", "Rbms3", "Slc24a2", "Pvalb", 
        "Calb1", "Nos1", 
        "Tac2", "Cadm1", "Cxcl14",
    ])
    mammillary_marker_genes = pd.DataFrame({
        "gene_name": mammillary_heatmap_genes,
        "cluster_revised": ["1"]*4 + ["2"]*7 + ["3"]*4 + ["4"]*2 + ["5"] * 3,
        "garbage": np.arange(len(mammillary_heatmap_genes))
    })

    mammillary_heatmap_data, mammillary_col_colors, mammillary_row_colors = get_heatmap_data_tmp(
        mammillary, mammillary_marker_genes, "cluster_revised", neuronal_palette
    )

    cg = sns.clustermap(
        mammillary_heatmap_data, z_score=0, vmin=-3, vmax=3, cmap=heatmap_cmap,
        xticklabels=False, yticklabels=mammillary_heatmap_genes,
        row_cluster=False, col_cluster=False,
        col_colors=mammillary_col_colors, 
        row_colors=mammillary_row_colors,
        robust=True, figsize=(4, 4), cbar_kws=dict(use_gridspec=True)
    )
    cg.ax_row_colors.set_ylabel("Genes", size="small")
    cg.ax_col_colors.set_title("Cells", size="small", zorder=100)
    cg.ax_heatmap.tick_params(axis="y", which="major", labelsize="xx-small")
    cg.ax_col_colors.xaxis.tick_top()

    fix_heatmap_colorbar(cg, mammillary, mammillary_heatmap_data, heatmap_cmap)
    cg.fig.subplots_adjust(right=0.78)

    cg.fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    fix_heatmap_annotations(cg)
    save_figure(cg.fig, "last_minute", "fig7_mammillary-heatmap", dpi=600, ext="png")


def figure_7c():
    genes = ["Col25a1", "Onecut2", "Calb1", "Alcam"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig7_mammillary-umap-subpopulation-markers-add1", figdir="figure_07", add_marker_edges=False,
        selection=neuronal_final.obs.cluster_revised.isin(list("12345")), label_va="bottom"
    )
    genes = ["Gpr83", "Spock3", "Serpini1", "Cxcl14"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig7_mammillary-umap-subpopulation-markers-add2", figdir="figure_07", add_marker_edges=False,
        selection=neuronal_final.obs.cluster_revised.isin(list("12345")), label_va="bottom"
    )
    genes = ["Pvalb", "Slc24a2", "Nos1", "Calb1"]
    scatter_2by2(
        neuronal_final, genes, 
        filename="fig7_mammillary-umap-subpopulation-markers-slc24a2", figdir="last_minute", add_marker_edges=False,
        selection=neuronal_final.obs.cluster_revised.isin(list("12345")), label_va="bottom"
    )





def fix_heatmap_colorbar(cg, heatmap_data):
    # remove shitty stuff
    cg.cax.remove()
    cg.ax_col_dendrogram.remove()
    cg.ax_row_dendrogram.remove()
    cg.fig.tight_layout()
    
    # here's the hack
    big_bbox = cg.ax_heatmap.get_position()
    height = big_bbox.y1 - big_bbox.y0
    
    cg.fig.subplots_adjust(right=0.9)
    new_cbar_ax = cg.fig.add_axes([0.91, 0.125, 0.025, height])
    _fig, _ax = plt.subplots(figsize=(0.1, 0.1))
    sns.heatmap(heatmap_data, vmin=-3, vmax=3, cmap="vlag", ax=_ax, cbar_ax=new_cbar_ax)
    _ax.remove()
    
    new_cbar_ax.tick_params(axis="y", labelsize="xx-small", length=2, width=0.3)
    new_cbar_ax.set_ylabel("Z-score transformed log(UMI + 1)", size="xx-small")





# ## LAST MINUTE FIGURES
# 
# 1. For Fig. 4b, could I please get a version that also includes Hrh3 and Hcrtr2?
# 2. For Fig. 6b, could I get a version that also includes Cartpt?
# 3. For Fig. 7b, Could I get a version that includes Calb1 to add to cluster 4.
# 4. for Fig. 7c (UMAP for cluster 3), could I get a panel for Slc24a2? (The one I have, Onecut2 doesn’t have ISH data)
# 5. For new suppl Fig. S2, could you make versions that have Ariel font please?  Otherwise I’m changing all the text by hand. 
# 6. for suppl. Fig. S7a, could I have UMAPs (same proportion) for Nxph1 and Ar?

# In[59]:


# FIGURE 4B
tuberomammillary_violin_genes = pd.Index([
    "Slc32a1", "Slc17a6", "Gad1", 
    "Slc18a2", "Hdc", "Wif1",
    "Maoa", "Maob", "Msrb2", "Itm2a", "Hcrtr2", "Hrh3", 
    "Bsx", "Tspan12", "Prph", 
    "Sncg",
])

marker_violins(
    neuronal_final, tuberomammillary_violin_genes, "cluster_revised", neuronal_palette,
    filename="fig4_tuberomammillary-marker-violins", figdir="last_minute"
)

#FIGURE 6B
mammillary_violin_genes = pd.Index([
    "Slc32a1", "Slc17a6",
    "Ctxn3", "Gpr83", "Rprm",  
    "Cck", "Adcyap1", "Fam19a1", "Slc24a2", "Foxb1",
    "Cpne9", "Cox6a2", "Cnih3", "Pvalb", "Cartpt",
    "Nts", "Onecut2", "Tac2",  "Cxcl14",
])

marker_violins(
    neuronal_final, mammillary_violin_genes, "cluster_revised", neuronal_palette,
    filename="fig7_mammillary-marker-violins", figdir="last_minute"
)

genes = ["Pvalb", "Slc24a2", "Nos1", "Calb1"]
scatter_2by2(
    neuronal_final, genes, 
    filename="fig7_mammillary-umap-subpopulation-markers-slc24a2", figdir="last_minute", add_marker_edges=False,
    selection=neuronal_final.obs.cluster_revised.isin(list("12345")), label_va="bottom"
)

genes = ["Cck", "Synpr", "Nxph1", "Tac1", "Nos1", "Ar"]
scatter_2by3(
    neuronal_final, genes, 
    filename="fig4_premammillary-umap-markers", figdir="last_minute", add_marker_edges=False
)
