#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from io import StringIO
from pathlib import Path
from operator import add, sub
from string import ascii_uppercase

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo

from matplotlib.patches import Rectangle
from scipy.ndimage.filters import convolve
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from scanpy.plotting._anndata import _prepare_dataframe

from utils.load_data import save_figure
from plotting.palettes import *
from plotting.settings import *


def fix_aspect_scatter_with_legend(fig):
    changed_ylims = None
    for ax in fig.get_axes():
        ax.set_aspect("equal", "box", share=True)
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        xmid, ymid = add(*xlims) / 2, add(*ylims) / 2
        xlen, ylen = sub(*xlims[::-1]), sub(*ylims[::-1])
        if xlen < ylen:
            diff = (ylen - xlen) / 2
            ax.set_xlim(xlims[0] - diff, xlims[1] + diff)
        else:
            diff = (xlen - ylen) / 2
            changed_ylims = (ylims[0] - diff, ylims[1] + diff)
            ax.set_ylim(*changed_ylims)
    fig.subplots_adjust(left=0, right=1)


def fix_aspect_scatter_with_cbar(fig):
    changed_ylims = None
    for ax in fig.get_axes():
        if not ax.is_last_col():
            ax.set_aspect("equal", "box", share=True)
            xlims = ax.get_xlim()
            ylims = ax.get_ylim()
            xmid, ymid = add(*xlims) / 2, add(*ylims) / 2
            xlen, ylen = sub(*xlims[::-1]), sub(*ylims[::-1])
            if xlen < ylen:
                diff = (ylen - xlen) / 2
                ax.set_xlim(xlims[0] - diff, xlims[1] + diff)
            else:
                diff = (xlen - ylen) / 2
                changed_ylims = (ylims[0] - diff, ylims[1] + diff)
                ax.set_ylim(*changed_ylims)
        elif changed_ylims:
            ax.set_ylim(*changed_ylims)
    fig.subplots_adjust(left=0)  # , right=1)


def add_scatter_borders(ax):
    for child in ax.get_children()[::-1]:
        if not isinstance(child, matplotlib.collections.PathCollection):
            continue
        child.set_edgecolor(["0.7"])
        child.set_linewidth([0.1])
    legend = ax.get_legend()
    if legend:
        for legend_handle in legend.legendHandles:
            legend_handle.set_edgecolors("0.6")
            legend_handle.set_linewidths([0.2])


def scatter_1by2(
    adata,
    genes,
    filename=None,
    figdir=None,
    selection=None,
    add_marker_edges=True,
    label_va="top",
    label_ha="left",
    markersize=3,
):
    gs = plt.GridSpec(1, 2, wspace=0, hspace=0)
    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    axs = (ax1, ax2)

    plot_data = adata
    plot_params = dict(
        use_raw=False, show=False, title="", size=markersize, color_map=red_colormap
    )
    if selection is not None:
        plot_data = adata[selection]
        plot_params["size"] = 12

    x = dict(left=0.05, right=0.95)[label_ha]
    y = dict(bottom=0.05, top=0.95)[label_va]
    for ax, gene in zip(axs, genes):
        sc.pl.umap(plot_data, color=gene, ax=ax, **plot_params)
        ax.text(
            x, y, gene, ha=label_ha, va=label_va, size="small", transform=ax.transAxes
        )

    filtered_children = list(
        filter(lambda c: isinstance(c, matplotlib.axes.Axes), fig.get_children())
    )
    for k, child in enumerate(filtered_children[len(axs) :]):
        fig.delaxes(child)

    for ax in axs:
        ax.set_xlabel("")
        ax.set_ylabel("")
        fix_aspect_scatter_with_legend(fig)
        if add_marker_edges:
            add_scatter_borders(ax)
    ax1.set_ylabel("UMAP2")
    ax1.set_xlabel("UMAP1")
    ax2.set_xlabel("UMAP1")

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.0)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    if (filename is not None) and (figdir is not None):
        save_figure(fig, figdir, filename)


def scatter_2by2(
    adata,
    genes,
    filename=None,
    figdir=None,
    selection=None,
    add_marker_edges=True,
    label_va="top",
    label_ha="left",
    markersize=3,
):
    gs = plt.GridSpec(2, 2, wspace=0, hspace=0)
    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axs = (ax1, ax2, ax3, ax4)

    plot_data = adata
    plot_params = dict(
        use_raw=False, show=False, title="", size=markersize, color_map=red_colormap
    )
    if selection is not None:
        plot_data = adata[selection]
        plot_params["size"] = 12

    x = dict(left=0.05, right=0.95)[label_ha]
    y = dict(bottom=0.05, top=0.95)[label_va]
    for ax, gene in zip(axs, genes):
        sc.pl.umap(plot_data, color=gene, ax=ax, **plot_params)
        ax.text(
            x, y, gene, ha=label_ha, va=label_va, size="small", transform=ax.transAxes
        )

    filtered_children = list(
        filter(lambda c: isinstance(c, matplotlib.axes.Axes), fig.get_children())
    )
    for k, child in enumerate(filtered_children[len(axs) :]):
        fig.delaxes(child)

    for ax in axs:
        ax.set_xlabel("")
        ax.set_ylabel("")
        fix_aspect_scatter_with_legend(fig)
        if add_marker_edges:
            add_scatter_borders(ax)
    ax1.set_ylabel("UMAP2")
    ax3.set_ylabel("UMAP2")
    ax3.set_xlabel("UMAP1")
    ax4.set_xlabel("UMAP1")

    if len(genes) < 4:
        for k in range(4 - len(genes)):
            axs[::-1][k].set_axis_off()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.0)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    if (filename is not None) and (figdir is not None):
        save_figure(fig, figdir, filename)


def scatter_2by3(
    adata,
    genes,
    filename=None,
    figdir=None,
    selection=None,
    add_marker_edges=True,
    label_va="top",
    label_ha="left",
    markersize=3,
):
    gs = plt.GridSpec(2, 3, wspace=0, hspace=0)
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    axs = (ax1, ax2, ax3, ax4, ax5, ax6)

    plot_data = adata
    plot_params = dict(
        use_raw=False, show=False, title="", size=markersize, color_map=red_colormap
    )
    if selection is not None:
        plot_data = adata[selection]
        plot_params["size"] = 12

    x = dict(left=0.05, right=0.95)[label_ha]
    y = dict(bottom=0.05, top=0.95)[label_va]
    for ax, gene in zip(axs, genes):
        sc.pl.umap(plot_data, color=gene, ax=ax, **plot_params)
        ax.text(
            x, y, gene, ha=label_ha, va=label_va, size="small", transform=ax.transAxes
        )

    filtered_children = list(
        filter(lambda c: isinstance(c, matplotlib.axes.Axes), fig.get_children())
    )
    for k, child in enumerate(filtered_children[len(axs) :]):
        fig.delaxes(child)

    for ax in axs:
        ax.set_xlabel("")
        ax.set_ylabel("")
        fix_aspect_scatter_with_legend(fig)
        if add_marker_edges:
            add_scatter_borders(ax)
    ax1.set_ylabel("UMAP2")
    ax4.set_ylabel("UMAP2")
    ax4.set_xlabel("UMAP1")
    ax5.set_xlabel("UMAP1")
    ax6.set_xlabel("UMAP1")

    if len(genes) < 6:
        for k in range(6 - len(genes)):
            axs[::-1][k].set_axis_off()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.0)
    fig.tight_layout(rect=(0.025, 0.025, 0.975, 0.975))
    if (filename is not None) and (figdir is not None):
        save_figure(fig, figdir, filename)


def marker_violins(
    adata,
    genes,
    cluster_key,
    palette,
    cluster_name="Cluster",
    filename=None,
    figdir=None,
):
    L = len(genes)
    N = len(adata.obs[cluster_key].unique())
    fig, axs = plt.subplots(L, 1, figsize=(0.5 + 0.25*N, 0.32 * L))
    violin_data = pd.DataFrame(adata[:, genes].X.toarray(), columns=genes)
    violin_data[cluster_name] = adata.obs[cluster_key].values

    for k, (ax, gene) in enumerate(zip(axs.flat, genes)):
        sns.violinplot(
            data=violin_data,
            x=cluster_name,
            y=gene,
            palette=palette,
            ax=ax,
            inner=None,
            linewidth=0.4,
            scale="width",
            cut=0,
        )
        ax.set_ylabel(
            gene, rotation=0, ha="right", va="top", size="small", style="italic"
        )
        ax.set_yticks([])

        ax.set_ylim(0, ax.get_ylim()[1])

        if k < (L - 1):
            ax.set_xticks([])
            ax.set_xlabel(None)
            ax.xaxis.label.set_visible(False)
        else:
            ax.xaxis.set_tick_params(labelsize="small")
        sns.despine(ax=ax, left=True)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0, wspace=0)
    if (filename is not None) and (figdir is not None):
        save_figure(fig, figdir, filename)


def genes_umi_violins(
    adata, cluster_key, palette, filename=None, figdir=None, ylims=[50000, 10000]
):
    obs_to_plot = ["total_counts", "n_genes_by_counts"]
    obs_names = ["UIMs", "Genes"]
    yticks = [np.linspace(0, ylims[0], 6), np.linspace(0, ylims[1], 6)]
    L = len(obs_names)
    fig, axs = plt.subplots(2, 1, figsize=(4, 2))
    violin_data = pd.DataFrame(adata.obs[obs_to_plot].values, columns=obs_names)
    violin_data["Cluster"] = adata.obs[cluster_key].values

    for k, (ax, obs, ylim, ytick) in enumerate(zip(axs.flat, obs_names, ylims, yticks)):
        ax.grid(axis="y", lw=0.25, color="0.6")

        sns.violinplot(
            data=violin_data,
            x="Cluster",
            y=obs,
            palette=palette,
            ax=ax,
            inner=None,
            linewidth=0.4,
            scale="width",
            cut=0,
        )
        ax.set_ylabel(obs, rotation=90, ha="center", va="center", size="medium")
        ax.set_yticks(ytick)
        ax.set_ylim(0, ylim)
        ax.yaxis.set_tick_params(labelsize="small")

        if k < (L - 1):
            ax.set_xticks([])
            ax.set_xlabel("")
        else:
            ax.xaxis.set_tick_params(labelsize="small")
            ax.set_xlabel("Cluster", size="medium")
        sns.despine(ax=ax, left=False)
        ax.set_axisbelow(True)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.25, wspace=0)
    if (filename is not None) and (figdir is not None):
        save_figure(fig, figdir, filename)


def get_heatmap_data(adata, markers, cluster_key, palette):
    grpd = markers[markers.AUROC > 0.8].groupby(cluster_key).head(25)
    genes = grpd.gene_name  # .unique()
    duplicated = genes.duplicated()
    clusters = grpd[cluster_key].astype(int).values - 1
    genes = genes[~duplicated]
    clusters = clusters[~duplicated]

    # heatmap_data = nonneuronal_final[nonneuronal_final.obs.sort_values("cluster_final").index, :]
    data = adata[:, genes].X.toarray()
    sort_order = adata.obs[cluster_key].reset_index(drop=True).sort_values().index
    data = data[sort_order, :].T
    col_colors = np.array(palette)[adata.obs[cluster_key].cat.codes[sort_order]]
    row_colors = np.array(palette)[clusters]

    weights = np.array(
        [
            [0, 0, 1, 0, 0],
            [0, 2, 4, 2, 0],
            [1, 4, 8, 4, 1],
            [0, 2, 4, 2, 0],
            [0, 0, 1, 0, 0],
        ],
        dtype=np.float,
    )
    weights = weights / np.sum(weights[:])

    smoothed = convolve(data, weights, mode="constant")
    return smoothed, col_colors, row_colors


def get_matrix_data(adata, markers, cluster_key, palette):
    grpd = markers.groupby(cluster_key).head(25)
    # grpd = grpd[~grpd.gene_name.duplicated()]

    n_clusters = len(adata.obs[cluster_key].cat.categories)
    data = np.zeros((n_clusters, n_clusters))
    tmp = adata[:, grpd.gene_name].X.toarray()
    for i, clusteri in enumerate(adata.obs[cluster_key].cat.categories):
        gene_ind = (grpd[cluster_key] == clusteri).values
        for j, clusterj in enumerate(adata.obs[cluster_key].cat.categories):
            cell_ind = adata.obs[cluster_key].isin([clusterj])
            data[i, j] = tmp[cell_ind, :][:, gene_ind].mean()

    col_colors = np.array(palette)
    row_colors = np.array(palette)

    return data, col_colors, row_colors


def fix_heatmap_colorbar(cg, adata, heatmap_data, cmap="vlag"):
    # remove shitty stuff
    cg.cax.remove()
    cg.ax_col_dendrogram.remove()
    cg.ax_row_dendrogram.remove()
    # cg.fig.tight_layout()

    # here's the hack
    big_bbox = cg.ax_heatmap.get_position()
    width = (big_bbox.x1 - big_bbox.x0) * 0.8
    x0, x1 = big_bbox.x1 - width, big_bbox.x1

    cg.fig.subplots_adjust(bottom=0.075)
    new_cbar_ax = cg.fig.add_axes([x0, 0.05, width, 0.01])
    _fig, _ax = plt.subplots(figsize=(0.1, 0.1))
    sns.heatmap(
        heatmap_data,
        vmin=-3,
        vmax=3,
        cmap=cmap,
        ax=_ax,
        cbar_ax=new_cbar_ax,
        cbar_kws={"orientation": "horizontal"},
    )
    _ax.remove()

    new_cbar_ax.tick_params(axis="x", labelsize=4, length=2, width=0.3)
    new_cbar_ax.set_xlabel("Z-score transformed log(UMI + 1)", size=4)


def fix_heatmap_annotations(cg):
    pos = cg.ax_row_colors.get_position()
    pos.x0 = (pos.x0 + pos.x1) / 2
    cg.ax_row_colors.set_position(pos)

    pos = cg.ax_col_colors.get_position()
    pos.y1 = (pos.y0 + pos.y1) / 2
    cg.ax_col_colors.set_position(pos)


def get_cluster_proportions(
    adata, cluster_key="cluster_final", sample_key="replicate", drop_values=None
):

    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)

    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props


def get_sample_proportions(
    adata, cluster_key="cluster_final", sample_key="replicate", drop_values=None
):

    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([sample_key, cluster_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=cluster_key, index=sample_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)

    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props


def plot_sample_proportions(cluster_props, cluster_palette=None, xlabel_rotation=0):
    fig, ax = plt.subplots(dpi=300, figsize=(6, 4))
    fig.patch.set_facecolor("white")

    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, n_colors=len(cluster_palette), as_cmap=True
        )

    cluster_props.plot(kind="bar", stacked=True, ax=ax, legend=None, colormap=cmap)

    ax.legend(
        bbox_to_anchor=(1.01, 1), frameon=False, title="Sample", ncol=1, columnspacing=1
    )
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Proportion")
    fig.tight_layout()

    return fig


def plot_celltype_proportions(cluster_props, cluster_palette=None, xlabel_rotation=0):
    fig, ax = plt.subplots(figsize=(5, 3), dpi=300)
    fig.patch.set_facecolor("white")

    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, n_colors=len(cluster_palette), as_cmap=True
        )

    cluster_props.iloc[::-1,].plot(
        kind="barh", width=0.97, stacked=True, ax=ax, legend=None, colormap=cmap
    )

    sns.despine(fig, ax, bottom=True, left=True)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])

    ax.legend(
        bbox_to_anchor=(0.5, 0.05),
        loc="upper center",
        handletextpad=0.2,
        columnspacing=1,
        frameon=False,
        title="",
        ncol=10,
    )
    fig.tight_layout()
    return fig
