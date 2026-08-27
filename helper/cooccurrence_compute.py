"""Cell-type vs. condtions spatial co-occurrence/distance analysis and plotting helpers."""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq

def compute_celltype_cooccurrence(adata, cell_types, anchors, interval, sample_col="punch_id", condition_col="Condition", label_col="annotation", min_cells=1):
    """Per-sample squidpy co-occurrence ratio of each cell type relative to each anchor, across distance bins.

    Parameters
    ----------
    adata : anndata.AnnData
        Object containing the `sample_col` and `condition_col` columns and the categorical
        label column given by `label_col`.
    cell_types : list of str
        Target cell-type labels to evaluate (the "j" side of the ratio).
    anchors : list of str
        Anchor labels to condition on (the "i" side of the ratio), e.g. transf classes.
    interval : numpy.ndarray
        Distance bin edges (µm) passed to ``squidpy.gr.co_occurrence``.
    sample_col : str, default "punch_id"
        ``adata.obs`` column identifying the sample; one squidpy run per unique value.
    condition_col : str, default "Condition"
        ``adata.obs`` column holding the experimental condition (assumed constant within
        each sample).
    label_col : str, default "annotation"
        ``adata.obs`` column holding the cell-type/anchor labels.
    min_cells : int, default 1
        Minimum cells a (sample, label) combination must have to be included; combinations
        below this are skipped entirely, never zero-filled.

    Returns
    -------
    pandas.DataFrame
        One row per (sample, anchor, cell_type, distance) with columns ``sample``,
        ``condition``, ``anchor``, ``cell_type``, ``distance``, ``ratio``. Empty if no sample
        has both an anchor and a cell type reaching `min_cells`.
    """
    records = []
    for s in adata.obs[sample_col].unique():
        ad_s = adata[adata.obs[sample_col] == s]
        condition = ad_s.obs[condition_col].iloc[0]
        ann = ad_s.obs[label_col].astype(str)
        counts = ann.value_counts()
        present = [lbl for lbl in cell_types + list(anchors) if counts.get(lbl, 0) >= min_cells]
        if not any(a in present for a in anchors) or not any(ct in present for ct in cell_types) or len(present) < 2:
            continue
        sub = ad_s[ann.isin(present)].copy()
        sub.obs[label_col] = sub.obs[label_col].astype(str).astype("category")
        occ, dist = sq.gr.co_occurrence(
            sub, cluster_key=label_col, interval=interval, copy=True,
            n_jobs=16, show_progress_bar=False,
        )
        cats = list(sub.obs[label_col].cat.categories)
        x = dist[1:]
        for anchor in anchors:
            if anchor not in cats:
                continue
            i = cats.index(anchor)
            for ct in cell_types:
                if ct not in cats or ct == anchor:
                    continue
                j = cats.index(ct)
                for d, ratio in zip(x, occ[i, j, :]):
                    records.append({
                        "sample": s, "condition": condition,
                        "anchor": anchor, "cell_type": ct, "distance": d, "ratio": ratio,
                    })
    return pd.DataFrame(records)


def plot_condition_cooccurrence(condition, summary_df, anchor_order, cell_type_colors, pub_rc):
    """Plot mean +/- SEM co-occurrence ratio vs. distance for one condition, faceted by anchor.

    Draws a 1 x len(anchor_order) grid; each panel overlays every key of `cell_type_colors`
    as a solid line (mean) with a shaded SEM band, plus a dashed reference line at ratio = 1.

    Parameters
    ----------
    condition : str
        Condition to plot (must match a value in ``summary_df["condition"]``).
    summary_df : pandas.DataFrame
        Aggregated co-occurrence summary with columns ``cell_type``, ``anchor``, ``condition``,
        ``distance``, ``mean``, ``sem``, ``n``.
    anchor_order : list of str
        Anchor labels and left-to-right panel order.
    cell_type_colors : dict
        Mapping of cell-type name -> line/fill color; also sets the overlay order.
    pub_rc : dict
        Matplotlib rcParams applied via ``mpl.rc_context`` for publication styling.

    Returns
    -------
    None
        Displays the figure with ``plt.show()``.
    """
    fig, axes = plt.subplots(1, len(anchor_order), figsize=(2.0 * len(anchor_order), 2.3),
                              sharex=True, sharey=True, facecolor="white")
    with mpl.rc_context(pub_rc):
        for i, anchor in enumerate(anchor_order):
            ax = axes[i]
            sub_anchor = summary_df[(summary_df["anchor"] == anchor) & (summary_df["condition"] == condition)]
            for cell_type, color in cell_type_colors.items():
                sub = sub_anchor[sub_anchor["cell_type"] == cell_type].sort_values("distance")
                if sub.empty:
                    continue
                ax.plot(sub["distance"], sub["mean"], color=color, lw=1.4, label=cell_type)
                ax.fill_between(sub["distance"], sub["mean"] - sub["sem"], sub["mean"] + sub["sem"],
                                 color=color, alpha=0.25, linewidth=0)

            ax.axhline(1, ls="--", lw=0.8, color="#888888", alpha=0.7, zorder=0)
            ax.set_ylim(0, 3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_title(anchor, fontsize=7, fontweight="bold")
            ax.set_xlabel("Distance (\u00b5m)", fontsize=7)
            if i == 0:
                ax.set_ylabel("Co-occurrence ratio", fontsize=7)

        # Shared legend below the panels (proxy handles so every cell type shows even if absent from a panel)
        legend_handles = [mpl.lines.Line2D([0], [0], color=color, lw=1.4) for color in cell_type_colors.values()]
        fig.legend(legend_handles, list(cell_type_colors.keys()), fontsize=6, frameon=False,
                   loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=min(len(cell_type_colors), 4))

        fig.suptitle(condition, fontsize=9, fontweight="bold", y=1.03)
        fig.tight_layout(w_pad=1.5, rect=(0, 0.16, 1, 1))
        plt.show()


def plot_all_conditions_cooccurrence(summary_df, condition_order, anchor_order, cell_type_colors, pub_rc):
    """Plot mean +/- SEM co-occurrence ratio vs. distance for every condition in one figure.

    Draws a len(condition_order) x len(anchor_order) grid (rows = conditions, columns =
    anchors); each panel overlays every key of `cell_type_colors` as a solid line (mean) with
    a shaded SEM band, plus a dashed reference line at ratio = 1. Anchor names title the top
    row only; condition names label the left column only.

    Parameters
    ----------
    summary_df : pandas.DataFrame
        Aggregated co-occurrence summary with columns ``cell_type``, ``anchor``, ``condition``,
        ``distance``, ``mean``, ``sem``, ``n``.
    condition_order : list of str
        Conditions and top-to-bottom row order.
    anchor_order : list of str
        Anchor labels and left-to-right column order.
    cell_type_colors : dict
        Mapping of cell-type name -> line/fill color; also sets the overlay order.
    pub_rc : dict
        Matplotlib rcParams applied via ``mpl.rc_context`` for publication styling.

    Returns
    -------
    None
        Displays the figure with ``plt.show()``.
    """
    n_rows, n_cols = len(condition_order), len(anchor_order)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.0 * n_cols, 2.1 * n_rows),
                              sharex=True, sharey=True, facecolor="white", squeeze=False)
    with mpl.rc_context(pub_rc):
        for r, condition in enumerate(condition_order):
            for c, anchor in enumerate(anchor_order):
                ax = axes[r][c]
                sub_anchor = summary_df[(summary_df["anchor"] == anchor) & (summary_df["condition"] == condition)]
                for cell_type, color in cell_type_colors.items():
                    sub = sub_anchor[sub_anchor["cell_type"] == cell_type].sort_values("distance")
                    if sub.empty:
                        continue
                    ax.plot(sub["distance"], sub["mean"], color=color, lw=1.4, label=cell_type)
                    ax.fill_between(sub["distance"], sub["mean"] - sub["sem"], sub["mean"] + sub["sem"],
                                     color=color, alpha=0.25, linewidth=0)

                ax.axhline(1, ls="--", lw=0.8, color="#888888", alpha=0.7, zorder=0)
                ax.set_ylim(0, 3)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                if r == 0:
                    ax.set_title(anchor, fontsize=7, fontweight="bold")
                if r == n_rows - 1:
                    ax.set_xlabel("Distance (\u00b5m)", fontsize=7)
                if c == 0:
                    ax.set_ylabel(condition, fontsize=7.5, fontweight="bold")

        # Shared legend below all panels (proxy handles so every cell type shows even if absent from a panel)
        legend_handles = [mpl.lines.Line2D([0], [0], color=color, lw=1.4) for color in cell_type_colors.values()]
        fig.legend(legend_handles, list(cell_type_colors.keys()), fontsize=6, frameon=False,
                   loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=min(len(cell_type_colors), 4))

        fig.suptitle("Co-occurrence ratio", fontsize=9, fontweight="bold", y=1.01)
        fig.tight_layout(w_pad=1.5, h_pad=1.2, rect=(0, 0.05, 1, 1))
        plt.show()


def plot_violin_matrix(values_df, target_cell_types, transf_types, week_order, celltype_colors, transf_colors, week_colors, near_zone_max, pub_rc):
    """Plot a matrix of violin plots of near-anchor co-occurrence ratios.

    Rows = `target_cell_types`, columns = `transf_types` x `week_order`. Row and column group
    labels are drawn as solid colored header segments with the label text inside (not as tick
    text or a title). Each panel is a violin of the raw ``ratio`` values in `values_df` for that
    (cell type, transf class, week) combination; a panel with 0 points is drawn flat at 0 (never
    zero-filled into a violin), and a panel with exactly 1 point is drawn as a single dot.

    Parameters
    ----------
    values_df : pandas.DataFrame
        Raw (unaveraged) co-occurrence values with columns ``cell_type``, ``anchor``, ``week``,
        ``ratio`` (e.g. the output of `compute_near_transf_values`).
    target_cell_types : list of str
        Row order (cell types).
    transf_types : list of str
        Outer column grouping (transf classes).
    week_order : list of str
        Inner column grouping within each transf class.
    celltype_colors : dict
        Mapping of cell type -> row-header background color.
    transf_colors : dict
        Mapping of transf class -> column-header and violin-body color.
    week_colors : dict
        Mapping of week -> inner column-header background color.
    near_zone_max : float
        Upper distance (µm) used only in the figure suptitle.
    pub_rc : dict
        Matplotlib rcParams applied via ``mpl.rc_context`` for publication styling.

    Returns
    -------
    None
        Displays the figure with ``plt.show()``.
    """
    n_rows = len(target_cell_types)
    cols = [(transf, week) for transf in transf_types for week in week_order]
    n_cols = len(cols)
    n_weeks = len(week_order)

    with mpl.rc_context(pub_rc):
        fig = plt.figure(figsize=(max(4.2, 0.5 * n_cols + 2.0), 0.55 * n_rows + 1.6), facecolor="white")
        outer_gs = fig.add_gridspec(2, 1, height_ratios=[0.75, n_rows], hspace=0.04)
        header_gs = outer_gs[0].subgridspec(2, n_cols + 2, width_ratios=[1.6, 0.9] + [1] * n_cols, hspace=0.15, wspace=0.15)
        data_gs = outer_gs[1].subgridspec(n_rows, n_cols + 2, width_ratios=[1.6, 0.9] + [1] * n_cols, hspace=0.6, wspace=0.15)

        for li, transf in enumerate(transf_types):
            ax_h = fig.add_subplot(header_gs[0, 2 + li * n_weeks: 2 + (li + 1) * n_weeks])
            ax_h.set_facecolor(transf_colors[transf])
            ax_h.set_xticks([]); ax_h.set_yticks([])
            for spine in ax_h.spines.values():
                spine.set_visible(False)
            ax_h.text(0.5, 0.5, transf, ha="center", va="center", fontsize=6.5, fontweight="bold",
                      color=text_color_for_bg(transf_colors[transf]), transform=ax_h.transAxes)

        for ci, (transf, week) in enumerate(cols):
            ax_h = fig.add_subplot(header_gs[1, 2 + ci])
            ax_h.set_facecolor(week_colors[week])
            ax_h.set_xticks([]); ax_h.set_yticks([])
            for spine in ax_h.spines.values():
                spine.set_visible(False)
            ax_h.text(0.5, 0.5, week, ha="center", va="center", fontsize=5.5, fontweight="bold",
                      color=text_color_for_bg(week_colors[week]), transform=ax_h.transAxes)

        for ri, cell_type in enumerate(target_cell_types):
            ax_h = fig.add_subplot(data_gs[ri, 0])
            ax_h.set_facecolor(celltype_colors[cell_type])
            ax_h.set_xticks([]); ax_h.set_yticks([])
            for spine in ax_h.spines.values():
                spine.set_visible(False)
            ax_h.text(0.5, 0.5, cell_type, ha="center", va="center", fontsize=6, fontweight="bold",
                      color=text_color_for_bg(celltype_colors[cell_type]), transform=ax_h.transAxes)

        data_axes = np.empty((n_rows, n_cols), dtype=object)
        for ri, cell_type in enumerate(target_cell_types):
            for ci, (transf, week) in enumerate(cols):
                ax = fig.add_subplot(data_gs[ri, 2 + ci], sharey=data_axes[ri, 0] if ci > 0 else None)
                data_axes[ri, ci] = ax

                sub = values_df[
                    (values_df["cell_type"] == cell_type) & (values_df["anchor"] == transf)
                    & (values_df["week"] == week)
                ]
                pts = sub["ratio"].to_numpy()

                if len(pts) == 0:
                    ax.axhline(0, lw=1.4, color="#888888", zorder=2)
                elif len(pts) == 1:
                    ax.scatter([0], pts, s=12, color=transf_colors[transf], zorder=3)
                else:
                    try:
                        parts = ax.violinplot(pts, widths=0.9, bw_method=0.5, showmeans=True, showextrema=False)
                        for body in parts["bodies"]:
                            body.set_facecolor(transf_colors[transf])
                            body.set_edgecolor("black")
                            body.set_linewidth(0.5)
                            body.set_alpha(0.85)
                        parts["cmeans"].set_color("black")
                        parts["cmeans"].set_linewidth(0.9)
                    except Exception:
                        ax.scatter(np.zeros_like(pts), pts, s=8, color=transf_colors[transf], zorder=3)

                ax.axhline(1, ls="--", lw=0.5, color="#888888", alpha=0.7, zorder=0)
                ax.set_xticks([])
                ax.tick_params(axis="y", labelsize=5)
                if ci > 0:
                    ax.tick_params(axis="y", labelleft=False)
                ax.spines[["top", "right", "bottom"]].set_visible(False)

        # fig.suptitle(
        #     f"co-occurrence ratio distribution in the near-transf zone (0-{near_zone_max} µm)",
        #     fontsize=8, fontweight="bold", y=0.90,
        # )
        plt.show()