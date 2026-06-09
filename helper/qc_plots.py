"""QC histogram grid — one row per sample, four columns:
  total_counts, n_genes_by_counts, log(total_counts), log(n_genes_by_counts).
Bars are grey, KDE overlay is black.

Sizing and fonts follow Nature publication requirements:
  - Font: Arial, minimum 7 pt labels / 8 pt axis titles
  - Panel width ~45 mm (4 panels = 180 mm ≈ Nature double-column 183 mm)
  - Panel height ~38 mm
  - Line widths ≥ 0.5 pt
  - Export at 300 dpi (colour) / 600 dpi (line art)
  - PDF font type 42 (TrueType outlines)
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ---- Nature-compliant rcParams applied locally per figure -------------------
_NATURE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,            # base (Nature minimum for figure text)
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

# Panel dimensions in inches (Nature double-column = 183 mm = 7.2 in)
# 4 panels across → each ~1.77 in (45 mm); height ~1.5 in (38 mm).
_PANEL_W = 45 / 25.4   # mm → in
_PANEL_H = 38 / 25.4


def plot_qc_histograms(
    adata,
    sample_col="sample_key",
    counts_col="total_counts",
    genes_col="n_genes_by_counts",
    bins=60,
    figsize_per_panel=(_PANEL_W, _PANEL_H),
    save_path=None,
    dpi=300,
):
    """Create a multi-row histogram grid (one row per sample) with 4 columns:
    raw counts, raw genes, log(counts), log(genes).

    Parameters
    ----------
    adata : AnnData
        Must have `counts_col` and `genes_col` in `.obs`.
    sample_col : str
        Column used to split rows (default "sample_key").
    counts_col, genes_col : str
        obs columns for total counts and detected genes.
    bins : int
        Number of histogram bins.
    figsize_per_panel : tuple
        (width, height) in inches per individual subplot.
    save_path : str or Path, optional
        If given, saves the figure (PNG + PDF).
    dpi : int
        Resolution for saved PNG (PDF is vector).

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    obs = adata.obs.copy()
    samples = sorted(obs[sample_col].unique())
    n_samples = len(samples)

    ncols = 2
    fw, fh = figsize_per_panel

    with mpl.rc_context(_NATURE_RC):
        fig, axes = plt.subplots(
            n_samples, ncols,
            figsize=(fw * ncols, fh * n_samples),
            constrained_layout=True,
        )
        if n_samples == 1:
            axes = axes[np.newaxis, :]

        col_labels = [
            (counts_col, False, "Reads Count\\Cells"),
            (genes_col,  False, "Genes Count\\Cells"),
            # (counts_col, True,  "Reads Count\\Cells (log)"),
            # (genes_col,  True,  "Genes Count\\Cells (log)"),
        ]

        for row_idx, sample in enumerate(samples):
            sub = obs[obs[sample_col] == sample]
            for col_idx, (col, use_log, xlabel) in enumerate(col_labels):
                ax = axes[row_idx, col_idx]
                vals = sub[col].values.astype(float)
                if use_log:
                    vals = np.log1p(vals)

                sns.histplot(
                    vals, bins=bins, ax=ax,
                    color="#B0B0B0", edgecolor="white", linewidth=0.3,
                    stat="count", kde=True,
                    line_kws={"color": "black", "linewidth": 0.8},
                )
                ax.set_xlabel(xlabel)
                if col_idx == 0:
                    ax.set_ylabel(sample, fontsize=8, fontweight="bold")
                else:
                    ax.set_ylabel("")

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(save_path), dpi=dpi, bbox_inches="tight")
        # Also save a PDF for vector figures (Nature prefers vector)
        pdf_path = save_path.with_suffix(".pdf")
        fig.savefig(str(pdf_path), bbox_inches="tight")

    return fig


# ---------------------------------------------------------------------------
# Box-plot panel: reads/cell, genes/cell, cells/gene
# ---------------------------------------------------------------------------

def plot_qc_boxplots(
    adata,
    sample_col="sample_key",
    counts_col="total_counts",
    genes_col="n_genes_by_counts",
    save_path=None,
    dpi=300,
    order_num = True,
    limit_y_label = True,
):
    """Three-panel box plot (one per metric) grouped by sample.

    Box shows 5th, 25th, 50th, 75th, 95th percentiles.
    Black frame, grey fill — same publication style as the histogram grid.

    Metrics:
      1. Reads (UMI) count per cell
      2. Genes detected per cell
      3. Cells per gene (number of cells expressing each gene)

    Parameters
    ----------
    adata : AnnData
    sample_col : str
    counts_col, genes_col : str
    save_path : str or Path, optional
    dpi : int

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import matplotlib as mpl

    obs = adata.obs.copy()
    if order_num:
        order_int = sorted(obs[sample_col].unique().astype(int).tolist())
        samples = [str(i) for i in order_int]
    else:
        samples = sorted(obs[sample_col].unique())

    # Compute cells-per-gene per sample
    from scipy.sparse import issparse
    X = adata.layers.get("counts", adata.X)
    cells_per_gene_data = []
    sample_arr = obs[sample_col].values
    for s in samples:
        mask = (sample_arr == s)
        Xs = X[mask]
        cpg = np.asarray((Xs > 0).sum(axis=0)).ravel()
        cells_per_gene_data.append(pd.DataFrame({
            "value": cpg,
            sample_col: s,
            "metric": "Cells\\Gene",
        }))
    cpg_df = pd.concat(cells_per_gene_data, ignore_index=True)

    # Build long-form dataframes for the other two metrics
    reads_df = obs[[sample_col, counts_col]].rename(columns={counts_col: "value"})
    reads_df["metric"] = "Reads Count\\Cell"
    genes_df = obs[[sample_col, genes_col]].rename(columns={genes_col: "value"})
    genes_df["metric"] = "Genes Count\\Cell"

    metrics_order = ["Reads Count\\Cell", "Genes Count\\Cell", "Cells\\Gene"]
    long = pd.concat([reads_df, genes_df, cpg_df], ignore_index=True)

    # Nature-style sizing: 3 panels across double-column (183 mm)
    panel_w = 60 / 25.4   # 60 mm each → 180 mm total
    panel_h = 55 / 25.4   # 55 mm height

    with mpl.rc_context(_NATURE_RC):
        fig, axes = plt.subplots(
            3, 1, figsize=(panel_w * 3, panel_h),
            constrained_layout=True,
        )
        for ax, metric in zip(axes, metrics_order):
            sub = long[long["metric"] == metric]

            # Custom percentile box: whis = [5, 95]
            bp = ax.boxplot(
                [sub.loc[sub[sample_col] == s, "value"].values for s in samples],
                labels=samples,
                whis=[5, 95],       # whiskers at 5th and 95th percentile
                showfliers=False,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=0.8),
                boxprops=dict(facecolor="#D0D0D0", edgecolor="black", linewidth=0.6),
                whiskerprops=dict(color="black", linewidth=0.6),
                capprops=dict(color="black", linewidth=0.6),
            )
            ax.set_title(metric, fontsize=6, fontweight="bold")
            ax.set_ylabel("")
            ax.tick_params(axis="x", rotation=0)
            
            if limit_y_label:
                values = sub["value"].dropna()
                #ax.set_ylim(0, values.max())
                ax.set_yticks([0, values.median()]) 
            else:
                ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=10))
                ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
                
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(save_path), dpi=dpi, bbox_inches="tight")
        fig.savefig(str(save_path.with_suffix(".pdf")), bbox_inches="tight")

    return fig
