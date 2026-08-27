import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea

# Matches the PUB_RC style used across report_30_supp_5_distance_celltype_plot_noncrop.ipynb.
PUB_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}


class InteractionMatrixHeatmap:
    """Compute and plot squidpy neighborhood-enrichment z-scores between cell types.

    Parameters
    ----------
    andata : anndata.AnnData
        Annotated data object containing spatial neighbors graph and cell-type
        labels in `.obs`. A filtered copy (cell types present in `colormap`
        only) is stored on `self.andata`; the input object is not modified.
    colormap : dict of str to str, optional
        Mapping from cell-type/category name to hex color, defining both the
        plotting order and label/boundary colors. Defaults to `DEFAULT_COLOR_MAP`,
        a built-in liver cell-type palette (non-hepatocyte, hepatocyte zonation,
        and lesion/transformed zonation groups).
    celltype_col : str, default "cell_type_fine"
        Column in `andata.obs` holding the categorical cell-type/annotation
        labels to correlate.
    title : str, default ""
        Default title displayed above the heatmap; can be overridden per call
        via `plot(title=...)`.
    mask_diagonal : bool, default True
        If True, blank out the diagonal (self-pair) cells, which otherwise
        dominate the shared color scale.
    lower_triangle : bool, default True
        If True, show only the lower triangle of the (symmetric) z-score
        matrix, masking the upper triangle.
    vmax_percentile : float, default 95
        Percentile of absolute z-scores (excluding masked cells) used to set
        symmetric color-scale limits (`vmin`/`vmax`).
    n_perms : int, default 1000
        Number of permutations for the `squidpy.gr.nhood_enrichment`
        permutation test.
    n_jobs : int, optional
        Number of parallel jobs passed to `squidpy.gr.nhood_enrichment`.
    backend : str, default "threading"
        Joblib backend passed to `squidpy.gr.nhood_enrichment`. Squidpy's own
        default ("loky") spawns subprocesses that make the permutation arrays
        read-only, which crashes its numba callback whenever `n_jobs` > 1;
        "threading" shares memory instead and avoids this.

    Examples
    --------
    >>> heatmap = InteractionMatrixHeatmap(
    ...     ad,
    ...     celltype_col="annotation_lesion_class_transformed",
    ...     vmax_percentile=80,
    ...     mask_diagonal=False,
    ...     n_perms=10_000,
    ... )
    >>> heatmap.compute()
    >>> heatmap.cell_counts
    >>> heatmap.plot(show_counts=True)
    """

    DEFAULT_COLOR_MAP = {
        # Neuroblastoma — adrenergic states: blues
        "Adrenergic_Proliferating_NB": "#084594",
        "Adrenergic_Cycling_NB": "#2171B5",
        "Adrenergic_StemLike_NB": "#4292C6",
        "Adrenergic_MYCNHigh_Metabolic_NB": "#6BAED6",
    
        # Neuroblastoma — transitional / mesenchymal states: greens
        "Transitional_NB": "#238B45",
        "Metabolic_NB": "#41AE76",
        "Mesenchymal_NB": "#74C476",
    
        # Myeloid — reds
        "M2_TAM": "#B2182B",
        "MonocyteDerived_TissueResident_Macrophage": "#FF5F59",
    
        # Lymphoid — orange
        "CD8_Cytotoxic_TCell_PartiallyExhausted": "#F28E2B",
    
        # Stromal / adrenal — purples
        "Schwann_CAF": "#6A1B9A",
        "Endothelial": "#9C6BC4",
        "Adrenocortical_Steroidogenic": "#CE93D8",
    }

    def __init__(self, andata, colormap=None, celltype_col="cell_type_fine", title="", mask_diagonal=True,
                 lower_triangle=True, vmax_percentile=95, n_perms=1000, n_jobs=None, backend="threading"):
        self.celltype_col = celltype_col
        self.title = title
        self.mask_diagonal = mask_diagonal
        self.lower_triangle = lower_triangle
        self.vmax_percentile = vmax_percentile
        self.n_perms = n_perms
        self.n_jobs = n_jobs
        self.backend = backend

        color_map = self.DEFAULT_COLOR_MAP if colormap is None else colormap
        self.color_map = color_map
        self.non_hep_keys, self.hep_keys, self.transf_keys = (
            list(color_map)[:14], list(color_map)[14:20], list(color_map)[20:]
        )

        present_celltypes = andata.obs[celltype_col].dropna().astype(str).unique()
        self.ordered_celltypes = [ct for ct in color_map.keys() if ct in present_celltypes]
        self.andata = andata[andata.obs[celltype_col].astype(str).isin(self.ordered_celltypes)].copy()
        self.andata.obs[celltype_col] = pd.Categorical(
            self.andata.obs[celltype_col].astype(str), categories=self.ordered_celltypes, ordered=True,
        )
        self.andata.uns[f"{celltype_col}_colors"] = [color_map.get(ct, "#bdbdbd") for ct in self.ordered_celltypes]

        self.zscore = None
        self.mask = None
        self.vmax = None

    @property
    def cell_counts(self):
        """pandas.Series: Cell counts for `ordered_celltypes` only, in plot order (0 if a type is missing)."""
        counts = self.andata.obs[self.celltype_col].value_counts()
        return counts.reindex(self.ordered_celltypes, fill_value=0)

    def compute(self, stat=True):
        """Run squidpy neighborhood-enrichment (or interaction_matrix) on the filtered data.

        Parameters
        ----------
        stat : bool, default True
            If True, run `squidpy.gr.nhood_enrichment` and store the z-score
            matrix/mask/vmax used by `plot()`. If False, only compute
            `squidpy.gr.interaction_matrix` (normalized) and return.
        """
        if not stat:
            sq.gr.interaction_matrix(self.andata, cluster_key=self.celltype_col, normalized=True)
            return

        sq.gr.nhood_enrichment(self.andata, cluster_key=self.celltype_col, show_progress_bar=True,
                                n_perms=self.n_perms, n_jobs=self.n_jobs, backend=self.backend, seed=1980)
        zscore = self.andata.uns[f"{self.celltype_col}_nhood_enrichment"]["zscore"].astype(float).copy()
        n = len(self.ordered_celltypes)

        # The matrix is symmetric - the upper triangle is redundant, so mask it (lower_triangle),
        # and mask the diagonal (mask_diagonal) since self-pairs get much larger z-scores than
        # cross-type pairs and would otherwise dominate a shared linear color scale.
        mask = np.triu(np.ones((n, n), dtype=bool), k=1) if self.lower_triangle else np.zeros((n, n), dtype=bool)
        if self.mask_diagonal:
            np.fill_diagonal(mask, True)

        self.zscore = zscore
        self.mask = mask
        self.vmax = np.nanpercentile(np.abs(zscore[~mask]), self.vmax_percentile)

    def plot(self, title=None, show_counts=False):
        """Draw the z-score heatmap prepared by `compute(stat=True)`.

        Parameters
        ----------
        title : str, optional
            Title displayed above the heatmap. Defaults to the `title` given
            to the constructor.
        show_counts : bool, default False
            If True, append each cell type's cell count to its y-axis label
            only, e.g. "T cells (n=64)", with the count in black.

        Returns
        -------
        None
            Displays the heatmap via `plt.show()`.
        """
        if self.zscore is None:
            raise RuntimeError("Call compute(stat=True) before plot().")
        title = self.title if title is None else title

        zscore, mask, vmax = self.zscore, self.mask, self.vmax
        ordered_celltypes, color_map = self.ordered_celltypes, self.color_map
        n = len(ordered_celltypes)

        with mpl.rc_context(PUB_RC):
            fig, ax = plt.subplots(figsize=(max(3.5, 0.35 * n + 1.6), 0.19 * n + 3.0))
            sns.heatmap(zscore, mask=mask, xticklabels=ordered_celltypes, yticklabels=ordered_celltypes,
                        cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                        linewidths=0.05, linecolor="black",
                        cbar_kws={"label": "Neighborhood Enrichment Score", "shrink": 0.5, "aspect": 15}, ax=ax)

            # Masked cells (upper triangle / diagonal) should show nothing - blank them out
            # completely instead of relying on "no fill", which still leaves faint grid lines.
            for i, j in zip(*np.where(mask)):
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="white", edgecolor="none", zorder=10))

            cbar = ax.collections[0].colorbar
            cbar.outline.set_edgecolor("black")
            cbar.outline.set_linewidth(1)

            for tick_label, ct in zip(ax.get_xticklabels(), ordered_celltypes):
                tick_label.set_color(color_map[ct])
                tick_label.set_fontweight("bold")

            if show_counts:
                # A single Text can't mix colors, so pack "<name>" (colored) + "(n=<count>)" (black)
                # as separate boxes and hide the plain y-tick text seaborn already drew.
                counts = self.cell_counts
                for tick_label in ax.get_yticklabels():
                    tick_label.set_visible(False)
                for i, ct in enumerate(ordered_celltypes):
                    name_box = TextArea(ct, textprops=dict(color=color_map[ct], fontweight="bold", fontsize=7))
                    count_box = TextArea(f" (n={counts[ct]})", textprops=dict(color="black", fontsize=7))
                    packed = HPacker(children=[name_box, count_box], align="center", pad=0, sep=1)
                    ab = AnnotationBbox(packed, (0, i + 0.5), xycoords=("axes fraction", "data"),
                                         box_alignment=(1, 0.5), xybox=(-2, 0), boxcoords="offset points",
                                         frameon=False, annotation_clip=False)
                    ax.add_artist(ab)
            else:
                for tick_label, ct in zip(ax.get_yticklabels(), ordered_celltypes):
                    tick_label.set_color(color_map[ct])
                    tick_label.set_fontweight("bold")

            boundary1 = sum(ct in self.non_hep_keys for ct in ordered_celltypes)
            boundary2 = boundary1 + sum(ct in self.hep_keys for ct in ordered_celltypes)
            for boundary in (boundary1, boundary2):
                if 0 < boundary < n:
                    ax.axhline(boundary, color="black", lw=1.2)
                    ax.axvline(boundary, color="black", lw=1.2)

            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color("black")
                spine.set_linewidth(0.01)

            ax.set_title(title, fontsize=9, fontweight="bold")
            plt.show()


def plotInteraction_matrix(andata,colormap = None,celltype_col = "cell_type_fine",title = "",stat = True,mask_diagonal = True,lower_triangle = True,vmax_percentile = 95,n_perms = 1000,n_jobs = None,backend = "threading",show_counts = False):
    """Plot a lower-triangle heatmap of squidpy neighborhood-enrichment z-scores between cell types.

    Thin backward-compatible wrapper around `InteractionMatrixHeatmap`; see that
    class to compute stats once and reuse `cell_counts`/`plot` separately.

    Parameters
    ----------
    andata : anndata.AnnData
        Annotated data object containing spatial neighbors graph and cell-type
        labels in `.obs`. A filtered copy (cell types present in `colormap`
        only) is used internally; the input object is not modified.
    colormap : dict of str to str, optional
        Mapping from cell-type/category name to hex color, defining both the
        plotting order and label/boundary colors. Defaults to a built-in liver
        cell-type palette (non-hepatocyte, hepatocyte zonation, and
        lesion/transformed zonation groups).
    celltype_col : str, default "cell_type_fine"
        Column in `andata.obs` holding the categorical cell-type/annotation
        labels to correlate.
    title : str, default ""
        Title displayed above the heatmap.
    stat : bool, default True
        If True, run `squidpy.gr.nhood_enrichment` and plot the z-score
        heatmap. If False, only compute `squidpy.gr.interaction_matrix`
        (normalized) and return without plotting.
    mask_diagonal : bool, default True
        If True, blank out the diagonal (self-pair) cells, which otherwise
        dominate the shared color scale.
    lower_triangle : bool, default True
        If True, show only the lower triangle of the (symmetric) z-score
        matrix, masking the upper triangle.
    vmax_percentile : float, default 95
        Percentile of absolute z-scores (excluding masked cells) used to set
        symmetric color-scale limits (`vmin`/`vmax`).
    n_perms : int, default 1000
        Number of permutations for the `squidpy.gr.nhood_enrichment`
        permutation test.
    n_jobs : int, optional
        Number of parallel jobs passed to `squidpy.gr.nhood_enrichment`.
    backend : str, default "threading"
        Joblib backend passed to `squidpy.gr.nhood_enrichment`; "threading"
        avoids a squidpy/numba crash that occurs with its own "loky" default
        when `n_jobs` > 1 (see `InteractionMatrixHeatmap`).
    show_counts : bool, default False
        If True, append each cell type's cell count to its axis label,
        e.g. "T cells (n=64)".

    Returns
    -------
    None
        Displays the heatmap via `plt.show()` (when `stat=True`) or returns
        after updating `andata.uns` in place (when `stat=False`). Nothing is
        returned in either case.

    Examples
    --------
    >>> plotInteraction_matrix(
    ...     ad,
    ...     celltype_col="annotation_lesion_class_transformed",
    ...     vmax_percentile=80,
    ...     mask_diagonal=False,
    ...     n_perms=10_000,
    ... )
    """
    heatmap = InteractionMatrixHeatmap(andata, colormap=colormap, celltype_col=celltype_col, title=title,
                                        mask_diagonal=mask_diagonal, lower_triangle=lower_triangle,
                                        vmax_percentile=vmax_percentile, n_perms=n_perms, n_jobs=n_jobs,
                                        backend=backend)
    heatmap.compute(stat=stat)
    if not stat:
        return
    heatmap.plot(title=title, show_counts=show_counts)