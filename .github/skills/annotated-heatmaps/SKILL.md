---
name: annotated-heatmaps
description: Create gene (or feature) x condition heatmaps with colored group text labels, colored condition/time strips, significance-star annotations, and a separate legend, in the project's established publication style. Use whenever building a new heatmap for grouped genes/features across samples, batches, weeks, or conditions.
---

# Annotated heatmaps (gene x condition, grouped by category)

This skill captures the exact heatmap style established in
`report_27_supp_3_gsea_natali.ipynb` (Section 5, `plot_realbatch_heatmap`). Apply it
whenever creating a new heatmap of curated genes/features grouped into categories
(pathways, cell types, programs, etc.) across ordered conditions (weeks, batches,
samples).

## 1. Always combine with `publication-grade-scientific-plots`

Use `mpl.rc_context(PUBLICATION_RC)` (see that skill) for every heatmap figure. Do not
duplicate font settings inline except where this skill calls for specific overrides
below.

## 2. Data shape

- Rows: individual genes/features, grouped by category and kept in a fixed, meaningful
  order (do not alphabetize or auto-sort).
- Columns: ordered conditions, optionally nested by a secondary grouping, e.g.
  `(week, batch)` tuples in the order `[(w2,b1), (w2,b2), (w4,b1), (w4,b2), ...]`.
- Value: a single continuous metric (log2FC, z-score, etc.) pivoted to a `gene x column`
  DataFrame with `.reindex(index=gene_order, columns=pd.MultiIndex.from_tuples(col_tuples))`
  to guarantee row/column order regardless of pivot-table default sorting.
- Significance/annotation: a separate `sig` DataFrame of the same shape, mapped to star
  strings (`***`/`**`/`*`/blank), shown as cell annotations instead of (or in addition
  to) numeric values, depending on what the user wants visible.

## 3. Color scale

- Diverging colormap (`RdBu_r` or similar) centered at 0 for signed metrics (log2FC,
  z-scores). Never use a sequential/rainbow colormap for signed data.
- Clip `vmin`/`vmax` to the **5th-95th percentile of `abs(values)`** (not the true
  min/max) for stronger, more readable contrast:
  ```python
  vmax = np.nanpercentile(np.abs(values.to_numpy()), 95)
  ```
  This is the default starting point; tighten further (e.g. 90th) if the user asks for
  more contrast, or use fixed symmetric bounds if reproducibility across figures matters
  more than adaptive contrast.
- Every heatmap cell gets a **black cell border**: `linewidths=0.8, linecolor="black"`
  in `sns.heatmap(...)`.
- When the value is a log2 fold change, the colorbar label is always
  `cbar_kws={"label": r"$\log_2(\mathrm{fold\ change})$", "shrink": 0.5}`. If the figure
  shows several contrasts, name the reference in the caption, not in the colorbar label.
  See the fold-change naming rule in `concise-project-coding`.

## 4. Row/category labeling - colored text, no separate color-bar axes

- Do **not** draw a separate row color-bar `Axes` next to the heatmap (e.g. via
  `plt.subplots(1, 2, ...)` with a narrow side axes). This has repeatedly caused the
  side-bar Axes to render **on top of and hide** the heatmap's y-tick label text due to
  axes z-order/draw-order interactions.
- Instead, color the y-tick label **text itself** to match each gene's category color,
  and put the category-name -> color mapping in a separate figure legend (not inline
  text next to the rows):
  ```python
  ax_heat.set_yticks(np.arange(n_genes) + 0.5)
  ytick_texts = ax_heat.set_yticklabels(gene_order, rotation=0, fontsize=7, fontweight="bold")
  for tick_label, g in zip(ytick_texts, gene_order):
      tick_label.set_color(category_colors[category_of_gene[g]])
  ```
  Always call `set_yticks()` immediately before `set_yticklabels()` (a `FixedLocator` +
  `FixedFormatter` pair) - calling `set_yticklabels()` alone, or reading back
  `get_yticklabels()` afterward, can silently produce empty/misaligned labels.
- Draw horizontal black separator lines between category blocks directly on the heatmap
  axes: `ax_heat.axhline(boundary, color="black", lw=1.2)`.

## 5. Column/condition labeling - stacked color strips, no cluttered x-tick text

- Hide default x-tick labels (`ax_heat.set_xticks([])`) when there is a nested column
  grouping (e.g. batch within week) that is better shown as color strips than as long
  tick-label text.
- Draw one colored strip per grouping level as `Rectangle` patches directly on the same
  `ax_heat` (with `clip_on=False` so they sit outside the 0..n_genes data area), stacked
  vertically above the heatmap, most granular level closest to the data:
  ```python
  # Innermost strip (e.g. batch/repeat): color only.
  for i, color in enumerate(col_batch_colors):
      ax_heat.add_patch(mpatches.Rectangle((i, -0.6), 1, 0.5, color=color, lw=0, clip_on=False))
  # Outer strip (e.g. week): color block + short text label centered above its span.
  for wi, group in enumerate(outer_order):
      x0 = wi * n_inner_per_outer
      ax_heat.add_patch(mpatches.Rectangle((x0, -1.15), n_inner_per_outer, 0.5,
                                            color=outer_palette[group], lw=0, clip_on=False))
      ax_heat.text(x0 + n_inner_per_outer / 2, -1.25, outer_label(group), ha="center",
                   va="bottom", fontsize=7, fontweight="bold")
  ```
- Keep the innermost grouping unlabeled by text when it is shown only for visual
  reference/context (e.g. technical replicates) rather than as the primary reported
  statistic - state this explicitly in a caveat markdown cell above the plotting code if
  the color strip encodes a *descriptive* rather than model-based quantity.

## 6. Legend, not inline text

- Category-name -> color mapping and significance-star meaning both belong in one
  `fig.legend(...)` at the bottom of the figure, not as text scattered around the axes:
  ```python
  legend_handles = (
      [mpatches.Patch(facecolor="none", edgecolor="none", label=lbl)
       for lbl in ["*** padj<0.001", "** padj<0.01", "* padj<0.05", "(blank) not significant"]]
      + [mpatches.Patch(facecolor=category_colors[c], edgecolor="none", label=c) for c in categories]
  )
  fig.legend(handles=legend_handles, loc="lower center", ncol=4, frameon=False,
             fontsize=7, bbox_to_anchor=(0.55, -0.12))
  ```

## 7. Sizing

- Figure width: `max(SINGLE_COL_WIDTH, 0.35 * n_columns + 1.6)` inches - grows with the
  number of columns but never shrinks below the single-column publication width.
- Figure height: `0.19 * n_genes + 0.7` inches - grows with the number of rows.
- Title (`fig.suptitle`) at 9 pt bold; all tick labels, cell-border stars, and legend
  text at 7 pt (6 pt only for the small in-cell significance stars if space is tight).

## 8. Statistical honesty in captions/docstrings

- If a heatmap mixes a rigorously modeled statistic (e.g. sample-level DESeq2 padj) with
  a merely descriptive per-subgroup value (e.g. per-technical-replicate log2FC not
  independently tested), say so explicitly in the function docstring and in a markdown
  cell above the plotting code. Never imply a subgroup-level p-value exists when it was
  not computed.

## 9. Do not save figures/tables unless asked

Per `concise-project-coding`, end the plotting cell with `plt.show()` only; do not add
`savefig`, `to_csv`, or similar export calls unless the user explicitly requests them.
