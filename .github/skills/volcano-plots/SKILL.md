---
name: volcano-plots
description: Create differential-expression volcano plots (log2 fold change vs -log10 adjusted P) in the project's established style, with non-overlapping gene labels. Use whenever building a new volcano plot from DESeq2/GSEA-style results, unless the user explicitly asks for a different style.
---

# Volcano plots (log2FC vs -log10 adjusted P)

This skill captures the volcano-plot style established in
`report_30_supp_4_gseapy_deg_lesion_pc_vs_lesion_pp.ipynb`. Apply it by default any
time a new volcano plot is requested from a differential-expression results table
(columns like `log2FoldChange`, `padj`, `gene`), unless the user explicitly asks for a
different look.

## 1. Always combine with `publication-grade-scientific-plots`

Use the generic biomedical defaults (fonts, dpi, inline display, no autosave) from that
skill for every volcano figure.

## 2. Data prep

```python
y = -np.log10(df["padj"].clip(lower=1e-300))  # clip avoids log10(0) = inf

sig_up = (df["padj"] < padj_threshold) & (df["log2FoldChange"] > log2fc_threshold)
sig_down = (df["padj"] < padj_threshold) & (df["log2FoldChange"] < -log2fc_threshold)
```

- Default thresholds: `padj_threshold = 0.05`, `log2fc_threshold = 1`. Use the user's
  thresholds if given instead.

## 3. Colors and layout

- Three point classes, always in this order (draw non-significant first so significant
  points sit on top):
  - Non-significant: `OTHER_COLOR = "#BDBDBD"`, `s=8`.
  - Down-significant: `DOWN_COLOR = "#5DADE2"`, `s=10`.
  - Up-significant: `UP_COLOR = "#C0392B"`, `s=10`.
- Dashed grey threshold guides (`color="grey", ls="--", lw=0.6`): vertical lines at
  `+log2fc_threshold` and `-log2fc_threshold`, horizontal line at
  `-np.log10(padj_threshold)`.
- Figure size: `figsize=(5, 4)` as the default starting point for a single volcano panel.

## 4. Axis labels and title

Always use these exact axis labels unless the user asks otherwise:

```python
ax.set_xlabel(r"$\log_2(\mathrm{fold\ change})$")
ax.set_ylabel(r"$-\log_{10}(\mathrm{adjusted}\ P)$")
```

Never write the x-label as plain `"log2 fold change"`, `"log2FC"`, or
`r"$\Delta z$ (log2FC)"`. See the fold-change naming rule in `concise-project-coding`.

Title format: `"<comparison context>: <group A> vs <group B>"`, e.g.
`"Week 2: lesion_pc vs lesion_pp"`.

## 5. Gene labels - always label, never overlapping

- Always label genes on the plot by default (do not skip labeling unless the user says
  otherwise), even if the request doesn't mention it explicitly.
- Select labels from the significant points only (`sig_up | sig_down`), ranked by
  significance: `top_idx = y[sig_up | sig_down].nlargest(n_labels).index`, with
  `n_labels` defaulting to 8-10 (increase modestly if the user asks for more genes
  labeled, but keep it small enough to stay readable).
- Use `adjustText.adjust_text` to resolve overlaps - this is required, not optional:
  ```python
  from adjustText import adjust_text

  texts = [
      ax.text(df.loc[i, "log2FoldChange"], y.loc[i], df.loc[i, "gene"], fontsize=7)
      for i in top_idx
  ]
  adjust_text(texts, ax=ax)
  ```
- Never manually offset label positions with fixed `xytext` shifts as a substitute for
  `adjust_text` - it doesn't generalize across different point densities.

## 6. Do not save figures unless asked

End the plotting cell with `plt.tight_layout()` then `plt.show()` only; do not add
`savefig` or similar export calls unless the user explicitly requests them.
