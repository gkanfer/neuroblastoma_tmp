---
name: publication-grade-scientific-plots
description: Create, revise, or review scientific plots so they are publication-ready and sized for journal submission. Use whenever working on figures, plots, charts, heatmaps, spatial plots, Scanpy/Seaborn/Matplotlib visualizations, or figure-export code in Python or Jupyter notebooks.
---

# Publication-grade scientific plots

Apply these rules whenever creating, modifying, or reviewing a scientific figure.

## 1. Priority of requirements

Use this order of precedence:

1. Explicit requirements in the user's current request.
2. Requirements of the journal named by the user.
3. Existing project-specific plotting conventions.
4. The generic biomedical defaults in this skill.

Never claim that a figure satisfies a named journal unless the journal profile was explicitly selected or its current author instructions were provided or verified.

## 2. Generic biomedical defaults

When no journal is specified, use:

- Typeface: Arial. If Arial is unavailable, use Liberation Sans or another metrically similar sans-serif font and mention the fallback briefly.
- Final-size font: 7 pt for axes, ticks, legends, and annotations.
- Panel labels, only when explicitly requested (see the hard rules below): 9 pt, bold.
- Figure width:
  - Single column: 3.5 inches.
  - Double column: 7.0 inches.
- Figure height: choose the smallest height that preserves readability; normally 2.5–5.5 inches and never make it unnecessarily tall.

### Compactness (applies to every publication figure)

Publication panels are small. Default to compact, tightly packed panels rather than the widths above whenever the figure is a single panel that will be assembled with others.

- Match the size of the sibling panels of the same figure. Check their `figsize` first and reuse it; a mismatched panel looks out of proportion in Illustrator.
- When no sibling exists, start around 1.7–2.5 inches wide and 1.8–2.2 inches tall for a simple panel (few categories, one axes) and only grow it if something is clipped or unreadable.
- Shrink the panel, never the text. Fonts stay at 7 pt; if content does not fit, reduce whitespace, shorten labels, or drop redundant elements.
- Remove wasted space: trim axis limits to the data (`set_xlim`/`set_ylim`), use `constrained_layout=True`, drop titles, drop axis labels whose tick labels are already self-explanatory, and keep margins minimal.
- Bars: keep the gaps between bars small. Set the category spacing just above the bar width (for example `width=0.6` with positions spaced `0.72`) and clip the x-limits to about half a bar beyond the outer bars.
- Legends: place them below the axes in a single horizontal row (`loc="upper center"`, `bbox_to_anchor=(0.5, -0.12)`, `ncol=n_categories`, `frameon=False`) with small `handlelength`, `handletextpad`, and `columnspacing` so they do not widen the panel.
- Raster export: 300 dpi for color or grayscale images.
- Rasterized line art or figures dominated by sharp lines/text: 600 dpi minimum.
- Prefer PDF or SVG for plots, diagrams, and other vector-compatible artwork.
- Use a white background unless the user or journal requests otherwise.
- Use accessible, colorblind-safe colors; do not rely only on red versus green.
- Use consistent fonts, line widths, marker sizes, and color mappings across related figures.

Do not use 12 pt as the default figure font. Use it only when the requested journal, output dimensions, or presentation context specifically calls for it.

### Hard rules: never do these unless explicitly asked

These two override every default in this skill. If you think one should be broken for a specific figure, say so and ask first; do not do it silently.

1. **Never add panel letters** (`a`, `b`, `c`, `d` / `A`, `B`, `C`, `D`) to plots, subplots, or figures. Do not use `ax.text()`, `ax.set_title()`, `fig.text()`, or `enumerate()` over axes to place them. The user adds panel letters when assembling the final figure in Illustrator. Add them only on an explicit request for panel labels.
2. **Never make axis labels bold.** Axis labels (`set_xlabel` / `set_ylabel`), tick labels, and colorbar labels use the regular weight. Do not pass `fontweight="bold"` to those calls and do not set `"axes.labelweight": "bold"` in the rc dict. Journals reserve bold for panel letters, so bold axis text reads as a draft plot. Use bold only when the user asks for it on a specific element.

Related Matplotlib detail: `fontweight` does not apply inside `$...$` mathtext. A label such as `r"$\log_2(\mathrm{fold\ change}) \pm$ s.e.m."` rendered with `fontweight="bold"` comes out half bold (only the text outside the math), which is another reason to leave axis labels at regular weight.

## 3. Journal profiles

When the user names a journal, use the corresponding profile below unless newer or more specific instructions are supplied.

### Nature

- Width: 3.54 inches for single-column or 7.09 inches for double-column figures.
- Maximum height: 6.69 inches.
- Typeface: Arial or Helvetica.
- Final-size text: 5–7 pt; use 7 pt by default for readability.
- Raster images: at least 300 dpi.

### Cell Press

- Final-size text: approximately 6–8 pt; use Arial 7 pt by default.
- Color or grayscale raster figures: at least 300 dpi.
- Black-and-white or line-focused raster figures: use at least 500–600 dpi.
- Prefer editable vector output when possible.

### Science / Science Advances

- Common widths: 3.5, 5.0, or 7.3 inches, depending on layout.
- Use a sans-serif typeface. Prefer the journal-requested face; otherwise use Arial.
- Default text: 7–8 pt; never below the journal minimum.
- Raster images: at least 300 dpi; use higher resolution for line art.

### PLOS

- Width: 2.63–7.5 inches depending on layout.
- Typeface: Arial, Times, or Symbol; prefer Arial for consistency.
- Final-size text: 8–12 pt.
- Raster resolution: 300–600 dpi.

## 4. Jupyter notebook behavior

When code is being executed in a Jupyter notebook or `.ipynb` file:

- Display the figure inline.
- Do not save or export any figure unless the user explicitly asks for a file.
- Do not add `savefig`, `ggsave`, `write_image`, or equivalent export calls by default.
- End the plotting cell with `plt.show()` or the appropriate display call.
- Avoid duplicate rendering: do not both call `plt.show()` and leave a bare figure object as the final expression.
- Use approximately 120–150 dpi for the notebook preview if needed. Publication export dpi is irrelevant until a raster file is actually saved.

When editing a Python script rather than a notebook, return or expose the `fig` and `ax` objects where practical, and save only when requested or when the surrounding script clearly implements an export workflow.

## 5. Matplotlib implementation rules

Prefer local styling with `matplotlib.rc_context` so plotting code does not unexpectedly alter the entire notebook or project.

Use a base style similar to:

```python
import matplotlib as mpl
import matplotlib.pyplot as plt

PUBLICATION_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 130,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
}

with mpl.rc_context(PUBLICATION_RC):
    fig, ax = plt.subplots(
        figsize=(3.5, 2.8),
        constrained_layout=True,
    )
    # Plot data here.
    ax.set_xlabel("X label (unit)")
    ax.set_ylabel("Y label (unit)")
    plt.show()
```

Adapt the width, height, font size, and output settings to the selected journal profile.

Additional Matplotlib rules:

- Create figures at their intended final physical size rather than making oversized figures and relying on later downscaling.
- Keep `legend.fontsize` equal to the body font size (7 pt by default). Do not shrink legend text below the tick and axis label size.
- Prefer `constrained_layout=True` or carefully call `tight_layout()` to prevent clipping.
- When exact final dimensions matter, avoid `bbox_inches="tight"` during export because it can change the final page dimensions.
- For very large scatterplots saved as PDF/SVG, rasterize only the dense point layer while keeping text and axes as vectors, for example `scatter(..., rasterized=True)`.
- Remove top and right spines when appropriate for the plot type, but do not apply this blindly to heatmaps, images, spatial plots, or boxed scientific diagrams.
- Avoid heavy grid lines. Use no grid or subtle major-grid lines only when they improve quantitative reading.
- Use sentence case for axis labels and legends.
- Include units in axis labels when applicable.
- Use scale bars for microscopy or spatial images when physical scale is known.

## 6. Seaborn, Scanpy, and Squidpy

- Seaborn may be used when it is the clearest implementation, but preserve the publication dimensions and typography from this skill.
- Avoid applying a global Seaborn theme that unexpectedly changes unrelated notebook cells. Prefer a local context where possible.
- For Scanpy or Squidpy plotting functions, pass `show=False` when supported, capture the returned axes or figure, and then apply final styling before displaying.
- Do not accept package defaults when they create oversized titles, tiny legends, clipped labels, or inconsistent fonts.
- For categorical biological data, keep category colors stable across plots and samples.
- For dot plots and heatmaps, ensure gene labels remain readable at final size; increase figure height rather than shrinking text below the journal minimum.
- For spatial plots, preserve equal aspect ratio and do not distort tissue geometry.

## 7. Scientific integrity and readability

- Never hide data or alter conclusions for visual appearance.
- Do not truncate axes in a misleading way. Bar plots should normally include zero unless a scientifically justified exception is clearly indicated.
- Show raw observations when practical, especially for small sample sizes.
- Clearly distinguish biological replicates from technical replicates.
- Label error bars or state their meaning in code comments or nearby text when ambiguity is possible.
- Do not use rainbow/jet color maps. Use perceptually uniform sequential or diverging maps appropriate to the data.
- Make missing values visually distinct from low values.
- Avoid 3D effects, shadows, gradients, and decorative chart elements.
- Keep legends concise and place them where they do not cover data.

## 8. Export behavior when explicitly requested

When the user asks to save or export:

- Prefer PDF or SVG for vector-compatible figures.
- Use TIFF or PNG for raster images according to the journal's accepted formats.
- Use 300 dpi for color/grayscale raster images and 600 dpi or more for line art when required.
- Preserve the exact intended physical dimensions.
- Use descriptive filenames without spaces when no naming convention is supplied.
- Do not overwrite an existing figure unless the user asked to replace it.
- Report the output path after saving.

Example:

```python
fig.savefig(
    "figure_1.pdf",
    facecolor="white",
)

fig.savefig(
    "figure_1.png",
    dpi=300,
    facecolor="white",
)
```

## 9. Project label conventions

Lesion classes are always displayed with these short names, never with the raw obs values:

| Value in `obs` | Display label |
| --- | --- |
| `lesion_pc` | PC-Les |
| `lesion_midlobular` | Mid-Les |
| `lesion_pp` | PP-Les |

Use them in legends, tick labels, and annotations. Keep the established lesion colors alongside them: `lesion_pc` `#015F0E`, `lesion_midlobular` `#C6C159`, `lesion_pp` `#8F0683`.

## 10. Final quality check

Before finishing any plotting task, verify:

1. The figure uses the requested or selected journal dimensions.
2. Text is readable at final size and not below the journal minimum.
3. Fonts and categorical colors are consistent.
4. Axis labels, units, and legends are complete, at regular weight, and no panel letters were added unless explicitly requested.
5. Nothing is clipped or overlapping.
6. Color choices remain understandable for color-vision deficiencies and grayscale printing when relevant.
7. Jupyter figures are shown inline and not saved unless requested.
8. Export format and dpi match the artwork type rather than applying 300 dpi blindly to everything.
9. The panel is compact, matches its sibling panels in size, and has no wasted whitespace or oversized gaps.
10. Lesion classes use the PC-Les / Mid-Les / PP-Les labels.
