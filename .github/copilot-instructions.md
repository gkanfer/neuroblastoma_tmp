# Repository-wide Copilot instructions

For all Python, data-analysis, and Jupyter notebook work, use the `concise-project-coding` skill. Keep code short but readable, follow the style of the current file, prefer existing and common scientific packages, and make only the changes needed for the request.

Never install packages, create or replace an environment, create helper files, overwrite input files, or save outputs unless the user explicitly asks. In Jupyter notebooks, add a short informative Markdown heading for each new analysis section and display results inline.

Whenever creating, editing, or reviewing scientific plots or figure-generation code, also use the `publication-grade-scientific-plots` skill. Treat publication dimensions, typography, accessibility, scientific integrity, and Jupyter inline-display behavior as required rather than optional.

Whenever creating a gene/feature x condition heatmap grouped by category (pathways, cell types, batches, weeks, etc.), also use the `annotated-heatmaps` skill for color scale, category-colored labels, condition color strips, and legend conventions.

Whenever creating a volcano plot (log2 fold change vs adjusted P) from differential-expression results, also use the `volcano-plots` skill for colors, thresholds, axis labels, and non-overlapping gene labeling.

Whenever asked to create, build, or install a new Mamba/Conda environment, use the `mamba-environment` skill.
