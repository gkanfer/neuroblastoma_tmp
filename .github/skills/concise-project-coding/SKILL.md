---
name: concise-project-coding
description: Write or revise Python and Jupyter notebook code in the project's established style: concise, direct, readable, minimally dependent, and safe with files and environments. Use for coding, analysis, notebook editing, refactoring, debugging, and data-processing tasks.
---

# Concise project coding

Apply these rules whenever creating, editing, or reviewing Python code or Jupyter notebooks.

## 1. Priority

Follow this order:

1. The user's explicit request in the current prompt.
2. The conventions in the current file and nearby project code.
3. The rules in this skill.

Make the smallest change that completes the requested task. Do not refactor unrelated code.

## 2. Keep code short, clear, and direct

- Use the shortest clear solution, not code golf.
- Prefer a simple linear sequence of operations over unnecessary abstraction.
- Reuse existing variables, objects, functions, imports, and project conventions.
- Use descriptive `snake_case` names consistent with the current file.
- Prefer list or dictionary comprehensions when they remain easy to read.
- Avoid deeply nested loops, repeated blocks, hidden side effects, and long chains of unrelated operations.
- Group related steps together so the code is concise without becoming spaghetti code.
- Add only short comments that explain intent or a non-obvious decision.
- Do not add large docstrings, extensive type architecture, classes, tests, or validation unless requested or already required by the surrounding code.

## 3. Dependencies and environments

- Prefer Python's standard library and common scientific packages already used by the project, such as NumPy, pandas, Matplotlib, Seaborn, SciPy, Scanpy, AnnData, and pathlib.
- Do not introduce a rare or specialized package when the task can be completed clearly with packages already available in the project.
- Check existing imports and the active notebook kernel before suggesting another dependency.
- Do not run or propose `pip install`, `conda install`, `mamba install`, environment creation, kernel creation, or environment replacement unless the user explicitly asks.
- Use the currently active environment or Jupyter kernel.

## 4. Functions, classes, and helper files

- Do not create a new helper module, Python file, package, class, or utility directory unless the user explicitly asks.
- In a notebook, keep task-specific code in the current notebook.
- Define a small local function only when it clearly removes meaningful repetition, improves readability, or the user asks for reusable code.
- Do not split a short analysis into several functions merely to appear modular.
- When a suitable project function already exists, reuse it rather than creating a duplicate.

## 5. Jupyter notebook behavior

When working in a `.ipynb` file:

- Add a short, informative Markdown heading before a new analysis section, for example `## Compare cell-type proportions`.
- Keep the heading specific and brief; do not add decorative or overly long introductions.
- Run code in the current notebook and current kernel.
- Display tables, summaries, and plots inline.
- Do not create a separate script or helper file unless explicitly requested.
- Do not save figures, tables, models, intermediate objects, or notebook copies unless explicitly requested.
- Avoid duplicate display calls or unnecessary diagnostic output.

## 6. Protect input files and source data

- Treat input files as read-only unless the user explicitly asks to modify or replace them.
- Do not overwrite the source notebook, AnnData file, CSV, pickle, Zarr store, image, model, configuration, or other input artifact.
- Do not add write operations such as `write_h5ad`, `write_zarr`, `to_csv`, `to_excel`, `pickle.dump`, `joblib.dump`, or `savefig` unless saving is requested.
- Do not silently change an input path into an output path.
- When the user requests an output file, preserve the original and use a clear new filename unless overwrite behavior is explicitly requested.
- Make in-memory changes only when needed for the requested analysis; avoid irreversible mutation of the main object when a small derived object or column is sufficient.

## 7. Edge cases and defensive code

- Do not add broad edge-case handling, fallback systems, logging frameworks, or defensive wrappers by default.
- Handle an edge case when the user requests it, when the current data shows it, or when failure is obvious and likely, such as an empty selection, missing required column, or a one-panel axes-shape issue.
- Keep necessary checks local and brief.
- Do not suppress exceptions broadly or use bare `except` blocks.

## 8. Match the user's established style

Prefer the style of the current open file. When local context is limited, use these conventions observed in the user's GitHub projects:

- Standard scientific Python imports near the top of the file or notebook section.
- Direct pandas, NumPy, Matplotlib, Seaborn, Scanpy, and SpatialData operations.
- Clear `snake_case` variables and functions.
- Explicit lists and dictionaries for annotations, colors, paths, and parameters.
- Straightforward loops for samples, conditions, panels, or cell types.
- Small comprehensions for simple mappings and filtering.
- Visible intermediate variables when they make biological or analytical steps easier to inspect.
- Practical, readable code rather than heavy software-engineering patterns for one-off analyses.

Do not imitate old mistakes or inconsistent formatting from a repository. Match the useful local conventions while keeping new code readable and correct.

## 9. Fold-change naming and labeling

One name per concept, everywhere in this project:

- **Model output**: keep the PyDESeq2 column name `log2FoldChange` unchanged. When one
  table holds several contrasts, suffix it: `log2FC_<case>_vs_<reference>`, with matching
  `padj_*` / `stat_*`.
- **Not from a model**: any descriptive ratio (e.g. library-size-normalized group means)
  must be suffixed `_descriptive`, e.g. `log2FC_descriptive`, and flagged as descriptive
  in the markdown cell and the function docstring.
- **Thresholds**: `log2fc_threshold`, `padj_threshold`. Function parameters: `lfc_col`,
  `lfc_min`.
- **Axis, colorbar, and legend text**: always `r"$\log_2(\mathrm{fold\ change})$"`.
  Aggregates across genes: `r"Mean $\log_2(\mathrm{fold\ change}) \pm$ s.e.m."`.
  Never `"log2 fold change"`, `"log2FC"`, or `r"$\Delta z$ (log2FC)"`. Use mathtext, not
  the Unicode subscript `₂`, which is missing from some Arial builds and silently
  falls back to another font in the PDF.
- **Sign convention**: positive means higher in the case group named first in the
  contrast. State this in the figure caption, not on the axis.

## 10. Response and edit scope

- Return or edit only the code needed for the requested task.
- Avoid unnecessary alternatives, optional frameworks, and speculative additions.
- Preserve existing names and structure when practical.
- Briefly state any assumption that materially affects the result.
- When asked to edit code, do not change unrelated cells or files.
