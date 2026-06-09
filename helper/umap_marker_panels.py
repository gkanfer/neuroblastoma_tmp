from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import to_hex


def normalize_cluster_id(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0", text):
        text = text[:-2]
    return text


def select_expression_layer(adata, preferred_layer="log"):
    if preferred_layer in getattr(adata, "layers", {}):
        print(f"Expression layer used: adata.layers['{preferred_layer}']")
        return preferred_layer
    print("Expression layer used: adata.X")
    return None


def find_umap_basis(adata, preferred=("X_umap_60neig", "X_umap", "umap")):
    obsm_keys = set(adata.obsm.keys())
    for basis in preferred:
        if basis in obsm_keys:
            return basis
        if f"X_{basis}" in obsm_keys:
            return f"X_{basis}"
    umap_like = [key for key in adata.obsm.keys() if "umap" in key.lower()]
    if umap_like:
        return umap_like[0]
    raise KeyError(f"No UMAP coordinates found in adata.obsm. Available keys: {list(adata.obsm.keys())}")


def get_umap_coordinates(adata, basis=None):
    basis = basis or find_umap_basis(adata)
    key = basis if basis in adata.obsm else f"X_{basis}"
    if key not in adata.obsm:
        raise KeyError(f"UMAP basis '{basis}' was not found in adata.obsm")
    coords = np.asarray(adata.obsm[key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(f"adata.obsm['{key}'] must be a two-column coordinate matrix")
    return coords[:, :2], key


def match_gene_to_var_names(adata, gene):
    gene = str(gene).strip()
    if gene in adata.var_names:
        return gene
    lookup = {str(var).upper(): str(var) for var in adata.var_names}
    return lookup.get(gene.upper())


def get_expression_vector(adata, gene, layer=None):
    matched_gene = match_gene_to_var_names(adata, gene)
    if matched_gene is None:
        raise KeyError(f"Gene '{gene}' is not present in adata.var_names")
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' is not present in adata.layers")
        matrix = adata[:, matched_gene].layers[layer]
    else:
        matrix = adata[:, matched_gene].X
    if sp.issparse(matrix):
        values = matrix.toarray().ravel()
    else:
        values = np.asarray(matrix).ravel()
    return pd.Series(values, index=adata.obs_names, name=matched_gene), matched_gene


def make_binned_expression_column(adata, gene, n_bins=20, layer=None, prefix="__expr_bin"):
    expr, matched_gene = get_expression_vector(adata, gene, layer=layer)
    positive = expr[expr > 0].copy()
    column = f"{prefix}_{matched_gene}"

    if positive.empty:
        adata.obs[column] = pd.Categorical(pd.Series(pd.NA, index=adata.obs_names, dtype="object"))
        return column, matched_gene, {}, expr

    n_actual_bins = min(int(n_bins), positive.nunique())
    if n_actual_bins <= 1:
        labels = ["bin_1"]
        binned_positive = pd.Series("bin_1", index=positive.index, dtype="object")
    else:
        cats = pd.qcut(positive, q=n_actual_bins, duplicates="drop")
        n_actual_bins = len(cats.cat.categories)
        labels = [f"bin_{idx + 1}" for idx in range(n_actual_bins)]
        binned_positive = cats.cat.rename_categories(labels).astype(str)

    binned = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
    binned.loc[binned_positive.index] = binned_positive
    adata.obs[column] = pd.Categorical(binned, categories=labels, ordered=True)
    palette = dict(zip(labels, [to_hex(plt.cm.Reds(v)) for v in np.linspace(0.35, 0.95, len(labels))]))
    return column, matched_gene, palette, expr


def _best_column(columns, include_any, include_all=()):
    scored = []
    for column in columns:
        lower = str(column).lower()
        if include_all and not all(token in lower for token in include_all):
            continue
        score = sum(token in lower for token in include_any) + sum(token in lower for token in include_all)
        if score:
            scored.append((score, column))
    if not scored:
        return None
    return sorted(scored, key=lambda item: (-item[0], str(item[1])))[0][1]


def _detect_marker_pair_column(df):
    name_based = _best_column(
        df.columns,
        include_any=("marker", "gene", "support", "pair"),
        include_all=(),
    )
    plus_counts = {
        column: df[column].astype(str).str.contains(r"\+", regex=True, na=False).sum()
        for column in df.columns
    }
    value_based = max(plus_counts, key=plus_counts.get) if plus_counts else None
    if value_based is not None and plus_counts[value_based] > 0:
        if name_based is None or plus_counts[value_based] >= plus_counts.get(name_based, 0):
            return value_based
    return name_based


def _read_marker_sheet(excel_path):
    sheets = pd.read_excel(excel_path, sheet_name=None)
    best_sheet = None
    best_score = -1
    for sheet_name, df in sheets.items():
        marker_col = _detect_marker_pair_column(df)
        score = 0 if marker_col is None else df[marker_col].astype(str).str.contains(r"\+", regex=True, na=False).sum()
        if score > best_score:
            best_sheet = (sheet_name, df)
            best_score = score
    if best_sheet is None:
        raise ValueError(f"No sheets found in {excel_path}")
    return best_sheet


def parse_marker_pair(value):
    if pd.isna(value):
        return []
    markers = [part.strip() for part in re.split(r"\s*\+\s*", str(value))]
    return [marker for marker in markers if marker]


def load_marker_pairs_from_excel(excel_path, clusters_to_focus, max_pairs_per_cluster=3):
    excel_path = Path(excel_path)
    if not excel_path.exists():
        raise FileNotFoundError(f"Marker-support workbook was not found: {excel_path}")

    sheet_name, df = _read_marker_sheet(excel_path)
    cluster_col = _best_column(df.columns, include_any=("cluster",), include_all=())
    annotation_col = _best_column(df.columns, include_any=("annotation", "cell type", "cell_type", "proposed", "label"))
    marker_col = _detect_marker_pair_column(df)

    missing = [
        name
        for name, column in (("cluster", cluster_col), ("annotation", annotation_col), ("supporting marker pair", marker_col))
        if column is None
    ]
    if missing:
        raise ValueError(
            f"Could not detect required column(s): {missing}. "
            f"Columns in sheet '{sheet_name}': {list(df.columns)}"
        )

    focus = [normalize_cluster_id(cluster) for cluster in clusters_to_focus]
    work = df.copy()
    work["__cluster__"] = work[cluster_col].map(normalize_cluster_id)
    work = work[work["__cluster__"].isin(focus)].copy()

    marker_pairs = {}
    for cluster_id in focus:
        sub = work[work["__cluster__"] == cluster_id]
        pairs = []
        annotation = None
        for _, row in sub.iterrows():
            markers = parse_marker_pair(row[marker_col])
            if not markers:
                continue
            if annotation is None and pd.notna(row[annotation_col]):
                annotation = str(row[annotation_col]).strip()
            pairs.append({"markers": markers, "source_pair": str(row[marker_col]).strip()})
            if len(pairs) >= max_pairs_per_cluster:
                break
        marker_pairs[cluster_id] = {
            "annotation": annotation or "Unannotated",
            "pairs": pairs,
        }

    print(f"Loaded marker-support workbook: {excel_path}")
    print(f"Sheet used: {sheet_name}")
    print(f"Detected columns: cluster='{cluster_col}', annotation='{annotation_col}', marker_pair='{marker_col}'")
    return marker_pairs, df


def load_marker_pairs_from_ranked_markers(
    csv_path,
    clusters_to_focus,
    max_pairs_per_cluster=3,
    annotation_lookup=None,
    cluster_col="group",
    gene_col="names",
):
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Ranked-marker CSV was not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required = {cluster_col, gene_col}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s) in {csv_path}: {sorted(missing)}")

    annotation_lookup = annotation_lookup or {}
    focus = [normalize_cluster_id(cluster) for cluster in clusters_to_focus]
    work = df.copy()
    work["__cluster__"] = work[cluster_col].map(normalize_cluster_id)
    marker_pairs = {}

    for cluster_id in focus:
        genes = (
            work.loc[work["__cluster__"] == cluster_id, gene_col]
            .dropna()
            .astype(str)
            .str.strip()
        )
        genes = [gene for gene in dict.fromkeys(genes) if gene]
        pairs = []
        for idx in range(0, min(len(genes), max_pairs_per_cluster * 2), 2):
            markers = genes[idx : idx + 2]
            if not markers:
                continue
            pairs.append({"markers": markers, "source_pair": " + ".join(markers)})
            if len(pairs) >= max_pairs_per_cluster:
                break
        marker_pairs[cluster_id] = {
            "annotation": annotation_lookup.get(cluster_id, "Ranked marker cluster"),
            "pairs": pairs,
        }

    print(f"Loaded fallback ranked-marker CSV: {csv_path}")
    print("Marker pairs were created from adjacent top-ranked marker genes because no usable marker-support Excel file was found.")
    return marker_pairs, df


def summarize_marker_availability(adata, marker_pairs):
    summary = {}
    for cluster_id, info in marker_pairs.items():
        available_genes = []
        missing_genes = []
        available_pairs = []
        for pair in info["pairs"]:
            present = []
            missing = []
            for gene in pair["markers"]:
                matched = match_gene_to_var_names(adata, gene)
                if matched is None:
                    missing.append(gene)
                    missing_genes.append(gene)
                else:
                    present.append(matched)
                    available_genes.append(matched)
            if present:
                available_pairs.append({"source_pair": pair["source_pair"], "available": present, "missing": missing})
            else:
                print(f"Warning: cluster {cluster_id} marker pair skipped; both markers missing: {pair['source_pair']}")

        summary[cluster_id] = {
            "annotation": info["annotation"],
            "available_genes": list(dict.fromkeys(available_genes)),
            "missing_genes": list(dict.fromkeys(missing_genes)),
            "available_pairs": available_pairs,
        }

    for cluster_id, info in summary.items():
        print(f"\nCluster {cluster_id} | {info['annotation']}")
        print("  Plot genes:", ", ".join(info["available_genes"]) if info["available_genes"] else "none")
        if info["missing_genes"]:
            print("  Missing genes:", ", ".join(info["missing_genes"]))
    return summary


def plot_cluster_overview_umap(adata, cluster_col, clusters_to_focus, basis=None, point_size=1.0):
    if cluster_col not in adata.obs:
        raise KeyError(f"Cluster column '{cluster_col}' is not present in adata.obs")
    coords, coord_key = get_umap_coordinates(adata, basis=basis)
    clusters = adata.obs[cluster_col].map(normalize_cluster_id).astype("category")
    categories = list(clusters.cat.categories)

    cmap = plt.get_cmap("tab20")
    color_map = {cat: to_hex(cmap(idx % 20)) for idx, cat in enumerate(categories)}
    if f"{cluster_col}_colors" in adata.uns:
        stored = list(adata.uns[f"{cluster_col}_colors"])
        color_map.update({cat: stored[idx] for idx, cat in enumerate(categories[: len(stored)])})

    fig, ax = plt.subplots(figsize=(5.2, 3.2), dpi=300)
    for cluster_id in categories:
        mask = clusters == cluster_id
        alpha = 0.95 if cluster_id in set(map(str, clusters_to_focus)) else 0.18
        size = point_size * (2.4 if cluster_id in set(map(str, clusters_to_focus)) else 1.0)
        ax.scatter(coords[mask, 0], coords[mask, 1], s=size, c=color_map[cluster_id], alpha=alpha, linewidths=0)

    for cluster_id in map(str, clusters_to_focus):
        mask = clusters == cluster_id
        if mask.sum() == 0:
            print(f"Warning: focused cluster {cluster_id} is not present in '{cluster_col}'")
            continue
        x, y = np.median(coords[mask, 0]), np.median(coords[mask, 1])
        ax.text(x, y, cluster_id, ha="center", va="center", fontsize=7, weight="bold",
                bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.8})

    ax.set_title(f"Cluster overview | {cluster_col}\nCoordinates: {coord_key}", fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_frame_on(False)
    plt.tight_layout()
    plt.show()
    return fig, ax


def plot_umap_gene_expression_panels(
    adata,
    genes,
    cluster_id=None,
    annotation=None,
    layer=None,
    basis=None,
    n_bins=20,
    genes_per_row=3,
    size_bg=1,
    size_expr=2,
):
    genes = list(dict.fromkeys([gene for gene in genes if gene]))
    if not genes:
        print(f"No available genes to plot for cluster {cluster_id}.")
        return None

    coords, _ = get_umap_coordinates(adata, basis=basis)
    ncols = min(int(genes_per_row), len(genes))
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.35 * ncols, 2.45 * nrows), dpi=300, squeeze=False)

    suptitle_bits = []
    if cluster_id is not None:
        suptitle_bits.append(f"Cluster {cluster_id}")
    if annotation:
        suptitle_bits.append(str(annotation))
    if suptitle_bits:
        fig.suptitle(" | ".join(suptitle_bits), fontsize=9, y=0.995)

    for idx, gene in enumerate(genes):
        ax = axes[idx // ncols, idx % ncols]
        expr, matched_gene = get_expression_vector(adata, gene, layer=layer)
        positive = expr > 0

        ax.scatter(coords[:, 0], coords[:, 1], s=size_bg, c="lightgrey", alpha=0.55, linewidths=0)
        if positive.sum() == 0:
            print(f"Skipping expression overlay for {matched_gene}: no positive values")
        else:
            values = expr.loc[positive]
            n_actual_bins = min(int(n_bins), values.nunique())
            if n_actual_bins <= 1:
                color_values = np.full(values.shape[0], 0.7)
            else:
                cats = pd.qcut(values, q=n_actual_bins, duplicates="drop", labels=False)
                denom = max(int(np.nanmax(cats)), 1)
                color_values = 0.35 + 0.60 * (np.asarray(cats, dtype=float) / denom)
            ax.scatter(
                coords[positive.values, 0],
                coords[positive.values, 1],
                s=size_expr,
                c=[to_hex(plt.cm.Reds(value)) for value in color_values],
                alpha=0.95,
                linewidths=0,
            )

        ax.set_title(matched_gene, fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_frame_on(False)

    for idx in range(len(genes), nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.94] if suptitle_bits else None)
    plt.show()
    return fig, axes