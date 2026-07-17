import argparse
import os
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd


project_dir = Path("/data/kanferg/Sptial_Omics/projects/NguyenLab/neuroblastoma")
sys.path.insert(0, str(project_dir))

from helper.high_res_plot import PunchPlotter


pathad = "/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig"
path_zarr = "/data/HiTIF/data/spatialomics/neuroblastoma/data/SCAF4480/PrimaryAnalysisOutput/SCAF4480_PA_xenium/nb_orig_label_v1.zarr"

annotation_colors = {
    "DISCARD": "#F2F2F2",

    "Adrenergic_Proliferating_NB": "#FF3B30",
    "Transitional_NB": "#00D1FF",
    "Adrenergic_StemLike_NB": "#39FF14",
    "Adrenergic_Cycling_NB": "#FF4FD8",
    "Adrenergic_MYCNHigh_Metabolic_NB": "#FFD60A",
    "Mesenchymal_NB": "#8E5CFF",

    "Schwann_CAF": "#C49A6C",
    "M2_TAM": "#4292C6",
    "CD8_Cytotoxic_TCell_PartiallyExhausted": "#D62728",
    "MonocyteDerived_TissueResident_Macrophage": "#08519C",
    "Endothelial": "#17BECF",
    "Adrenocortical_Steroidogenic": "#F16913",
}

nb_annotations = [
    "Adrenergic_Proliferating_NB",
    "Transitional_NB",
    "Adrenergic_StemLike_NB",
    "Adrenergic_Cycling_NB",
    "Adrenergic_MYCNHigh_Metabolic_NB",
    "Mesenchymal_NB",
]

non_nb_annotations = [
    "Schwann_CAF",
    "M2_TAM",
    "CD8_Cytotoxic_TCell_PartiallyExhausted",
    "MonocyteDerived_TissueResident_Macrophage",
    "Endothelial",
    "Adrenocortical_Steroidogenic",
]

nb_palette = {label: annotation_colors[label] for label in nb_annotations}
non_nb_palette = {label: annotation_colors[label] for label in non_nb_annotations}


def make_punch_table(sdata):
    sample_table = pd.read_csv(project_dir / "data_in/05082026_Sample_punch_id.csv")
    sample_table["Condition_clean"] = sample_table["Condition"].replace({
        "IL-15.CAR": "IL15",
        "IL-21.CAR": "IL21",
        "IL-15.IL-21.CAR": "Both",
    })

    punch_table = sample_table.dropna(subset=["punch_id"]).copy()
    punch_table["punch_id"] = punch_table["punch_id"].astype(str).str.split(",")
    punch_table = punch_table.explode("punch_id")
    punch_table["punch_id"] = punch_table["punch_id"].str.strip()
    punch_table = punch_table[punch_table["punch_id"].isin(list(sdata.shapes.keys()))].copy()

    return punch_table


def make_punch_titles(punch_table):
    punch_titles = {}
    for _, row in punch_table.iterrows():
        title_bits = [f"Punch {row['punch_id']}"]
        if pd.notna(row.get("Type")):
            title_bits.append(str(row["Type"]))
        if pd.notna(row.get("Ear tag")):
            title_bits.append(f"Mouse {row['Ear tag']}")
        punch_titles[row["punch_id"]] = " | ".join(title_bits)
    return punch_titles


def make_highres_panel(
    sdata,
    punches,
    active_labels,
    palette,
    panel_title,
    save_file,
    map_annotation,
    punch_titles,
    dpi=300,
    panel_dpi=300,
):
    n_punches = len(punches)
    ncols = min(4, n_punches)
    nrows = int(np.ceil(n_punches / ncols))

    with tempfile.TemporaryDirectory() as tmpdir:
        image_paths = []
        for punch_id in punches:
            print(f"  rendering {panel_title}: punch {punch_id}", flush=True)
            image_path = Path(tmpdir) / f"{panel_title}_{punch_id}.png"
            sdata_obj = PunchPlotter(
                sdata=sdata,
                element_name=str(punch_id),
                legend_loc=None,
                dpi=dpi,
            )
            ax = sdata_obj.plot_punch(
                colname="annotation_v1",
                colname_map=map_annotation,
                palette=palette,
                active_labels=active_labels,
                figsize=(6, 6),
                title=" ",
                legend_loc=None,
                show=False,
                save_path=image_path,
                dpi=dpi,
            )
            plt.close(ax.figure)
            image_paths.append(image_path)

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(4.0 * ncols, 4.4 * nrows),
            dpi=panel_dpi,
            squeeze=False,
        )
        axes_flat = axes.ravel()

        for ax, punch_id, image_path in zip(axes_flat, punches, image_paths):
            ax.imshow(mpimg.imread(image_path))
            ax.set_title(punch_titles.get(str(punch_id), f"Punch {punch_id}"), fontsize=9, pad=4)
            ax.set_axis_off()

        for ax in axes_flat[len(punches):]:
            ax.set_axis_off()

        handles = [Patch(facecolor=palette[label], edgecolor="none", label=label) for label in active_labels]
        legend_ncol = 2 if len(active_labels) <= 6 else 3
        fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=legend_ncol,
            frameon=False,
            fontsize=8,
        )
        fig.suptitle(panel_title, fontsize=12, y=0.995)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.91, bottom=0.22, wspace=0.04, hspace=0.12)
        fig.savefig(save_file, dpi=panel_dpi, bbox_inches="tight")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--conditions", nargs="+", default=["CAR", "IL15", "IL21", "Both"])
    parser.add_argument("--max-punches", type=int, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--output-dir", type=Path, default=project_dir / "figures")
    args = parser.parse_args()

    figures_dir = args.output_dir
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("loading AnnData", flush=True)
    andata = sc.read_h5ad(os.path.join(pathad, "adata_nb_processed_v1.h5ad"))
    map_annotation = dict(zip(andata.obs["cell_id"], andata.obs["annotation_v1"]))

    print("loading SpatialData", flush=True)
    sdata = sd.read_zarr(path_zarr)

    punch_table = make_punch_table(sdata)
    punch_titles = make_punch_titles(punch_table)

    for condition in args.conditions:
        punches = punch_table.loc[punch_table["Condition_clean"] == condition, "punch_id"].tolist()
        if args.max_punches is not None:
            punches = punches[:args.max_punches]
        if not punches:
            print(f"skipping {condition}: no punches found", flush=True)
            continue

        print(f"{condition}: {len(punches)} punches", flush=True)
        make_highres_panel(
            sdata=sdata,
            punches=punches,
            active_labels=nb_annotations,
            palette=nb_palette,
            panel_title=f"{condition} punches | NB states",
            save_file=figures_dir / f"highres_{condition}_NB.png",
            map_annotation=map_annotation,
            punch_titles=punch_titles,
            dpi=args.dpi,
            panel_dpi=args.dpi,
        )
        make_highres_panel(
            sdata=sdata,
            punches=punches,
            active_labels=non_nb_annotations,
            palette=non_nb_palette,
            panel_title=f"{condition} punches | Non-NB compartments",
            save_file=figures_dir / f"highres_{condition}_nonNB.png",
            map_annotation=map_annotation,
            punch_titles=punch_titles,
            dpi=args.dpi,
            panel_dpi=args.dpi,
        )


if __name__ == "__main__":
    main()
