import tifffile
import os
from spatialdata.models import Image2DModel
import spatialdata_io
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
from spatialdata.models import ShapesModel
import anndata
import scanpy as sc
from spatialdata import SpatialData
from spatialdata import polygon_query
from spatialdata import bounding_box_query
from napari_spatialdata import Interactive
from spatialdata.models import TableModel
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from spatialdata import SpatialData
import spatialdata as sd
from scipy import ndimage as ndi
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from spatialdata.transformations import (
    BaseTransformation,
    Identity,
    Sequence,
    align_elements_using_landmarks,
    get_transformation,
    set_transformation,
    Affine
)
import pickle


class SpatioHEprocess:
    def __init__(self, andata_file, andata_path, sdata_path,sdata_file,data_set,cluster_cols, HE_image_path,image_only = False,HE_crop_path = None,morphology_focus= False,*args,
    **kwargs):
        self.andata_file = andata_file
        self.andata_path = andata_path
        self.sdata_path = sdata_path
        self.sdata_file = sdata_file
        self.data_set = data_set
        self.cluster_cols = cluster_cols #could be list e.g. ['spatialleiden_resolution_1_8,spatialleiden_resolution_1_5']
        self.HE_image_path = HE_image_path
        self.image_only = image_only
        self.HE_crop_path = HE_crop_path # '/data/HiTIF/data/spatialomics/cystic_duct/images_HE/slide_1_r3.zarr'
        self.morphology_focus = morphology_focus
        self.sdata, self.andata = self.data_prepare()
        
        
    def read_files(self):
        sdata = spatialdata_io.xenium(
                    path=os.path.join(self.sdata_path, self.sdata_file),
                    n_jobs=8,
                    cells_boundaries=True,
                    nucleus_boundaries=False,
                    morphology_focus=self.morphology_focus,
                    cells_as_circles=False,)
        andata = sc.read_h5ad(os.path.join(self.andata_path, self.andata_file))
        return sdata, andata
    
    @staticmethod
    def add_cluster_from_pickle(andata, pkl_file, new_col, match_col="orig_obs_name_sample"):
        """
        Add clustering results from a pickle dictionary to andata.obs.

        The pickle file should contain:
            key   = original obs name
            value = cluster label
        """

        with open(pkl_file, "rb") as f:
            cluster_dict = pickle.load(f)

        # make sure dictionary keys are strings
        cluster_dict = {str(k): str(v) for k, v in cluster_dict.items()}

        # map cluster labels into AnnData
        andata.obs[new_col] = andata.obs[match_col].map(cluster_dict)

        # report mapping quality
        n_total = andata.n_obs
        n_missing = andata.obs[new_col].isna().sum()

        print(f"{new_col}")
        print(f"  total cells:   {n_total}")
        print(f"  matched cells: {n_total - n_missing}")
        print(f"  missing cells: {n_missing}")

        # convert to category for plotting
        andata.obs[new_col] = andata.obs[new_col].astype("category")

        return andata
    
    
    def data_prepare(self):
        sdata,andata = self.read_files()
        andata = andata[andata.obs['sample'] == self.data_set].copy()
        andata.obs["orig_obs_name"] = andata.obs["orig_obs_name"].astype(str)
        andata.obs["sample"] = andata.obs["sample"].astype(str)
        andata.obs['orig_obs_name_sample'] = andata.obs['orig_obs_name']+'_'+andata.obs['sample']
        from pathlib import Path
        pkl_dir = Path(".")
        andata = self.add_cluster_from_pickle(andata,
            pkl_dir / "stlieden_1_8.pkl",
            "spatialleiden_resolution_1_8"
        )

        cluster_lookup = andata.obs[["orig_obs_name"] + self.cluster_cols].copy()

        # make sure each orig_obs_name appears only once
        cluster_lookup = cluster_lookup.drop_duplicates("orig_obs_name")

        cluster_lookup = cluster_lookup.set_index("orig_obs_name")

        # -----------------------------
        # match to SpatialData table
        # -----------------------------
        table = sdata["table"]

        # make sure sdata table index is string
        table.obs.index = table.obs.index.astype(str)

        # -----------------------------
        # add cluster columnsxw
        # -----------------------------
        for col in self.cluster_cols:
            table.obs[col] = table.obs.index.map(cluster_lookup[col])
            table.obs[col] = table.obs[col].astype("category")

            n_missing = table.obs[col].isna().sum()
            n_total = table.n_obs

            print(f"{col}")
            print(f"  total cells in sdata table: {n_total}")
            print(f"  matched cells:             {n_total - n_missing}")
            print(f"  missing cells:             {n_missing}")
        return sdata, andata

    def read_HandE_image(self):
        if self.HE_crop_path:
            sdata_small = sd.read_zarr(self.HE_crop_path)
            return sdata_small
        img_array = tifffile.imread(self.HE_image_path)
    
        # convert from yxc -> cyx
        img_array = np.transpose(img_array, (2, 0, 1))
        
        print("Original shape:", img_array.shape)
        print("Original dtype:", img_array.dtype)
        # img_array = img_array[:,:,:10_000]
        # create SpatialData object
        image_element_small = Image2DModel.parse(img_array, scale_factors=(2, 2, 2),rgb=True)
        # rename_coordinate
        if self.sdata is None:
            print("sdata is None")
            return
        else:
            sdata_temp = self.sdata
            #sdata_temp.rename_coordinate_systems({"global": "labels_coor"})
            if self.image_only:
                sdata_small = sd.SpatialData(images={'HE': image_element_small})
            else:
                sdata_small = sd.SpatialData(images={'HE': image_element_small},shapes={'cell_labels': ShapesModel.parse(sdata_temp['cell_boundaries'])})
            return sdata_small
    
    def build_cells_for_alignment(self):
        self.sdata.shapes['cell_boundaries'][self.cluster_cols[0]] = self.sdata["table"].obs[self.cluster_cols[0]]
        self.sdata.shapes['cell_boundaries_orig'] = self.sdata.shapes['cell_boundaries'].copy()
        self.sdata.shapes['cell_boundaries']['geometry'] = self.sdata.shapes['cell_boundaries'].geometry.centroid
        self.sdata.shapes['cell_boundaries']['radius'] = 1
        sdata_label = sd.SpatialData(shapes={'cell_labels': ShapesModel.parse(self.sdata['cell_boundaries'])})
        sdata_label.rename_coordinate_systems({"global": "label_coor"})
        sdata_label_orig = sd.SpatialData(shapes={'cell_labels': ShapesModel.parse(self.sdata['cell_boundaries_orig'])})
        sdata_label_orig.rename_coordinate_systems({"global": "label_coor"})
        return sdata_label,sdata_label_orig
    
    def build_id2cluster_mapping(self):
        id2cluster = dict(zip(self.sdata["table"].obs['cell_id'], self.sdata["table"].obs[self.cluster_cols[0]]))
        return id2cluster
            
            
            
            
def postpone_transformation(sdata: SpatialData,transformation: BaseTransformation,source_coordinate_system: str,target_coordinate_system: str):
        for element_type, element_name, element in sdata._gen_elements():
            old_transformations = get_transformation(element, get_all=True)
            if source_coordinate_system in old_transformations:
                old_transformation = old_transformations[source_coordinate_system]
                sequence = Sequence([old_transformation, transformation])
                set_transformation(element, sequence, target_coordinate_system)
                
                
def run_affine(sdata_align,HE_name = "HE",new_coordinate_system = "aligned"):
    affine = align_elements_using_landmarks(
    references_coords=sdata_align["lanmark_label"],
    moving_coords=sdata_align["landmark_HE"],
    reference_element=sdata_align["cell_labels"],
    moving_element=sdata_align[HE_name],
    reference_coordinate_system="label_coor",
    moving_coordinate_system="HE_coor",
    new_coordinate_system="aligned")
    return affine            