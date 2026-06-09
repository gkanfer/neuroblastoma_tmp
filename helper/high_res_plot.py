import tifffile
import os
from spatialdata.models import Image2DModel
import spatialdata_io
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
# from spatialdata.models import ShapesModels
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
from helper.HE_process import SpatioHEprocess

#andata_file, andata_path, sdata_path,sdata_file,data_set,cluster_cols, HE_image_path,image_only,HE_crop_path

class HighPlotStatic(SpatioHEprocess):
    def __init__(self,cluster_colors_dict,min_coordinate,max_coordinate,image_element = 'DAPI',dpi = 100,axis_off=False,legend_loc = None, *args,**kwargs):
        self.cluster_colors_dict = cluster_colors_dict
        self.min_coordinate = min_coordinate
        self.max_coordinate= max_coordinate
        self.image_element = image_element
        self.dpi = dpi
        self.axis_off = axis_off
        self.legend_loc = legend_loc # 'right margin'/"bottom margin"
        super().__init__(HE_image_path= None,image_only = None,HE_crop_path= None,morphology_focus= True,*args,**kwargs)
        
    
    def crop(self):
        sdata_crop = bounding_box_query(
            self.sdata,
            min_coordinate=self.min_coordinate,
            max_coordinate=self.max_coordinate,
            axes=("y", "x"),
            target_coordinate_system="global",)
        return sdata_crop
    
        
    def plotstatic(self):
        if self.min_coordinate is not None and self.max_coordinate is not None:
            sdata_curr = self.crop()
        else:
            sdata_curr = self.sdata
        cluster_palette = list(self.cluster_colors_dict.values())
        plt.style.use("dark_background")
        ax = sdata_curr.pl.render_images(
            element="morphology_focus", 
            channel='DAPI',cmap= 'grey', colorbar=False).pl.render_labels(element='cell_labels', color=self.cluster_cols[0],groups=list(self.cluster_colors_dict.keys()), palette=cluster_palette, na_color='default',contour_px = 3, outline_alpha=0.0, fill_alpha=1, scale=None, table_name='table', table_layer=None).pl.show(coordinate_systems="global",title=f" ",figsize=(10, 10),dpi = self.dpi,na_in_legend=False,legend_loc = self.legend_loc,return_ax=True)
        if self.axis_off:
            ax.set_axis_off()
        plt.show()