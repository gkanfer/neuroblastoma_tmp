# mamba activate spatiocore_spatial_env
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
# from napari_spatialdata import Interactive
from spatialdata.models import TableModel
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from spatialdata import SpatialData
from napari.utils.colormaps import Colormap
from napari_spatialdata import Interactive
import spatialdata as sd
from scipy import ndimage as ndi
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import pickle

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
import sys


sys.path.insert(
    0,
    "/data/kanferg/Sptial_Omics/projects/NguyenLab/neuroblastoma"
)

from helper.utils import add_cluster_from_pickle
from helper.high_res_plot import PunchPlotter



pathad = '/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig'
andata = sc.read_h5ad(os.path.join(pathad, "adata_nb_processed_v1.h5ad"))
