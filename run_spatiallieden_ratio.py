# mamba activate spatialleiden-cupy 
import scanpy as sc
import spatialleiden as sl
import squidpy as sq
import numpy as np
import os
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import random
import pandas as pd
import argparse
import pickle
import h5py
import os
import pickle

pathout = '/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig'
pathdict = '/data/HiTIF/data/spatialomics/neuroblastoma/data/misc/res_07_ratio1p5'

pathout = "/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig"
fn = os.path.join(pathout, "adata_nb_processed.h5ad")

with h5py.File(fn, "a") as f:
    for key in [
        "uns/log1p/base",
        "uns/pca/params/mask_var",
    ]:
        if key in f:
            del f[key]

andata = sc.read_h5ad(fn)




seed = 1948
andata.obsp['connectivities'] = andata.obsp['30neig_connectivities']
andata.obsp['distances'] = andata.obsp['30neig_distances']


radius = np.linspace(3,50,10,dtype = np.float64)
for r in radius:
    sq.gr.spatial_neighbors(
        andata,
        spatial_key="spatial",
        library_key="punch_id",
        coord_type="generic",
        delaunay=False,
        radius=r,
    )
    andata.obsp["spatial_connectivities"] = sl.distance2connectivity(andata.obsp["spatial_distances"])
    r_int = int(r)
    sl.spatialleiden(andata,resolution = (0.7,0.7), layer_ratio=1.5, directed=(False, False), seed=seed, key_added=f"spatialleiden_resolution_1p5_r{r_int}")
    dict_curr = dict(zip(andata.obs["cell_id"], andata.obs[f'spatialleiden_resolution_1p5_r{r_int}'].astype(str)))
    with open(os.path.join(pathdict, f"stlieden_ratio1p5_r{r_int}.pkl"), "wb") as f:
        pickle.dump(dict_curr, f)  