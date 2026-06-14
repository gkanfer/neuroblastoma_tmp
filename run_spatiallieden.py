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

pathout = '/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig'
pathdict = '/data/HiTIF/data/spatialomics/neuroblastoma/data/misc'

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

# adata_sub.obs['cluster_60neig_resolution_0.6'] = adata_sub.obs['cluster_60neig_resolution_0.6'].astype('str')
sq.gr.spatial_neighbors(
    andata,
    spatial_key="spatial",
    library_key="punch_id",
    coord_type="generic",
    delaunay=False,
    radius=30,
)
andata.obsp["spatial_connectivities"] = sl.distance2connectivity(
    andata.obsp["spatial_distances"]
)


sl.spatialleiden(andata,layer_ratio=1.0, directed=(False, False), seed=seed, key_added="spatialleiden_resolution_1")

sl.spatialleiden(andata, layer_ratio=1.3, directed=(False, False), seed=seed, key_added="spatialleiden_resolution_1_3")

sl.spatialleiden(andata, layer_ratio=1.5, directed=(False, False), seed=seed, key_added="spatialleiden_resolution_1_5")

sl.spatialleiden(andata, layer_ratio=1.8, directed=(False, False), seed=seed, key_added="spatialleiden_resolution_1_8")


import pickle
dict_1 = dict(zip(andata.obs["cell_id"], andata.obs['spatialleiden_resolution_1'].astype(str)))
dict_1_3 = dict(zip(andata.obs["cell_id"], andata.obs['spatialleiden_resolution_1_3'].astype(str)))
dict_1_5 = dict(zip(andata.obs["cell_id"], andata.obs['spatialleiden_resolution_1_5'].astype(str)))
dict_1_8 = dict(zip(andata.obs["cell_id"], andata.obs['spatialleiden_resolution_1_8'].astype(str)))

with open(os.path.join(pathdict, "stlieden_1.pkl"), "wb") as f:
    pickle.dump(dict_1, f)

with open(os.path.join(pathdict, "stlieden_1_3.pkl"), "wb") as f:
        pickle.dump(dict_1_3, f)

with open(os.path.join(pathdict, "stlieden_1_5.pkl"), "wb") as f:
    pickle.dump(dict_1_5, f)
    
with open(os.path.join(pathdict, "stlieden_1_8.pkl"), "wb") as f:
    pickle.dump(dict_1_8, f)    