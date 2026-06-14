# mamba activate spatiocore_spatial_env
import os
import re
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import random

pathout = '/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig'
andata = sc.read_h5ad(os.path.join(pathout, "adata_nb_processed.h5ad"))

sc.tl.umap(andata, neighbors_key='30neig')

pathout = '/data/HiTIF/data/spatialomics/neuroblastoma/data/andata_orig'
andata.write_h5ad(os.path.join(pathout, "adata_nb_processed.h5ad"))