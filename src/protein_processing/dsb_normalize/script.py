from warnings import warn
import os
import scanpy as sc
import muon as mu
import numpy as np
from mudata import MuData
import pandas as pd
from anndata import AnnData
from muon import prot as pt

## VIASH START
par = {
    "data_raw": "../dsb_index/mudata_raw.h5mu",
    "output": "dsb_processed",
    "cell_index":  "../dsb_index/dsb_processed/cell_idx.csv",
    "empty_index": "../dsb_index/dsb_processed/empty_idx.csv",
    "pseudocount": 10,
    "denoise_counts": True,
    "isotype_controls": None,
    "add_layer": False,
    "random_state": None
}
## VIASH END

if par['data_raw'] is not None:
    if par['data_raw'].endswith('h5ad'):
        raw_data = sc.read_h5ad(par['data_raw'])
    elif par['data_raw'].endswith('h5mu'):
        raw_data = mu.read_h5mu(par['data_raw'])
    elif par['data_raw'].endswith('h5'):
        raw_data = mu.read_10x_h5(par['data_raw'])
    else:
        raise TypeError("data_raw must be an AnnData or a MuData object with 'prot' modality")
    if "prot" not in raw_data.mod:
            raise TypeError("data_raw must be an AnnData or a MuData object with 'prot' modality")
else:
    raise ValueError( "Raw data is not available.")

if par['cell_index'] is None and par['empty_index'] is None:
    raise ValueError( "Given the unfiltered object data_raw, at least one index file must be "
                      "provided for foreground and background signals.")
                      
cells = None
empty = None

if par['empty_index'] is not None:
    empty_idx = pd.read_csv(par["empty_index"], header=None).iloc[:, 0].tolist()
    empty = raw_data[raw_data.obs_names.isin(empty_idx)]
if par['cell_index'] is not None:
    cell_idx = pd.read_csv(par["cell_index"], header=None).iloc[:, 0].tolist()
    cells = raw_data[raw_data.obs_names.isin(cell_idx)]

if empty is None:
    empty = raw_data[~raw_data.obs_names.isin(cell_idx)]
if cells is None:
    cells = raw_data[~raw_data.obs_names.isin(empty_idx)]

pt.pp.dsb(cells, empty, isotype_controls=par['isotype_controls'], pseudocount=par["pseudocount"],
    denoise_counts=par["denoise_counts"], add_layer=par["add_layer"],
    random_state=par["random_state"])

if not os.path.exists(par["output"]):
    os.makedirs(par["output"])
if isinstance(cells, MuData):
    cells.write(os.path.join(par["output"], "normalized.h5mu"))
else:
    cells.write(os.path.join(par["output"], "normalized.h5ad"))
    
