from warnings import warn
import scanpy as sc
import muon as mu
import numpy as np
from mudata import MuData
import pandas as pd
import os
from anndata import AnnData

## VIASH START
par = {
    "data_raw": "mudata_raw.h5mu",
    "output": "dsb_processed",
    "cell_index": None,
    "empty_counts_range": [1.5, 2.8],
    "cell_counts_range": None
}

## VIASH END

if par['data_raw'].endswith('h5mu'):
    raw_data = mu.read_h5mu(par['data_raw'])
elif par['data_raw'].endswith('h5'):
    raw_data = mu.read_10x_h5(par['data_raw'])
else:
    raise TypeError("data_raw must be a MuData object with 'prot' and 'rna' modalities")

if "prot" not in raw_data.mod or "rna" not in raw_data.mod:
    raise TypeError("Raw data does not contain 'prot' or 'rna' modalities")
if raw_data.mod["rna"].n_obs != raw_data.mod["prot"].n_obs:
    raise ValueError("different numbers of cells in 'rna' and 'prot' modalities.")

droplet_barcode = raw_data.mod["prot"].obs_names
if par["cell_index"] is not None:
    cell_barcode = pd.read_csv(par["cell_index"],header=None).iloc[:, 0].tolist()
    empty_barcode = list(set(droplet_barcode).difference(cell_barcode))
else:
    cell_barcode = None
    empty_barcode = None

log10umi = np.log10(np.asarray(raw_data.mod["rna"].X.sum(axis=1)).squeeze() + 1)

if par['empty_counts_range'] is not None:
    if len(par['empty_counts_range']) != 2:
        raise ValueError("Invalid count ranges provided for the empty droplets.")
    if par['cell_counts_range'] is not None and max(*par['empty_counts_range']) > min(*par['cell_counts_range']):
        raise ValueError("Overlapping count ranges")
    empty_idx = np.where(
        (log10umi >= min(*par['empty_counts_range'])) & (log10umi < max(*par['empty_counts_range'])))[0]
    empty_idx = droplet_barcode[empty_idx]
    if empty_barcode is not None:
        empty_barcode = list(set(empty_barcode) & set(empty_idx))
    else:
        empty_barcode = empty_idx


if par['cell_counts_range'] is not None:
    if len(par['cell_counts_range']) != 2:
        raise ValueError("Invalid count ranges provided for true cells.")
    cell_idx = np.where(
    (log10umi >= min(*par['cell_counts_range'])) & (log10umi < max(*par['cell_counts_range'])))[0]
    cell_idx = droplet_barcode[cell_idx]
    if cell_barcode is not None:
        cell_barcode = list(set(cell_barcode) & set(cell_idx))
    else:
        cell_barcode = cell_idx

if empty_barcode is None:
    if cell_barcode is None:
        raise ValueError("Neither cell_index nor counts ranges for empty droplets "
                         "or cells provided for filtering empty droplets.")
    empty_barcode = list(set(droplet_barcode).difference(cell_barcode))
elif cell_barcode is None:
    cell_barcode = list(set(droplet_barcode).difference(empty_barcode))

if not os.path.exists(par["output"]):
    os.makedirs(par["output"])
pd.DataFrame(cell_barcode).to_csv(os.path.join(par["output"], "cell_idx.csv"), header=None, index=None)
pd.DataFrame(empty_barcode).to_csv(os.path.join(par["output"], "empty_idx.csv"), header=None, index=None)
