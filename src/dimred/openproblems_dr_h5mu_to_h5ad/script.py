import mudata as mu
import anndata as ad

## VIASH START
par = {
  "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  "input_modality": "rna",
  "input_layer_counts": "log_normalized",
  "input_layer_normalized": "log_normalized",
  "input_var_hvg_score": None,
  "output": "output.h5mu"
}
## VIASH END

print("Reading h5mu file", flush=True)
mdata = mu.read_h5mu(par["input"])

print("Transforming to anndata", flush=True)
def get_matrix(mdata, modality, layer):
  if layer is None:
    return mdata.mod[modality].X
  return mdata.mod[modality].layers[layer]

# create var
var = mdata.mod[par["input_modality"]].var[[]]

if par["input_var_hvg_score"] is not None:
  assert par["input_var_hvg_score"] in var.columns, f"Variable {par['input_var_hvg_score']} not found in var"
  var["hvg_score"] = mdata.mod[par["input_modality"]].var[par["input_var_hvg_score"]]

# create anndata
adata = ad.AnnData(
  layers={
    "counts": get_matrix(mdata, par["input_modality"], par["input_layer_counts"]),
    "normalized": get_matrix(mdata, par["input_modality"], par["input_layer_normalized"])
  },
  obs=mdata.mod[par["input_modality"]].obs[[]],
  var=var,
  uns={
    "dataset_id": "dummy",
    "normalization_id": "dummy"
  }
)

print("Writing h5ad file", flush=True)
adata.write_h5ad(par["output"], compression="gzip")
