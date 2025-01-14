import mudata as mu
import anndata as ad

## VIASH START
par = {
  "input_dataset": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  "input_output": "work/90/56062df29c88150755a174b63bdb82/_viash_par/input_output_1/run.pymde.output.h5ad",
  "input_modality": "rna",
  "output_obsm_key": "X_dr",
  "output": "output.h5mu"
}
## VIASH END

print("Reading h5mu file", flush=True)
mdata = mu.read_h5mu(par["input_dataset"])
adata = ad.read_h5ad(par["input_output"])

adata_dest = mdata.mod[par["input_modality"]]
adata_dest.obsm[par["output_obsm_key"]] = adata.obsm["X_emb"]

print("Writing h5ad file", flush=True)
mdata.write_h5mu(par["output"])
