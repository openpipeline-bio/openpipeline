import mudata as mu

##VIASH START
par = {
  "input": "output_test/split_modalities/foo_types.csv",
  "mod_dir": "output_test/split_modalities/h5mu"
}

meta = {
    "resources_dir": "resources_test/pbmc_1k_protein_v3"
}

##VIASH END

print ("Loading data", flush=True)
data = mu.read_h5mu(par["input"])

assert "X_umap" in data.mod["rna"].obsm, "X_umap not found in .obsm"
assert data.mod["rna"].obsm["X_umap"].shape[1] == 2, f"X_umap has wrong shape expected 2 n_comp but got {data.mod['rna'].obsm['X_umap'].shape[1]}"


print("Test successful!", flush=True)