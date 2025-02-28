import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

##VIASH START
par = {"input": "output.h5mu", "orig_input": "input.5mu"}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"

# Check atac modality
assert_annotation_objects_equal(
    input.mod["atac"], output.mod["atac"], promote_precision=True
)

# Check rna modality
assert "X_umap" in output.mod["rna"].obsm, "X_umap not found in .obsm"
assert (
    output.mod["rna"].obsm["X_umap"].shape[1] == 2
), f"X_umap has wrong shape expected 2 n_comp but got {output.mod['rna'].obsm['X_umap'].shape[1]}"
assert "pca_variance" in output.mod["rna"].uns
assert "pca_loadings" in output.mod["rna"].varm


print("Test successful!", flush=True)
