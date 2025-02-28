import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal
from openpipelinetest_utils.utils import remove_annotation_column

##VIASH START
par = {
    "input": "resources_test/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu",
    "orig_input": "test.h5mu",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"

# Check vdj_t modality
# Allow X_umap to be overwritten
input_vdj = input.mod["vdj_t"]
# del input_vdj.obsm['X_umap']
output_vdj = output.mod["vdj_t"]
# del output_vdj.obsm['X_umap']
assert_annotation_objects_equal(input_vdj, output_vdj, promote_precision=True)

# Check prot modality
# Ignore the PCA layer and its derivatives, as its allowed to be overwritten for this test.
input_prot = input.mod["prot"]
del input_prot.varm["pca_loadings"]
del input_prot.obsm["X_pca"]
del input_prot.obsm["X_umap"]
output_prot = output.mod["prot"]
del output_prot.varm["pca_loadings"]
del output_prot.obsm["X_pca"]
del output_prot.obsm["X_umap"]
assert_annotation_objects_equal(input_prot, output_prot, promote_precision=True)


# Check rna modality
# Allow the highly variable genes and PCA + derivatives to be overwritten
input_rna = input.mod["rna"]
input_rna = remove_annotation_column(input_rna, "filter_with_hvg", "var")
del input_rna.varm["pca_loadings"]
del input_rna.obsm["X_pca"]
del input_rna.obsm["X_umap"]
del input_rna.layers["log_normalized"]
output_rna = output.mod["rna"]
output_rna = remove_annotation_column(output_rna, "filter_with_hvg", "var")
del output_rna.obsm["X_pca"]
del output_rna.varm["pca_loadings"]
del output_rna.obsm["X_umap"]
del output_rna.layers["log_normalized"]
assert_annotation_objects_equal(input_rna, output_rna, promote_precision=True)


print("Test successful!", flush=True)
