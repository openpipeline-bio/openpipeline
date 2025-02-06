from scanpy._utils import check_nonnegative_integers
import mudata
import scvi

### VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "input_layer": None,
    "var_gene_names": "ensemblid",
    "obs_batch": "donor_id",
    "obs_labels": None,
    "obs_size_factor": None,
    "obs_categorical_covariate": None,
    "obs_continuous_covariate": None,
    "var_input": None,
    "output": "foo.h5mu",
    "obsm_output": "X_scvi_integrated",
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "reduce_lr_on_plateau": True,
    "lr_factor": 0.6,
    "lr_patience": 30,
    "max_epochs": 10,
    "n_obs_min_count": 10,
    "n_var_min_count": 10,
    "output_model": "TS_blood_reference",
    "output_compression": "gzip",
    "n_hidden_nodes": 128,
    "n_dimensions_latent_space": 30,
    "n_hidden_layers": 2,
    "dropout_rate": 0.1,
    "dispersion": "gene",
    "gene_likelihood": "nb",
    "use_layer_normalization": "both",
    "use_batch_normalization": "none",
    "encode_covariates": True,
    "deeply_inject_covariates": False,
    "use_observed_lib_size": False,
}

meta = {"resources_dir": "src/utils"}
### VIASH END

import sys

sys.path.append(meta["resources_dir"])

from subset_vars import subset_vars
from set_var_index import set_var_index
from compress_h5mu import write_h5ad_to_h5mu_with_compression


# TODO: optionally, move to qa
# https://github.com/openpipelines-bio/openpipeline/issues/435
def check_validity_anndata(adata, layer, obs_batch, n_obs_min_count, n_var_min_count):
    assert check_nonnegative_integers(
        adata.layers[layer] if layer else adata.X
    ), "Make sure input adata contains raw_counts"

    assert len(set(adata.var_names)) == len(
        adata.var_names
    ), "Dataset contains multiple genes with same gene name."

    # Ensure every obs_batch category has sufficient observations
    assert (
        min(adata.obs[[obs_batch]].value_counts()) > n_obs_min_count
    ), f"Anndata has fewer than {n_obs_min_count} cells."

    assert (
        adata.n_vars > n_var_min_count
    ), f"Anndata has fewer than {n_var_min_count} genes."


def main():
    adata = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])

    if par["var_input"]:
        # Subset to HVG
        adata_subset = subset_vars(adata, subset_col=par["var_input"]).copy()
    else:
        adata_subset = adata.copy()

    # Sanitize gene names and set as index of the AnnData object
    adata_subset = set_var_index(adata_subset, par["var_gene_names"])

    check_validity_anndata(
        adata_subset,
        par["input_layer"],
        par["obs_batch"],
        par["n_obs_min_count"],
        par["n_var_min_count"],
    )
    # Set up the data
    scvi.model.SCVI.setup_anndata(
        adata_subset,
        batch_key=par["obs_batch"],
        layer=par["input_layer"],
        labels_key=par["obs_labels"],
        size_factor_key=par["obs_size_factor"],
        categorical_covariate_keys=par["obs_categorical_covariate"],
        continuous_covariate_keys=par["obs_continuous_covariate"],
    )

    # Set up the model
    vae_uns = scvi.model.SCVI(
        adata_subset,
        n_hidden=par["n_hidden_nodes"],
        n_latent=par["n_dimensions_latent_space"],
        n_layers=par["n_hidden_layers"],
        dropout_rate=par["dropout_rate"],
        dispersion=par["dispersion"],
        gene_likelihood=par["gene_likelihood"],
        use_layer_norm=par["use_layer_normalization"],
        use_batch_norm=par["use_batch_normalization"],
        encode_covariates=par[
            "encode_covariates"
        ],  # Default (True) is for better scArches performance -> maybe don't use this always?
        deeply_inject_covariates=par[
            "deeply_inject_covariates"
        ],  # Default (False) for better scArches performance -> maybe don't use this always?
        use_observed_lib_size=par[
            "use_observed_lib_size"
        ],  # When size_factors are not passed
    )

    plan_kwargs = {
        "reduce_lr_on_plateau": par["reduce_lr_on_plateau"],
        "lr_patience": par["lr_patience"],
        "lr_factor": par["lr_factor"],
    }

    # Train the model
    vae_uns.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_monitor=par["early_stopping_monitor"],
        early_stopping_patience=par["early_stopping_patience"],
        early_stopping_min_delta=par["early_stopping_min_delta"],
        plan_kwargs=plan_kwargs,
        check_val_every_n_epoch=1,
        accelerator="auto",
    )
    # Note: train_size=1.0 should give better results, but then can't do early_stopping on validation set

    # Get the latent output
    adata.obsm[par["obsm_output"]] = vae_uns.get_latent_representation()

    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )
    if par["output_model"]:
        vae_uns.save(par["output_model"], overwrite=True)


if __name__ == "__main__":
    main()
