import mudata as mu
import tempfile
import subprocess
import os
import sys
import numpy as np
from scipy.sparse import csr_matrix
from cellbender.remove_background.downstream import anndata_from_h5
## VIASH START
file_input = "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

par = {
    # inputs
    "input": file_input,
    "modality": "rna",
    # outputs
    "output": "output.h5mu",
    "layer_output": "corrected",
    "obs_background_fraction": "background_fraction",
    "obs_cell_probability": "cell_probability",
    "obs_cell_size": "cell_size",
    "obs_droplet_efficiency": "droplet_efficiency",
    "obs_latent_scale": "latent_scale",
    "var_ambient_expression": "ambient_expression",
    "obsm_gene_expression_encoding": "gene_expression_encoding",
    # args
    "expected_cells_from_qc": False,
    "expected_cells": 1000,
    "total_droplets_included": 25000,
    "force_cell_umi_prior": None,
    "force_empty_umi_prior": None,
    "model": "full",
    "epochs": 150,
    "low_count_threshold": 5,
    "z_dim": 64,
    "z_layers": [512],
    "training_fraction": 0.9,
    "empty_drop_training_fraction": 0.2,
    "ignore_features": [],
    "fpr": [0.01],
    "exclude_feature_types": [],
    "projected_ambient_count_threshold": 0.1,
    "learning_rate": 1.0E-4,
    "final_elbo_fail_fraction": None,
    "epoch_elbo_fail_fraction": None,
    "num_training_tries": 1,
    "learning_rate_retry_mult": 0.2,
    "posterior_batch_size": 128,
    "posterior_regulation": None,
    "alpha": None,
    "q": None,
    "estimator": "mckp",
    "estimator_multiple_cpu": False,
    "constant_learning_rate": True,
    "debug": False,
    "cuda": False
}
meta = {
    "temp_dir": os.getenv("VIASH_TEMP"),
    "resources_dir": "src/correction/cellbender_remove_background"
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()


logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Performing log transformation on modality %s", mod)
data = mdata.mod[mod]

# import pathlib
# with pathlib.Path(os.path.dirname(par["output"])) / "cellbender" as temp_dir:
#     os.mkdir(temp_dir)
with tempfile.TemporaryDirectory(prefix="cellbender-", dir=meta["temp_dir"]) as temp_dir:
    # construct paths within tempdir
    input_file = os.path.join(temp_dir, "input.h5ad")
    output_file = os.path.join(temp_dir, "output.h5")

    logger.info("Creating AnnData input file for CellBender: '%s'", input_file)
    data.write_h5ad(input_file)

    logger.info("Constructing CellBender command")
    cmd_pars = [
        "cellbender", "remove-background",
        "--input", input_file,
        "--output", output_file,
        # don't create checkpoints because they're not used / returned anyways
        "--checkpoint-mins", "99999999"
    ]

    if meta.get("cpus") is not None:
        cmd_pars += ["--cpu-threads", str(meta["cpus"])]

    extra_args = [
        ("--expected-cells", "expected_cells", True),
        ("--total-droplets-included", "total_droplets_included", True),
        ("--force-cell-umi-prior", "force_cell_umi_prior", True),
        ("--force-empty-umi-prior", "force_empty_umi_prior", True),
        ("--model", "model", True),
        ("--epochs", "epochs", True),
        ("--low-count-threshold", "low_count_threshold", True),
        ("--z-dim", "z_dim", True),
        ("--z-layers", "z_layers", True),
        ("--training-fraction", "training_fraction", True),
        ("--empty-drop-training-fraction", "empty_drop_training_fraction", True),
        ("--ignore-features", "ignore_features", True),
        ("--fpr", "fpr", True),
        ("--exclude-feature-types", "exclude_feature_types", True),
        ("--projected-ambient-count-threshold", "projected_ambient_count_threshold", True),
        ("--learning-rate", "learning_rate", True),
        ("--final-elbo-fail-fraction", "final_elbo_fail_fraction", True),
        ("--epoch-elbo-fail-fraction", "epoch_elbo_fail_fraction", True),
        ("--num-training-tries", "num_training_tries", True),
        ("--learning-rate-retry-mult", "learning_rate_retry_mult", True),
        ("--posterior-batch-size", "posterior_batch_size", True),
        ("--posterior-regulation", "posterior_regulation", True),
        ("--alpha", "alpha", True),
        ("--q", "q", True),
        ("--estimator", "estimator", True),
        ("--estimator-multiple-cpu", "estimator_multiple_cpu", False),
        ("--constant-learning-rate", "constant_learning_rate", False),
        ("--debug", "debug", False),
        ("--cuda", "cuda", False),
    ]
    for (flag, name, is_kwarg) in extra_args:
        if par[name]:
            values = par[name] if isinstance(par[name], list) else [par[name]]
            cmd_pars += [flag] + [str(val) for val in values] if is_kwarg else [flag]

    if par["expected_cells_from_qc"] and "metrics_cellranger" in data.uns:
        assert par["expected_cells"] is None, "If min_counts is defined, expected_cells should be undefined"
        assert par["total_droplets_included"] is None, "If min_counts is defined, expected_cells should be undefined"
        met = data.uns["metrics_cellranger"]
        col_name = "Estimated Number of Cells"
        assert col_name in met.columns, "%s should be a column in .obs[metrics_cellranger]"
        est_cells = met[col_name].values[0]
        logger.info("Selecting --expected-cells %d and --total-droplets-included %d", est_cells, est_cells * 5)
        cmd_pars += ["--expected-cells", str(est_cells), "--total-droplets-included", str(5*est_cells)]

    logger.info("Running CellBender: '%s'", ' '.join(cmd_pars))
    out = subprocess.check_output(cmd_pars).decode("utf-8")

    logger.info("Reading CellBender 10xh5 output file: '%s'", output_file)
    adata_out = anndata_from_h5(output_file, analyzed_barcodes_only=False)

    logger.info("CellBender output format:", adata_out)

    # AnnData object with n_obs x n_vars = 6794880 x 33538
    #     obs: 'cellbender_analyzed'
    #     var: 'ambient_expression', 'feature_type', 'genome', 'gene_id', 'cellbender_analyzed'
    #     uns: 'background_fraction', 'barcode_indices_for_latents', 'cell_probability', 'cell_size', 'droplet_efficiency', 'gene_expression_encoding', 
    #          'cell_size_lognormal_std', 'empty_droplet_size_lognormal_loc', 'empty_droplet_size_lognormal_scale', 'swapping_fraction_dist_params', 
    #          'barcodes_analyzed', 'barcodes_analyzed_inds', 'estimator', 'features_analyzed_inds', 'fraction_data_used_for_testing', 'learning_curve_learning_rate_epoch', 
    #          'learning_curve_learning_rate_value', 'learning_curve_test_elbo', 'learning_curve_test_epoch', 'learning_curve_train_elbo', 'learning_curve_train_epoch', 
    #          'target_false_positive_rate'

    logger.info("Copying X output to MuData")
    data.layers[par["layer_output"]] = adata_out.X

    logger.info("Copying .obs output to MuData")
    obs_store = {
        "obs_background_fraction": "background_fraction",
        "obs_cell_probability": "cell_probability",
        "obs_cell_size": "cell_size",
        "obs_droplet_efficiency": "droplet_efficiency",
        "obs_latent_scale": "latent_scale"
    }
    for to_name, from_name in obs_store.items():
        if par[to_name]:
            if from_name in adata_out.obs:
                data.obs[par[to_name]] = adata_out.obs[from_name]
            # when using unfiltered data, the values will be in uns instead of obs
            elif from_name in adata_out.uns and "barcode_indices_for_latents" in adata_out.uns:
                vec = np.zeros(data.n_obs)
                vec[adata_out.uns["barcode_indices_for_latents"]] = adata_out.uns[from_name]
                data.obs[par[to_name]] = vec

    logger.info("Copying .var output to MuData")
    var_store = { "var_ambient_expression": "ambient_expression" }
    for to_name, from_name in var_store.items():
        if par[to_name]:
            data.var[par[to_name]] = adata_out.var[from_name]

    logger.info("Copying obsm_gene_expression_encoding output to MuData")
    obsm_store = { "obsm_gene_expression_encoding": "gene_expression_encoding" }
    for to_name, from_name in obsm_store.items():
        if par[to_name]:
            if from_name in adata_out.obsm:
                 data.obsm[par[to_name]] = adata_out.obsm[from_name]
            elif from_name in adata_out.uns and "barcode_indices_for_latents" in adata_out.uns:
                matrix_to_store = adata_out.uns[from_name]
                number_of_obs = data.X.shape[0]
                latent_space_sparse = csr_matrix((number_of_obs, par["z_dim"]),
                                                 dtype=adata_out.uns[from_name].dtype)
                obs_rows_in_space_representation = adata_out.uns["barcode_indices_for_latents"]
                latent_space_sparse[obs_rows_in_space_representation] = adata_out.uns[from_name]
                data.obsm[par[to_name]] = latent_space_sparse
            else:
                raise RuntimeError("Requested to save latent gene encoding, but the data is either missing "
                                   "from cellbender output or in an incorrect format.")


logger.info("Writing to file %s", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
