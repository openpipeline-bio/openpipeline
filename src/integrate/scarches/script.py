import sys
import mudata
import scvi
from torch.cuda import is_available as cuda_is_available

try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "input_layer": None,
    "input_obs_batch": "sample_id",
    "input_obs_label": None,
    "input_var_gene_names": None,
    "input_obs_categorical_covariate": None,
    "input_obs_continuous_covariate_keys": None,
    "unkown_celltype_label": "Unkown",
    "reference": "TS_blood_reference/scvi",
    "modality": "rna",
    "output": "foo.h5mu",
    "model_output": "./hlca_query_model",
    # Other
    "obsm_output": "X_integrated_scanvi",
    "early_stopping": None,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "max_epochs": 10,
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from set_var_index import set_var_index
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _read_model_name_from_registry(model_path) -> str:
    """Read registry with information about the model, return the model name"""
    registry = scvi.model.base.BaseModelClass.load_registry(model_path)
    return registry["model_name"]


def _detect_base_model(model_path):
    """Read from the model's file which scvi_tools model it contains"""

    names_to_models_map = {
        "AUTOZI": scvi.model.AUTOZI,
        "CondSCVI": scvi.model.CondSCVI,
        "DestVI": scvi.model.DestVI,
        "LinearSCVI": scvi.model.LinearSCVI,
        "PEAKVI": scvi.model.PEAKVI,
        "SCANVI": scvi.model.SCANVI,
        "SCVI": scvi.model.SCVI,
        "TOTALVI": scvi.model.TOTALVI,
        "MULTIVI": scvi.model.MULTIVI,
        "AmortizedLDA": scvi.model.AmortizedLDA,
        "JaxSCVI": scvi.model.JaxSCVI,
    }

    return names_to_models_map[_read_model_name_from_registry(model_path)]


def _align_query_with_registry(adata_query, model, model_path):
    registry = model.load_registry(model_path)["setup_args"]
    base_model = _detect_base_model(model_path)

    # Sanitize gene names and set as index of the AnnData object
    # all scArches VAE models expect gene names to be in the .var index
    adata_query = set_var_index(adata_query, par["input_var_gene_names"])

    # align layers
    if not registry["layer"] == par["input_layer"]:
        if registry["layer"]:
            adata_query.layers[registry["layer"]] = (
                adata_query.layers[par["input_layer"]]
                if par["input_layer"]
                else adata_query.X
            )

        else:
            adata_query.X = (
                adata_query.layers[par["input_layer"]]
                if par["input_layer"]
                else adata_query.X
            )

    # align batch .obs field
    # relevant for AUTOZI, LinearSCVI, PEAKVI, SCANVI, SCVI, TOTALVI, MULTIVI, JaxSCVI
    if registry["batch_key"] and not registry["batch_key"] == par["input_obs_batch"]:
        adata_query.obs[registry["batch_key"]] = adata_query.obs[par["input_obs_batch"]]

    # align celltype .obs field
    # relevant for AUTOZI, CondSCVI, LinearSCVI, PEAKVI, SCANVI, SCVI
    if registry["labels_key"]:
        adata_query.obs[registry["labels_key"]] = (
            adata_query.obs[par["input_obs_label"]]
            if par["input_obs_label"]
            else par["unkown_celltype_label"]
        )

    # align categorical_covariate_kes  .obs field
    # relevant for PEAKVI, SCANVI, SCVI, TOTALVI, MULTIVI

    # align continuous_covariate_keys .obs field
    # relevant for PEAKVI, SCANVI, SCVI, TOTALVI, MULTIVI

    # unkown_key
    # relevant for SCANVI

    # size_factor_key
    # relevant for SCANVI, SCVI, TOTALVI, MULTIVI

    if base_model == "TOTALVI" or base_model == "MULTIVI":
        pass
    # protein_expression_obsm_key
    # relevant for TOTALVI, MULTIVI

    # protein_names_uns_key
    # relevant for TOTALVI, MULTIVI

    return adata_query


def map_to_existing_reference(adata_query, model_path, check_val_every_n_epoch=1):
    """
    A function to map the query data to the reference atlas

    Input:
        * adata_query: An AnnData object with the query
        * model_path: The reference model directory

    Output:
        * vae_query: the trained scvi_tools model
        * adata_query: The AnnData object with the query preprocessed for the mapping to the reference
    """
    model = _detect_base_model(model_path)

    # Keys of the AnnData query object need to match the exact keys in the reference model registry
    adata_query = _align_query_with_registry(adata_query, model, model_path)

    try:
        model.prepare_query_anndata(adata_query, model_path)
    except ValueError:
        logger.warning(
            "ValueError thrown when preparing adata for mapping. Clearing .varm field to prevent it"
        )
        adata_query.varm.clear()
        model.prepare_query_anndata(adata_query, model_path)

    # Load query data into the model
    vae_query = model.load_query_data(adata_query, model_path, freeze_dropout=True)

    # Train scArches model for query mapping
    vae_query.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_monitor=par["early_stopping_monitor"],
        early_stopping_patience=par["early_stopping_patience"],
        early_stopping_min_delta=par["early_stopping_min_delta"],
        check_val_every_n_epoch=check_val_every_n_epoch,
        use_gpu=(cuda_is_available() or mps_is_available()),
    )

    return vae_query, adata_query


def _convert_object_dtypes_to_strings(adata):
    """Convert object dtypes in .var and .obs to string to prevent error when saving file"""

    def convert_cols(df):
        object_cols = df.columns[df.dtypes == "object"]
        for col in object_cols:
            df[col] = df[col].astype(str)
        return df

    adata.var = convert_cols(adata.var)
    adata.obs = convert_cols(adata.obs)

    return adata


def _get_model_path(model_path: str):
    """Obtain path to the directory with reference model. If the proposed `model_path` is a .zip archive, unzip it. If nesessary, convert model to the new format

    Parameters
    ----------
    model_path : str
        Path to a directory, where to search for the model or to a zip file containing the model

    Returns
    -------
    Path to a directory with reference model in format of scvi-tools>=0.15
    """
    import os
    import zipfile
    import tempfile
    from pathlib import Path

    if os.path.isdir(model_path) and "model.pt" in os.listdir(model_path):
        # Probably, the `model_path` already contains model in the output format of scvi-tools>=0.15
        return model_path

    # The model either has old format or is a zip file downloaded from Zenodo
    new_directory = Path(tempfile.TemporaryDirectory().name)

    if zipfile.is_zipfile(model_path):
        with zipfile.ZipFile(model_path) as archive:
            archive.extractall(new_directory)
            model_dir = next(new_directory.glob("**/*.pt")).parent

    else:
        model_dir = next(Path(model_path).glob("**/*.pt")).parent

    if "model_params.pt" in os.listdir(model_dir):
        # The model is in the `directory`, but it was generated with scvi-tools<0.15
        # TODO: for new references (that could not be SCANVI based), we need to check the base class somehow. Reading registry does not work with models generated by scvi-tools<0.15
        # Here I assume that the reference model is for HLCA and thus is SCANVI based
        converted_model_path = os.path.join(model_dir, "converted")
        scvi.model.SCANVI.convert_legacy_save(model_dir, converted_model_path)
        return converted_model_path

    elif "model.pt" in os.listdir(model_dir):
        # Archive contained model in the new format, so just return the directory
        return model_dir

    else:
        raise ValueError(
            "Cannot find model in the provided reference path. Please, provide a path or a link to the directory with reference model. For HLCA use https://zenodo.org/record/6337966/files/HLCA_reference_model.zip"
        )


def main():
    logger.info("Reading %s, modality %s", par["input"], par["modality"])
    adata = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])
    adata_query = adata.copy()

    model_path = _get_model_path(par["reference"])
    vae_query, adata_query = map_to_existing_reference(
        adata_query, model_path=model_path
    )
    model_name = _read_model_name_from_registry(model_path)

    # Save info about the used model
    adata.uns["integration_method"] = model_name

    logger.info("Trying to write latent representation")
    output_key = par["obsm_output"].format(model_name=model_name)
    adata.obsm[output_key] = vae_query.get_latent_representation()

    logger.info("Converting dtypes")
    adata = _convert_object_dtypes_to_strings(adata)

    logger.info(
        "Saving h5mu file to %s with compression %s",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )

    logger.info("Saving model")
    vae_query.save(par["model_output"], overwrite=True)


if __name__ == "__main__":
    main()
