import sys
import mudata as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": None,
    "input_var_gene_names": None,
    "input_obs_batch": "sample_id",
    "input_obs_label": None,
    "input_id": "query",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "reference_layer": None,
    "reference_var_gene_names": "ensemblid",
    "reference_obs_batch": "donor_id",
    "reference_obs_label": "cell_ontology_class",
    "reference_id": "reference",
    "output_query": "pbmc_1k_protein_v3_mms_ligned.h5mu",
    "output_reference": "TS_Blood_filtered_aligned.h5mu",
    "output_layer": "_counts",
    "output_var_gene_names": "_gene_names",
    "output_obs_batch": "_sample_id",
    "output_obs_label": "_cell_type",
    "output_obs_id": "_dataset",
    "input_reference_gene_overlap": 100,
    "unkown_celltype_label": "Unknown",
    "overwrite_existing_key": False,
    "output_compression": None,
    "preserve_var_index": False,
    "output_var_index": "_ori_var_index",
    "output_var_common_genes": "_common_vars",
    "align_layers_raw_counts": True,
    "align_layers_lognormalized_counts": False,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression
from cross_check_genes import cross_check_genes

logger = setup_logger()


def copy_layer(adata, input_layer, output_layer, overwrite=False):
    if output_layer not in adata.layers.keys() and input_layer:
        logger.info(f"Copying layer from `{input_layer}` to `{output_layer}`")
        adata.layers[output_layer] = adata.layers[input_layer].copy()

    elif output_layer not in adata.layers.keys() and not input_layer:
        logger.info(f"Copying layer from .X to `{output_layer}`")
        adata.layers[output_layer] = adata.X.copy()

    else:
        if not overwrite:
            raise ValueError(
                f"Layer `{output_layer}` already exists. Data can not be copied."
            )
        logger.warning(f"Layer `{output_layer}` already exists. Overwriting data.")
        adata.layers[output_layer] = (
            adata.layers[input_layer].copy() if input_layer else adata.X.copy()
        )

    return adata


def copy_obs(adata, input_obs_key, output_obs_key, overwrite=False, fill_value=None):
    if output_obs_key not in adata.obs.keys() and input_obs_key:
        logger.info(f"Copying .obs field from `{input_obs_key}` to `{output_obs_key}`")
        adata.obs[output_obs_key] = adata.obs[input_obs_key].copy()

    elif output_obs_key not in adata.obs.keys() and not input_obs_key:
        logger.info(f"Assigning fill value `{fill_value}` to .obs {output_obs_key}")
        adata.obs[output_obs_key] = fill_value

    else:
        if not overwrite:
            raise ValueError(
                f".obs key `{output_obs_key}` already exists. Data can not be copied."
            )

        logger.warning(f".obs key `{output_obs_key}` already exists. Overwriting data.")
        adata.obs[output_obs_key] = (
            adata.obs[input_obs_key].copy() if input_obs_key else fill_value
        )

    return adata


def copy_and_sanitize_var_gene_names(
    adata,
    input_var_key,
    output_var_key,
    overwrite=par["overwrite_existing_key"],
    preserve_index=par["preserve_var_index"],
    var_index_field="ori_var_index",
):
    if output_var_key not in adata.var.keys() and input_var_key:
        logger.info(f"Copying .var field from `{input_var_key}` to `{output_var_key}`")
        adata.var[output_var_key] = adata.var[input_var_key].str.replace(
            "\\.[0-9]+$", "", regex=True
        )

    elif output_var_key not in adata.var.keys() and not input_var_key:
        logger.info(f"Copying .var index to `{output_var_key}`")
        adata.var[output_var_key] = adata.var.index.str.replace(
            "\\.[0-9]+$", "", regex=True
        )

    else:
        if not overwrite:
            raise ValueError(
                f".var key `{output_var_key}` already exists. Data can not be copied."
            )

        logger.warning(f".var key `{output_var_key}` already exists. Overwriting data.")
        adata.var[output_var_key] = (
            adata.var[input_var_key].str.replace("\\.[0-9]+$", "", regex=True)
            if input_var_key
            else adata.var.index.str.replace("\\.[0-9]+$", "", regex=True)
        )

    if not preserve_index:
        logger.info("Replacing .var index with sanitized gene names...")
        adata.var[var_index_field] = adata.var.index
        adata.var.index = list(adata.var[output_var_key])

    return adata


def main():
    logger.info("Reading query and reference data")
    input_adata = mu.read_h5ad(par["input"], mod=par["modality"])
    input_modality = input_adata.copy()

    reference_adata = mu.read_h5ad(par["reference"], mod=par["modality"])
    reference_modality = reference_adata.copy()

    # Aligning layers
    logger.info("### Aligning layers")

    if par["align_layers_lognormalized_counts"] and par["align_layers_raw_counts"]:
        if par["input_layer"] == par["input_layer_lognormalized"]:
            raise ValueError(
                "Layer names for raw and lognormalized counts in the query data can not be identical."
            )

        if par["reference_layer"] == par["reference_layer_lognormalized"]:
            raise ValueError(
                "Layer names for raw and lognormalized counts in the reference data can not be identical."
            )

    if par["align_layers_raw_counts"]:
        logger.info("## Copying query layer raw counts...")
        input_modality = copy_layer(
            input_modality,
            par["input_layer"],
            par["output_layer"],
            overwrite=par["overwrite_existing_key"],
        )
        logger.info("## Copying reference layer raw counts...")
        reference_modality = copy_layer(
            reference_modality,
            par["reference_layer"],
            par["output_layer"],
            overwrite=par["overwrite_existing_key"],
        )

    if par["align_layers_lognormalized_counts"]:
        logger.info("## Copying query layer lognormalized counts...")
        input_modality = copy_layer(
            input_modality,
            par["input_layer_lognormalized"],
            par["output_layer_lognormalized"],
            overwrite=par["overwrite_existing_key"],
        )
        logger.info("## Copying reference layerlognormalized counts...")
        reference_modality = copy_layer(
            reference_modality,
            par["reference_layer_lognormalized"],
            par["output_layer_lognormalized"],
            overwrite=par["overwrite_existing_key"],
        )

    # Aligning batch labels
    logger.info("### Aligning batch labels")
    logger.info("## Copying query .obs batch field...")
    input_modality = copy_obs(
        input_modality,
        par["input_obs_batch"],
        par["output_obs_batch"],
        par["overwrite_existing_key"],
        fill_value=None,
    )
    logger.info("## Copying reference .obs batch field...")
    reference_modality = copy_obs(
        reference_modality,
        par["reference_obs_batch"],
        par["output_obs_batch"],
        overwrite=par["overwrite_existing_key"],
        fill_value=None,
    )

    # Aligning cell type labels
    logger.info("### Aligning cell type labels")
    logger.info("## Copying query .obs cell type label field...")
    input_modality = copy_obs(
        input_modality,
        par["input_obs_label"],
        par["output_obs_label"],
        par["overwrite_existing_key"],
        fill_value=par["unkown_celltype_label"],
    )
    logger.info("## Copying reference .obs cell type label field...")
    reference_modality = copy_obs(
        reference_modality,
        par["reference_obs_label"],
        par["output_obs_label"],
        overwrite=par["overwrite_existing_key"],
        fill_value=par["unkown_celltype_label"],
    )

    # Aligning and sanitizing gene names
    logger.info("### Aligning and sanitizing gene names")
    logger.info("## Copying query .var gene names field...")
    input_modality = copy_and_sanitize_var_gene_names(
        input_modality,
        par["input_var_gene_names"],
        par["output_var_gene_names"],
        overwrite=par["overwrite_existing_key"],
        preserve_index=par["preserve_var_index"],
        var_index_field=par["output_var_index"],
    )
    logger.info("## Copying reference .var gene names field...")
    reference_modality = copy_and_sanitize_var_gene_names(
        reference_modality,
        par["reference_var_gene_names"],
        par["output_var_gene_names"],
        overwrite=par["overwrite_existing_key"],
        preserve_index=par["preserve_var_index"],
        var_index_field=par["output_var_index"],
    )

    # Cross check genes
    logger.info("### Cross checking genes")
    common_vars = cross_check_genes(
        input_modality.var[par["output_var_gene_names"]],
        reference_modality.var[par["output_var_gene_names"]],
        min_gene_overlap=par["input_reference_gene_overlap"],
    )

    # Add common vars to the output
    input_modality.var[par["output_var_common_genes"]] = input_modality.var[
        par["output_var_gene_names"]
    ].isin(common_vars)
    reference_modality.var[par["output_var_common_genes"]] = reference_modality.var[
        par["output_var_gene_names"]
    ].isin(common_vars)

    # Adding an id to the query and reference datasets
    logger.info("### Adding an id to the query and reference datasets")
    logger.info(
        f"## Adding `{par['input_id']}` to the .obs `{par['output_obs_id']}` field of the query dataset..."
    )
    input_modality.obs[par["output_obs_id"]] = par["input_id"]
    logger.info(
        f"## Adding `{par['reference_id']}` to the .obs `{par['output_obs_id']}` field of the reference dataset..."
    )
    reference_modality.obs[par["output_obs_id"]] = par["reference_id"]

    # Writing output data
    logger.info("Writing output data")
    write_h5ad_to_h5mu_with_compression(
        par["output_query"],
        par["input"],
        par["modality"],
        input_modality,
        par["output_compression"],
    )

    write_h5ad_to_h5mu_with_compression(
        par["output_reference"],
        par["reference"],
        par["modality"],
        reference_modality,
        par["output_compression"],
    )


if __name__ == "__main__":
    main()
