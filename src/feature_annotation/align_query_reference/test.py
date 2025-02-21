import sys
import os
import pytest
import subprocess
import re
import numpy as np
import mudata as mu

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"


@pytest.fixture
def copy_layer(random_h5mu_path):
    def wrapper(input_mudata_file, output_layer_name):
        input_mudata = mu.read_h5mu(input_mudata_file)
        input_adata = input_mudata.mod["rna"]
        input_adata.layers[output_layer_name] = input_adata.X.copy()
        output_mudata_file = random_h5mu_path()
        input_mudata.write_h5mu(output_mudata_file)
        return output_mudata_file

    return wrapper


@pytest.fixture
def add_var(random_h5mu_path):
    def wrapper(input_mudata_file, var_field):
        input_mudata = mu.read_h5mu(input_mudata_file)
        input_adata = input_mudata.mod["rna"]
        input_adata.var[var_field] = "fill_value"
        output_mudata_file = random_h5mu_path()
        input_mudata.write_h5mu(output_mudata_file)
        return output_mudata_file

    return wrapper


@pytest.fixture
def add_obs(random_h5mu_path):
    def wrapper(input_mudata_file, obs_field):
        input_mudata = mu.read_h5mu(input_mudata_file)
        input_adata = input_mudata.mod["rna"]
        input_adata.obs[obs_field] = "fill_value"
        output_mudata_file = random_h5mu_path()
        input_mudata.write_h5mu(output_mudata_file)
        return output_mudata_file

    return wrapper


def test_simple_execution(run_component, random_h5mu_path):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--input_obs_batch",
            "sample_id",
            "--reference",
            reference_file,
            "--reference_obs_batch",
            "donor_id",
            "--reference_obs_label",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--output_query",
            output_query_file,
            "--output_reference",
            output_reference_file,
        ]
    )

    assert os.path.exists(output_query_file), "Output query file does not exist"
    assert os.path.exists(output_reference_file), "Output reference file does not exist"

    query_adata = mu.read_h5mu(output_query_file).mod["rna"]
    reference_adata = mu.read_h5mu(output_reference_file).mod["rna"]

    expected_layers = ["_counts"]
    expected_var = ["_gene_names", "_ori_var_index", "_common_vars"]
    expected_obs = ["_sample_id", "_cell_type", "_dataset"]

    # Evaluate presence of obs, var, layers
    assert all(
        key in list(query_adata.obs) for key in expected_obs
    ), f"Query obs columns should be: {expected_obs}, found: {query_adata.obs.keys()}."

    assert all(
        key in list(query_adata.var) for key in expected_var
    ), f"Query var columns should be: {expected_var}, found: {query_adata.var.keys()}."
    assert all(
        key in list(query_adata.layers) for key in expected_layers
    ), f"Query layers should be: {expected_layers}, found: {query_adata.layers.keys()}."

    assert all(
        key in list(reference_adata.obs) for key in expected_obs
    ), f"Reference obs columns should be: {expected_obs}, found: {reference_adata.obs.keys()}."

    assert all(
        key in list(reference_adata.var) for key in expected_var
    ), f"Reference var columns should be: {expected_var}, found: {reference_adata.var.keys()}."
    assert all(
        key in list(reference_adata.layers) for key in expected_layers
    ), f"Reference layers should be: {expected_layers}, found: {reference_adata.layers.keys()}."

    # Evaluate values in obs, var, layers
    assert np.all(
        query_adata.var["_gene_names"] == query_adata.var.index
    ), "Query .var _gene_names should be equal to query .var index"
    assert np.all(
        query_adata.obs["_sample_id"] == query_adata.obs["sample_id"]
    ), "Query .obs _sample_id should be equal to query .obs sample_id"
    assert np.all(
        query_adata.obs["_cell_type"] == "Unknown"
    ), "Query .obs _cell_type should have value Unkown"
    assert np.all(
        query_adata.obs["_dataset"] == "query"
    ), "Query .obs _dataset should have value query"

    assert np.all(
        reference_adata.var["_gene_names"] == reference_adata.var.index
    ), "Reference .var _gene_names should be equal to query .var index"

    assert np.all(
        reference_adata.var.index != reference_adata.var["ensemblid"]
    ), "Reference .var index should not be equal to reference .var ensemblid"

    assert np.all(
        reference_adata.var.index
        == [re.sub("\\.[0-9]+$", "", s) for s in reference_adata.var["ensemblid"]]
    ), "Reference .var index should be equal to sanitized reference .var ensemblid "

    assert np.any(
        reference_adata.var["_ori_var_index"] != reference_adata.var["ensemblid"]
    ), "Reference .var _ori_var_index should not be all equal to reference .var ensemblid"

    assert np.all(
        reference_adata.var["_ori_var_index"] != reference_adata.var.index
    ), "Reference .var _ori_var_index should be equal to reference .var index"

    assert np.all(
        reference_adata.obs["_sample_id"] == reference_adata.obs["donor_id"]
    ), "Reference .obs _sample_id should be equal to query .obs sample_id"
    assert np.all(
        reference_adata.obs["_cell_type"] == reference_adata.obs["cell_ontology_class"]
    ), "Reference .obs _cell_type should be equal to reference .obs cell_ontology_class"
    assert np.all(
        reference_adata.obs["_dataset"] == "reference"
    ), "Reference .obs _dataset should have value reference"


def test_align_multiple_layers(run_component, random_h5mu_path):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--align_layers_lognormalized_counts",
            "True",
            "--input_layer_lognormalized",
            "log_normalized",
            "--input_obs_batch",
            "sample_id",
            "--reference",
            reference_file,
            "--reference_layer_lognormalized",
            "log_normalized",
            "--reference_obs_batch",
            "donor_id",
            "--reference_obs_label",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--output_query",
            output_query_file,
            "--output_reference",
            output_reference_file,
        ]
    )

    query_adata = mu.read_h5mu(output_query_file).mod["rna"]
    reference_adata = mu.read_h5mu(output_reference_file).mod["rna"]

    expected_layers = ["_counts", "_log_normalized"]

    assert all(
        key in list(query_adata.layers) for key in expected_layers
    ), f"Query layers should be: {expected_layers}, found: {query_adata.layers.keys()}."
    assert all(
        key in list(reference_adata.layers) for key in expected_layers
    ), f"Reference layers should be: {expected_layers}, found: {reference_adata.layers.keys()}."

    args_identical_query_layers = [
        "--input",
        input_file,
        "--modality",
        "rna",
        "--align_layers_lognormalized_counts",
        "True",
        "--reference_layer_lognormalized",
        "--log_normalized",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_identical_query_layers)
    assert re.search(
        r"Layer names for raw and lognormalized counts in the query data can not be identical.",
        err.value.stdout.decode("utf-8"),
    )
    args_identical_query_layers = [
        "--input",
        input_file,
        "--modality",
        "rna",
        "--align_layers_lognormalized_counts",
        "True",
        "--input_layer_lognormalized",
        "--log_normalized",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_identical_query_layers)
    assert re.search(
        r"Layer names for raw and lognormalized counts in the reference data can not be identical.",
        err.value.stdout.decode("utf-8"),
    )


def test_copy_layer(run_component, random_h5mu_path, copy_layer):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()
    input_file_with_layer = copy_layer(input_file, "copied_counts")

    run_component(
        [
            "--input",
            input_file_with_layer,
            "--input_layer",
            "copied_counts",
            "--modality",
            "rna",
            "--input_obs_batch",
            "sample_id",
            "--reference",
            reference_file,
            "--reference_obs_batch",
            "donor_id",
            "--reference_obs_label",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--output_query",
            output_query_file,
            "--output_reference",
            output_reference_file,
        ]
    )

    query_adata = mu.read_h5mu(output_query_file).mod["rna"]

    assert (
        "_counts" in query_adata.layers.keys()
    ), f"_counts layer should be present, found: {query_adata.layers.keys()}."

    args_existing_layer = [
        "--input",
        input_file_with_layer,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_layer",
        "copied_counts",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_existing_layer)
    assert re.search(
        r"ValueError: Layer `copied_counts` already exists. Data can not be copied.",
        err.value.stdout.decode("utf-8"),
    )

    disable_raise_args = [
        "--input",
        input_file_with_layer,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_layer",
        "copied_counts",
        "--overwrite_existing_key",
        "True",
    ]

    run_component(disable_raise_args)


def test_overwrite_obs(run_component, random_h5mu_path, add_obs):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()

    input_file_with_obs = add_obs(input_file, "Obs")

    args = [
        "--input",
        input_file_with_obs,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_obs_batch",
        "Obs",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: .obs key `Obs` already exists. Data can not be copied.",
        err.value.stdout.decode("utf-8"),
    )

    disable_raise_args = [
        "--input",
        input_file,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_obs_batch",
        "Obs" "--overwrite_existing_key",
        "True",
    ]

    run_component(disable_raise_args)


def test_overwrite_var(run_component, random_h5mu_path, add_var):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()

    input_file_with_obs = add_var(input_file, "Var")

    args = [
        "--input",
        input_file_with_obs,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_var_gene_names",
        "Var",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: .var key `Var` already exists. Data can not be copied.",
        err.value.stdout.decode("utf-8"),
    )

    disable_raise_args = [
        "--input",
        input_file,
        "--modality",
        "rna",
        "--input_obs_batch",
        "sample_id",
        "--reference",
        reference_file,
        "--reference_obs_batch",
        "donor_id",
        "--reference_obs_label",
        "cell_ontology_class",
        "--reference_var_gene_names",
        "ensemblid",
        "--output_query",
        output_query_file,
        "--output_reference",
        output_reference_file,
        "--output_var_gene_names",
        "Var",
        "--overwrite_existing_key",
        "True",
    ]

    run_component(disable_raise_args)


def test_preserve_var_index(run_component, random_h5mu_path):
    output_query_file = random_h5mu_path()
    output_reference_file = random_h5mu_path()
    # input_file_transformed = copy_layer(input_file, "copied_counts")

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--input_obs_batch",
            "sample_id",
            "--reference",
            reference_file,
            "--reference_obs_batch",
            "donor_id",
            "--reference_obs_label",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--output_query",
            output_query_file,
            "--output_reference",
            output_reference_file,
            "--preserve_var_index",
        ]
    )

    ori_query_adata = mu.read_h5mu(input_file).mod["rna"]
    ori_reference_adata = mu.read_h5mu(reference_file).mod["rna"]
    query_adata = mu.read_h5mu(output_query_file).mod["rna"]
    reference_adata = mu.read_h5mu(output_reference_file).mod["rna"]

    assert "_ori_var_index" not in query_adata.var.keys()
    assert "_ori_var_index" not in reference_adata.var.keys()

    assert np.all(
        ori_query_adata.var.index == query_adata.var.index
    ), "Query .var index should be preserved"

    assert np.all(
        ori_reference_adata.var.index == reference_adata.var.index
    ), "Reference .var index should be preserved"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
