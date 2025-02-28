import os
import subprocess
import scanpy as sc
import mudata as mu
import sys
import pytest
import re
import pandas as pd

from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "resources_dir": "resources_test/",
    "config": "./src/feature_annotation/highly_variable_features_scanpy/config.vsh.yaml",
    "executable": "./target/executable/feature_annotation/highly_variable_features_scanpy/highly_variable_features_scanpy",
}
## VIASH END

sys.path.append(meta["resources_dir"])


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_data(input_path):
    mu_in = mu.read_h5mu(input_path)
    return mu_in


@pytest.fixture
def lognormed_test_data(input_data):
    rna_in = input_data.mod["rna"]
    assert "filter_with_hvg" not in rna_in.var.columns
    log_transformed = sc.pp.log1p(rna_in, copy=True)
    rna_in.layers["log_transformed"] = log_transformed.X
    rna_in.uns["log1p"] = log_transformed.uns["log1p"]
    return input_data


@pytest.fixture
def lognormed_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed.h5mu"
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture
def lognormed_batch_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed_batch.h5mu"
    rna_mod = lognormed_test_data.mod["rna"]
    rna_mod.obs["batch"] = "A"
    column_index = rna_mod.obs.columns.get_indexer(["batch"])
    rna_mod.obs.iloc[slice(rna_mod.n_obs // 2, None), column_index] = "B"
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture()
def filter_data_path(tmp_path, input_data):
    temp_h5mu = tmp_path / "filtered.h5mu"
    rna_in = input_data.mod["rna"]
    sc.pp.filter_genes(rna_in, min_counts=20)
    input_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture()
def common_vars_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed_var_input.h5mu"
    rna_in = lognormed_test_data.mod["rna"]
    rna_in.var["common_vars"] = False
    rna_in.var["common_vars"].iloc[:10000] = True
    rna_in.var["common_vars"] = rna_in.var["common_vars"].astype("boolean")
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


def test_filter_with_hvg(run_component, lognormed_test_data_path):
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_test_data_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del data["rna"].varm["hvg"]
    assert_annotation_objects_equal(lognormed_test_data_path, data)


def test_filter_with_hvg_var_input(run_component, common_vars_data_path):
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            common_vars_data_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
        ]
    )

    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            common_vars_data_path,
            "--output",
            "common_vars.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
            "--var_input",
            "common_vars",
        ]
    )

    mdata = mu.read_h5mu("output.h5mu")
    common_vars = mu.read_h5mu("common_vars.h5mu")

    # Assert detected HVG are different
    hvg = mdata.mod["rna"][:, mdata.mod["rna"].var["filter_with_hvg"]].var_names
    common_vars_hvg = common_vars.mod["rna"][
        :, common_vars.mod["rna"].var["filter_with_hvg"]
    ].var_names

    assert len(hvg) != len(
        common_vars_hvg
    ), "Number of HVG should be different when var_input is defined"

    # Assert original data is unchanged
    # Put the output data back into its original shape
    # so that we can compare it to the input
    mdata.mod["rna"].var = mdata.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    mdata.var = mdata.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del mdata["rna"].varm["hvg"]

    common_vars.mod["rna"].var = common_vars.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    common_vars.var = common_vars.var.drop(
        columns=["rna:filter_with_hvg"], errors="raise"
    )
    del common_vars["rna"].varm["hvg"]

    assert_annotation_objects_equal(mdata, common_vars)


def test_filter_with_hvg_batch_with_batch(
    run_component, lognormed_batch_test_data_path
):
    """
    Make sure that selecting a layer works together with obs_batch_key.
    https://github.com/scverse/scanpy/issues/2396
    """
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_batch_test_data_path,
            "--output",
            "output.h5mu",
            "--obs_batch_key",
            "batch",
            "--layer",
            "log_transformed",
        ]
    )
    assert os.path.exists("output.h5mu")
    output_data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in output_data.mod["rna"].var.columns

    # Check the contents of the output to check if the correct layer was selected
    input_mudata = mu.read_h5mu(lognormed_batch_test_data_path)
    input_data = input_mudata.mod["rna"].copy()
    input_data.X = input_data.layers["log_transformed"].copy()
    del input_data.layers["log_transformed"]
    input_data.uns["log1p"]["base"] = None
    expected_output = sc.pp.highly_variable_genes(
        input_data, batch_key="batch", inplace=False, subset=False
    )
    expected_output = expected_output.reindex(index=input_mudata.mod["rna"].var.index)
    pd.testing.assert_series_equal(
        expected_output["highly_variable"],
        output_data.mod["rna"].var["filter_with_hvg"],
        check_names=False,
    )
    # Put the output data back into its original shape
    # so that we can compare it to the input
    output_data.mod["rna"].var = output_data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    output_data.var = output_data.var.drop(
        columns=["rna:filter_with_hvg"], errors="raise"
    )
    assert_annotation_objects_equal(lognormed_batch_test_data_path, output_data)


def test_filter_with_hvg_seurat_v3_requires_n_top_features(run_component, input_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--flavor",
                "seurat_v3",  # Uses raw data.
                "--output",
                "output.h5mu",
            ]
        )
    assert re.search(
        "When flavor is set to 'seurat_v3', you are required to set 'n_top_features'.",
        err.value.stdout.decode("utf-8"),
    )


def test_filter_with_hvg_seurat_v3(run_component, input_path):
    run_component(
        [
            "--input",
            input_path,
            "--flavor",
            "seurat_v3",  # Uses raw data.
            "--output",
            "output.h5mu",
            "--n_top_features",
            "50",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    assert_annotation_objects_equal(input_path, data)


def test_filter_with_hvg_cell_ranger(run_component, filter_data_path):
    run_component(
        [
            "--input",
            filter_data_path,
            "--flavor",
            "cell_ranger",  # Must use filtered data.
            "--output",
            "output.h5mu",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del data["rna"].varm["hvg"]
    assert_annotation_objects_equal(filter_data_path, data)


def test_filter_with_hvg_cell_ranger_unfiltered_data_change_error_message(
    run_component, input_path
):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--flavor",
                "cell_ranger",  # Must use filtered data, but in this test we use unfiltered data
                "--output",
                "output.h5mu",
            ]
        )
    assert re.search(
        r"Scanpy failed to calculate hvg. The error "
        r"returned by scanpy \(see above\) could be the "
        r"result from trying to use this component on unfiltered data.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
