from mudata import read_h5mu
import pytest

# from openpipelinetest_utils.asserters import assert_annotation_objects_equal
import sys
import os
import shutil
from pathlib import Path

##VIASH START
par = {
    "input": "work/f7/bb4e1486addee1b3b8938616939f43/_viash_par/input_1/no_leiden_resolutions_test.knn.output.h5mu",
    "orig_input": "",
}

meta = {"resources_dir": "src/base/openpipelinetestutils"}

##VIASH END
sys.path.append(meta["resources_dir"])
from openpipelinetest_utils.asserters import assert_annotation_objects_equal


@pytest.fixture
def input_mudata():
    return read_h5mu(par["input"])


@pytest.fixture
def orig_input_mudata():
    return read_h5mu(par["orig_input"])


def test_mdata_fields(expected_fields, fields):
    assert all(
        key in fields for key in expected_fields
    ), f"Expected fields are: {expected_fields}, found: {fields}"


def test_run(input_mudata, orig_input_mudata):
    expected_obsm = [
        "gene_id_tokens",
        "values_tokenized",
        "padding_mask",
        "bin_edges",
        "X_integrated_scgpt",
        "binned_counts",
    ]
    expected_var = [
        "scgpt_filter_with_hvg",
        "scgpt_cross_checked_genes",
        "gene_symbols",
    ]
    expected_obs = ["cell_type_pred", "cell_type_probability", "batch_label"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."

    test_mdata_fields(expected_obsm, list(input_mudata.mod["rna"].obsm.keys()))
    test_mdata_fields(expected_var, list(input_mudata.mod["rna"].var.keys()))
    test_mdata_fields(expected_obs, list(input_mudata.mod["rna"].obs.keys()))

    # hvg subsetting is not exact - add 2% to allowed data shape
    assert (
        input_mudata.obsm["binned_counts"].shape[1]
        <= par["n_hvg"] + 0.02 * par["n_hvg"]
    ), f"Input shape should be lower or equal than --n_hvg {par['n_hvg']}, found: {input_mudata.shape[1]}."

    assert_annotation_objects_equal(
        input_mudata.mod["rna"], orig_input_mudata.mod["rna"], promote_precision=True
    )


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    from importlib import resources
    shutil.copyfile(
        resources.files("openpipeline_testutils").joinpath("conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
