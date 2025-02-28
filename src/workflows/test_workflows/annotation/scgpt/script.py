from mudata import read_h5mu
import shutil
import os
import sys
from pathlib import Path
import pytest

##VIASH START
par = {"input": "input.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_obsm = [
        "gene_id_tokens",
        "values_tokenized",
        "padding_mask",
        "bin_edges",
        "binned_counts",
    ]
    expected_var = ["scgpt_filter_with_hvg", "scgpt_cross_checked_genes"]
    expected_obs = ["scgpt_pred", "scgpt_probability"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."

    assert all(
        key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm
    ), f"Input mod['rna'] obs columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    assert all(
        key in list(input_mudata.mod["rna"].var) for key in expected_var
    ), f"Input mod['rna'] var columns should be: {expected_var}, found: {input_mudata.mod['rna'].var.keys()}."
    assert all(
        key in list(input_mudata.mod["rna"].obs) for key in expected_obs
    ), f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    # hvg subsetting is not exact - add 10% to allowed data shape
    assert (
        input_mudata.mod["rna"].obsm["binned_counts"].shape[1]
        <= par["n_hvg"] + 0.1 * par["n_hvg"]
    ), f"Input shape should be lower or equal than --n_hvg {par['n_hvg']}, found: {input_mudata.shape[1]}."


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    from importlib import resources
    shutil.copyfile(
        resources.files("openpipeline_testutils").joinpath("conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
