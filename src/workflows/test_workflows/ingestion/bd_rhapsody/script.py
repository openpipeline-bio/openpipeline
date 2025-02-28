from mudata import read_h5mu
import numpy as np
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
    expected_var = ["gene_name", "feature_type", "reference_file", "gene_ids"]
    expected_obs = ["run_id", "library_id", "cell_id"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    assert "prot" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    # assert list(input_mudata.var.columns) == expected_var, f"Input var columns should be: {expected_var}."
    assert all(
        key in list(input_mudata.mod["rna"].var.columns) for key in expected_var
    ), f"Input mod['rna'] var columns should be: {expected_var}, found: {input_mudata.mod["rna"].var.keys()}."
    assert all(
        key in list(input_mudata.mod["rna"].obs.columns) for key in expected_obs
    ), f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod["rna"].obs.keys()}."
    assert all(
        key in list(input_mudata.mod["prot"].var.columns) for key in expected_var
    ), f"Input mod['prot'] var columns should be: {expected_var}, found: {input_mudata.mod["prot"].var.keys()}."
    assert all(
        key in list(input_mudata.mod["prot"].obs.columns) for key in expected_obs
    ), f"Input mod ['prot'] obs columns should be: {expected_obs}, found: {input_mudata.mod["prot"].obs.keys()}."
    assert np.array_equal(
        input_mudata.mod["rna"].var["feature_type"].unique(), ["Gene Expression"]
    ), "Output X should only contain Gene Expression vars."
    assert np.array_equal(
        input_mudata.mod["prot"].var["feature_type"].unique(), ["Antibody Capture"]
    ), "Output X should only contain Gene Expression vars."


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    from importlib import resources
    shutil.copyfile(
        resources.files("openpipeline_testutils").joinpath("conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
