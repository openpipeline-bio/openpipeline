from mudata import read_h5mu
import numpy as np
import shutil
import os
import sys
from pathlib import Path
import pytest

##VIASH START
par = {
    "input": "bbknn_knn/output.h5mu"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_obsm = ["X_leiden_bbknn_umap"]
    expected_obs = ["cell_type_pred", "cell_type_probability"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm), f"Input mod['rna'] obs columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    assert all(key in list(input_mudata.mod["rna"].obs) for key in expected_obs), f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    shutil.copyfile(os.path.join(meta['resources_dir'], "openpipelinetestutils", "conftest.py"),
                    os.path.join(HERE_DIR, "conftest.py"))
    sys.exit(pytest.main(["--import-mode=importlib"]))