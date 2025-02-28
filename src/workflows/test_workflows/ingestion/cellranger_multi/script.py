from mudata import read_h5mu
from pathlib import Path
import shutil
import os
import sys
import pytest

##VIASH START
par = {"input": "input.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    for input_path in par["input"]:
        input_mudata = read_h5mu(input_path)

        assert list(input_mudata.mod.keys()) == ["rna", "prot", "vdj_t"]
        assert list(input_mudata.uns.keys()) == ["metrics_cellranger"]
        expected_metrics = [
            "Category",
            "Library Type",
            "Grouped By",
            "Group Name",
            "Metric Name",
            "Metric Value",
        ]
        assert (
            input_mudata.uns["metrics_cellranger"].columns.to_list() == expected_metrics
        )


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    from importlib import resources
    shutil.copyfile(
        resources.files("openpipeline_testutils").joinpath("conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
