import os
import pytest
import shutil
import sys
from pathlib import Path
from mudata import read_h5mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal
from openpipelinetest_utils.utils import remove_annotation_column


##VIASH START
par = {"input": "input.h5mu", "og_input": "og_input.h5mu"}

meta = {
    "resources_dir": "resources_test/concat_test_data",
}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["og_input"])
    output_mudata = read_h5mu(par["input"])

    assert (
        input_mudata.n_mod == output_mudata.n_mod
    ), "Number of modalities should be the same"
    assert (
        input_mudata.mod.keys() == output_mudata.mod.keys()
    ), "Modalities should be the same"
    assert list(output_mudata.mod.keys()) == [
        "rna",
        "atac",
    ], "Modalities should be rna and atac"

    obs_cols_to_remove = []
    for top_n_vars in ("50", "100", "200", "500"):
        obs_cols_to_remove.append(f"pct_of_counts_in_top_{top_n_vars}_vars")

    obs_cols_to_remove.extend(["total_counts",
                               "num_nonzero_vars",
                               "fraction_mitochondrial",
                               "fraction_ribosomal",
                               "total_counts_mitochondrial",
                               "total_counts_ribosomal",
                               "pct_mitochondrial",
                               "pct_ribosomal"])
    var_cols_to_remove = ["obs_mean", "total_counts", "num_nonzero_obs", "pct_dropout", "mitochondrial", "ribosomal"]

    assert set(obs_cols_to_remove).issubset(
        set(output_mudata.mod["rna"].obs.columns.to_list())
    )
    assert set(var_cols_to_remove).issubset(
        set(output_mudata.mod["rna"].var.columns.to_list())
    )

    initial_mudata = remove_annotation_column(
        output_mudata, obs_cols_to_remove, axis="obs", modality_name="rna"
    )
    initial_mudata = remove_annotation_column(
        initial_mudata, var_cols_to_remove, axis="var", modality_name="rna"
    )

    assert_annotation_objects_equal(input_mudata, initial_mudata)


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    from importlib import resources
    shutil.copyfile(
        resources.files("openpipeline_testutils").joinpath("conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
