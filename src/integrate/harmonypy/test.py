import sys
import pytest
import mudata
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/integrate/harmonypy/harmonypy",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "src/integrate/harmony/config.vsh.yaml",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


def test_harmonypy(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--obsm_input",
            "X_pca",
            "--obsm_output",
            "X_pca_int",
            "--obs_covariates",
            "harmony_integration_leiden_1.0",
            "--output",
            output_path,
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = mudata.read_h5mu(output_path)
    assert "X_pca_int" in output_data.mod["rna"].obsm
    assert (
        output_data.mod["rna"].obsm["X_pca_int"].shape
        == output_data.mod["rna"].obsm["X_pca"].shape
    )
    # Delete output slots in order to tests for equality with
    # input data
    del output_data["rna"].obsm["X_pca_int"]
    assert_annotation_objects_equal(input_file, output_data)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
