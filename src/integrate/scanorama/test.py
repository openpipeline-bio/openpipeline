import sys
import pytest
from mudata import read_h5mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/docker/integrate/scanorama/scanorama",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "src/integrate/scanorama/config.vsh.yaml",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def input_with_batch(random_h5mu_path):
    tmp_input_path = random_h5mu_path()

    input_data = read_h5mu(input_file)
    mod = input_data.mod["rna"]
    number_of_obs = mod.n_obs
    mod.obs["batch"] = "A"
    column_index = mod.obs.columns.get_indexer(["batch"])
    mod.obs.iloc[slice(number_of_obs // 2, None), column_index] = "B"
    input_data.write(tmp_input_path)

    return tmp_input_path, input_data


def test_simple_integration(run_component, input_with_batch, random_h5mu_path):
    tmp_input_path, _ = input_with_batch
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            tmp_input_path,
            "--output",
            output_path,
            "--obs_batch",
            "batch",
            "--obsm_input",
            "X_pca",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()

    # check output
    data = read_h5mu(output_path)
    assert "X_scanorama" in data.mod["rna"].obsm

    # Delete output data in order to compare with input
    del data["rna"].obsm["X_scanorama"]
    assert_annotation_objects_equal(tmp_input_path, data)


def test_obsm_output(run_component, input_with_batch, random_h5mu_path):
    tmp_input_path, _ = input_with_batch
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            str(tmp_input_path),
            "--output",
            str(output_path),
            "--obsm_output",
            "X_test",
            "--obs_batch",
            "batch",
            "--obsm_input",
            "X_pca",
        ]
    )
    assert output_path.is_file()

    # check output
    data = read_h5mu(output_path)
    assert "X_test" in data.mod["rna"].obsm

    # Delete output data in order to compare with input
    del data["rna"].obsm["X_test"]
    assert_annotation_objects_equal(tmp_input_path, data)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
