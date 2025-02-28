import mudata as mu
import sys
from pathlib import Path
import pytest
import numpy as np
from subprocess import CalledProcessError
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/filter/delimit_fraction/delimit_fraction",
    "resources_dir": "resources_test/",
    "config": "./src/filter/delimit_fraction/config.vsh.yaml",
}

## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def original_input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_h5mu(original_input_path):
    input_data = mu.read_h5mu(original_input_path)
    input_data.mod["rna"].obs["test_fraction"] = np.random.rand(
        input_data.mod["rna"].n_obs
    )
    return input_data


@pytest.fixture
def input_h5mu_string_data(original_input_path):
    input_data = mu.read_h5mu(original_input_path)
    string_data = ["these", "are", "random", "values"]
    input_data.mod["rna"].obs["test_fraction"] = np.random.choice(
        string_data, input_data.mod["rna"].n_obs
    )
    return input_data


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    output_path = random_h5mu_path()
    input_h5mu.write(output_path)
    return output_path


@pytest.fixture
def input_path_string_data(input_h5mu_string_data, random_h5mu_path):
    output_path = random_h5mu_path()
    input_h5mu_string_data.write(output_path)
    return output_path


def test_filter_nothing(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--min_fraction",
            "0",
            "--max_fraction",
            "1",
            "--output_compression",
            "gzip",
            "--obs_name_filter",
            "test_output",
            "--obs_fraction_column",
            "test_fraction",
        ]
    )
    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    assert "test_output" in mu_out.mod["rna"].obs
    assert mu_out.mod["rna"].obs["test_output"].all()

    mu_out.mod["rna"].obs = mu_out.mod["rna"].obs.drop(["test_output"], axis=1)
    mu_out.obs = mu_out.obs.drop(["rna:test_output"], axis=1)
    mu_out.update()
    assert_annotation_objects_equal(input_path, mu_out)


def test_filtering_a_little(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--min_fraction",
            "0.5",
            "--max_fraction",
            "0.7",
            "--output_compression",
            "gzip",
            "--obs_name_filter",
            "test_output",
            "--obs_fraction_column",
            "test_fraction",
        ]
    )

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    assert not mu_out.mod["rna"].obs["test_output"].all()
    assert mu_out.mod["rna"].obs["test_output"].any()

    mu_out.mod["rna"].obs = mu_out.mod["rna"].obs.drop(["test_output"], axis=1)
    mu_out.obs = mu_out.obs.drop(["rna:test_output"], axis=1)
    mu_out.update()
    assert_annotation_objects_equal(input_path, mu_out)


def test_filtering_wrong_data_raises(
    run_component, input_path_string_data, random_h5mu_path
):
    output_path = random_h5mu_path()
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path_string_data,
                "--output",
                output_path,
                "--min_fraction",
                "0.5",
                "--max_fraction",
                "0.7",
                "--output_compression",
                "gzip",
                "--obs_name_filter",
                "test_output",
                "--obs_fraction_column",
                "test_fraction",
            ]
        )
    assert (
        "Column 'test_fraction' does not contain float datatype."
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    exit(pytest.main([__file__]))
