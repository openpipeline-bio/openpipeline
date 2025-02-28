import sys
import pytest
import numpy as np
from mudata import read_h5mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "target/docker/scale/scale",
    "resources_dir": "resources_test",
    "config": "src/transform/scale/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def input_data(input_path):
    return read_h5mu(input_path)


def test_scaling_input_layer(
    run_component, input_data, write_mudata_to_file, random_h5mu_path
):
    """
    The component must select the correct input layer.
    """
    input_data.mod["rna"].layers["test_layer"] = input_data.mod["rna"].X.copy()
    del input_data.mod["rna"].X
    input_path = write_mudata_to_file(input_data)
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_file,
            "--input_layer",
            "test_layer",
            "--ouput_compression",
            "gzip",
        ]
    )

    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod["rna"].X
    mean = np.mean(output_x, axis=0, dtype=np.float64)
    variance = np.multiply(output_x, output_x).mean(axis=0, dtype=np.float64) - mean**2
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))
    # Remove the output layer in order to
    # test equality for the other data slots
    del output_data.mod["rna"].X
    assert_annotation_objects_equal(input_path, output_data)


def test_scaling_output_layer(run_component, random_h5mu_path, input_path):
    """
    Output data must create the specified output layer.
    """
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_file,
            "--output_layer",
            "scaled",
            "--ouput_compression",
            "gzip",
        ]
    )

    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    assert "scaled" in output_data.mod["rna"].layers
    output_scaled = output_data.mod["rna"].layers["scaled"]
    mean = np.mean(output_scaled, axis=0, dtype=np.float64)
    variance = (
        np.multiply(output_scaled, output_scaled).mean(axis=0, dtype=np.float64)
        - mean**2
    )
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))
    # Remove the output layer in order to
    # test equality for the other data slots
    del output_data.mod["rna"].layers["scaled"]
    assert_annotation_objects_equal(input_path, output_data)


def test_scaling(run_component, random_h5mu_path, input_path):
    """
    Output data must be centered around mean 0 and it has unit variance.
    """
    output_file = random_h5mu_path()

    run_component(
        ["--input", input_path, "--output", output_file, "--ouput_compression", "gzip"]
    )

    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod["rna"].X
    mean = np.mean(output_x, axis=0, dtype=np.float64)
    variance = np.multiply(output_x, output_x).mean(axis=0, dtype=np.float64) - mean**2
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))
    # Restore the output layer in order to
    # test equality for the other data slots
    input_x = read_h5mu(input_path)["rna"].X
    output_data["rna"].X = input_x
    assert_annotation_objects_equal(input_path, output_data)


def test_scaling_noncenter(run_component, random_h5mu_path, input_path):
    """
    Check if centering can be disabled.
    """
    output_file = random_h5mu_path()

    run_component(
        ["--input", input_path, "--output", str(output_file), "--zero_center", "false"]
    )
    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod["rna"].X
    mean = np.mean(output_x, axis=0, dtype=np.float64)
    assert not np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    # Restore the output layer in order to
    # test equality for the other data slots
    input_x = read_h5mu(input_path)["rna"].X
    output_data["rna"].X = input_x
    assert_annotation_objects_equal(input_path, output_data)


def test_scaling_maxvalue(run_component, random_h5mu_path, input_path):
    """
    Check if output data is clipped when using --max_value
    """
    output_file = random_h5mu_path()

    run_component(
        ["--input", input_path, "--output", str(output_file), "--max_value", "0.5"]
    )
    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod["rna"].X
    assert np.all(output_x <= 0.5)
    # Restore the output layer in order to
    # test equality for the other data slots
    input_x = read_h5mu(input_path)["rna"].X
    output_data["rna"].X = input_x
    assert_annotation_objects_equal(input_path, output_data)


def test_scaling_modality(run_component, random_h5mu_path, input_path):
    """
    Check if 'rna' modality remain untouched when using '--modality prot' argument.
    """
    output_file = random_h5mu_path()

    run_component(
        ["--input", input_path, "--output", str(output_file), "--modality", "prot"]
    )
    assert output_file.is_file()
    output_data = read_h5mu(output_file)

    output_prot = output_data.mod["prot"].X
    mean = np.mean(output_prot, axis=0, dtype=np.float64)
    variance = (
        np.multiply(output_prot, output_prot).mean(axis=0, dtype=np.float64) - mean**2
    )
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))
    # Restore the output layer in order to
    # test equality for the other data slots
    input_x = read_h5mu(input_path)["prot"].X
    output_data["prot"].X = input_x
    assert_annotation_objects_equal(input_path, output_data)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
