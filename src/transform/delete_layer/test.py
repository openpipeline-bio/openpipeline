import sys
import pytest

from mudata import read_h5mu
from subprocess import CalledProcessError
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "./target/executable/transform/delete_layer/delete_layer",
    "executable": "target/executable/transform/delete_layer/delete_layer",
    "config": "src/transform/delete_layer/config.vsh.yaml",
    "resources_dir": "./resources_test/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.mark.parametrize("output_compression", [None, "gzip", "lzf"])
def test_delete_layer(run_component, random_h5mu_path, output_compression):
    temp_input = random_h5mu_path()
    temp_output = random_h5mu_path()

    # create input file
    input = read_h5mu(input_file)
    new_layer = input.mod["rna"].X
    input.mod["rna"].layers["test"] = new_layer
    assert "test" in input.mod["rna"].layers.keys()
    input.write_h5mu(temp_input)

    arguments = [
        "--input",
        temp_input,
        "--modality",
        "rna",
        "--layer",
        "test",
        "--output",
        temp_output,
    ]
    if output_compression:
        arguments.extend(["--output_compression", output_compression])

    # run command
    run_component(arguments)

    # check if output is correct
    assert temp_output.is_file()
    output = read_h5mu(temp_output)
    assert "test" not in output.mod["rna"].layers.keys()
    assert set(output.mod) == {"rna", "prot"}

    # Test that other data from input remains untouched
    del input.mod["rna"].layers["test"]
    assert_annotation_objects_equal(input, output)


def test_missing_layer_raises(run_component, random_h5mu_path):
    output = random_h5mu_path()
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--modality",
                "rna",
                "--layer",
                "test",
                "--output",
                str(output),
            ]
        )
    assert not output.is_file()
    assert "Layer 'test' is not present in modality rna." in err.value.stdout.decode(
        "utf-8"
    )


def test_missing_layer_missing_ok(run_component, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--layer",
            "test",
            "--output",
            str(output),
            "--missing_ok",
        ]
    )
    assert output.is_file()
    output_data = read_h5mu(output)
    assert "test" not in output_data.mod["rna"].layers.keys()

    # Test that other data from input remains untouched
    assert_annotation_objects_equal(input_file, output)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
