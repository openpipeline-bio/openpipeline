import sys
import pytest
import subprocess
from mudata import MuData, read_h5mu
import pandas as pd
import numpy as np
import re
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/dataflow/merge/merge",
    "resources_dir": "./resources_test/merge_test_data/",
    "config": "./src/dataflow/merge/config.vsh.yml",
}
## VIASH END


@pytest.fixture
def mudata_non_overlapping_observations(request, random_h5mu_path):
    mudata_to_change_path = request.getfixturevalue(request.param)
    # mudata_to_change_path = small_mudata_mod1_path
    temp_h5mu = random_h5mu_path()
    mudata_to_change = read_h5mu(mudata_to_change_path)
    # Remove 1 observation
    removed_observation_name = mudata_to_change.obs.index[-1]
    mudata_to_change = mudata_to_change[: mudata_to_change.n_obs - 1]
    mudata_to_change.write(temp_h5mu, compression="gzip")
    return temp_h5mu, removed_observation_name


@pytest.fixture
def extra_var_column_value():
    return "bar"


@pytest.fixture
def extra_var_column_name():
    return "test"


@pytest.fixture
def mudata_with_extra_var_column(
    random_h5mu_path, request, extra_var_column_value, extra_var_column_name
):
    [sample1_path, sample2_path] = request.getfixturevalue(request.param)
    result = []
    for sample_path, column_value in (
        (sample1_path, extra_var_column_value),
        (sample2_path, np.nan),
    ):
        sample = read_h5mu(sample_path)
        mod_names = list(sample.mod.keys())
        assert len(mod_names) == 1
        mod_name = mod_names[0]
        sample.mod[mod_name].var[extra_var_column_name] = column_value
        sample.var = sample.var.convert_dtypes(
            infer_objects=True,
            convert_integer=True,
            convert_string=False,
            convert_boolean=True,
            convert_floating=False,
        )
        new_path = random_h5mu_path()
        sample.write(new_path)
        result.append(new_path)
    return result


def test_merge(
    run_component, random_h5mu_path, small_mudata_mod1_path, small_mudata_mod2_path
):
    """
    Test a simple merge with fully overlapping observations
    """
    output_path = random_h5mu_path()
    # input_sample1_path, input_sample2_path = split_small_mudata_path
    args = [
        "--input",
        small_mudata_mod1_path,
        "--input",
        small_mudata_mod2_path,
        "--output",
        output_path,
        "--output_compression",
        "gzip",
    ]
    run_component(args)

    assert output_path.is_file()
    concatenated_data = read_h5mu(output_path)
    data_sample1 = read_h5mu(small_mudata_mod1_path)
    data_sample2 = read_h5mu(small_mudata_mod2_path)

    expected_concatenated_data = MuData(
        {"mod1": data_sample1.mod["mod1"], "mod2": data_sample2.mod["mod2"]}
    )
    assert_annotation_objects_equal(concatenated_data, expected_concatenated_data)


@pytest.mark.parametrize(
    "mudata_non_overlapping_observations",
    ["small_mudata_mod1_path"],
    indirect=["mudata_non_overlapping_observations"],
)
def test_merge_non_overlapping_observations(
    run_component,
    mudata_non_overlapping_observations,
    random_h5mu_path,
    small_mudata_mod2_path,
):
    """
    Merge with differing observations in the samples
    """
    edited_h5mu_path, removed_observation_name = mudata_non_overlapping_observations
    output_path = random_h5mu_path()
    # Remove 1 observation
    run_component(
        [
            "--input",
            edited_h5mu_path,
            "--input",
            small_mudata_mod2_path,
            "--output",
            output_path,
        ]
    )

    assert output_path.is_file()
    concatenated_data = read_h5mu(output_path, backed=False)
    data_sample1 = read_h5mu(edited_h5mu_path, backed=False)
    data_sample2 = read_h5mu(small_mudata_mod2_path, backed=False)

    expected_concatenated_data = MuData(
        {"mod1": data_sample1.mod["mod1"], "mod2": data_sample2.mod["mod2"]}
    )

    assert set(concatenated_data.obs_names) == (
        set(data_sample1.obs_names) | set(data_sample2.obs_names)
    )
    assert concatenated_data[removed_observation_name:]["mod1"].n_obs == 0
    assert concatenated_data[removed_observation_name:]["mod2"].n_obs == 1

    np.testing.assert_equal(
        concatenated_data.copy()[removed_observation_name:]["mod2"].X.data,
        data_sample2.copy()[removed_observation_name:]["mod2"].X.data,
    )

    assert_annotation_objects_equal(concatenated_data, expected_concatenated_data)


@pytest.mark.parametrize(
    "extra_var_column_name,extra_var_column_value,expected",
    [
        ("test", "bar", "bar"),
        ("test", True, True),
        ("test", 0.1, 0.1),
        ("test", np.nan, pd.NA),
    ],
)
@pytest.mark.parametrize(
    "mudata_with_extra_var_column",
    ["split_small_mudata_path"],
    indirect=["mudata_with_extra_var_column"],
)
def test_boolean_and_na_types(
    run_component,
    mudata_with_extra_var_column,
    extra_var_column_name,
    expected,
    random_h5mu_path,
):
    """
    Test if merging booleans of NAs results in the .var .obs column being writeable
    """
    input_sample1_path, input_sample2_path = mudata_with_extra_var_column
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_sample1_path,
            "--input",
            input_sample2_path,
            "--output",
            output_path,
        ]
    )
    assert output_path.is_file()
    merged_data = read_h5mu(output_path, backed=False)
    first_sample_mod = list(read_h5mu(input_sample1_path).mod)[0]
    second_sample_mod = list(read_h5mu(input_sample2_path).mod)[0]

    expected_merged_data = MuData(
        {
            "mod1": read_h5mu(input_sample1_path).mod["mod1"],
            "mod2": read_h5mu(input_sample2_path).mod["mod2"],
        }
    )

    if not pd.isna(expected):
        assert merged_data.var.loc["var1"][extra_var_column_name] == expected
        assert (
            merged_data.mod[first_sample_mod].var.loc["var1"][extra_var_column_name]
            == expected
        )
    else:
        assert pd.isna(merged_data.var.loc["var1"][extra_var_column_name])
        assert pd.isna(
            merged_data.mod[first_sample_mod].var.loc["var1"][extra_var_column_name]
        )
    assert pd.isna(merged_data.var.loc["var4"][extra_var_column_name])
    assert pd.isna(
        merged_data.mod[second_sample_mod].var.loc["var4"][extra_var_column_name]
    )

    assert_annotation_objects_equal(merged_data, expected_merged_data)


def test_same_modalities_raises(
    run_component, random_h5mu_path, split_small_mudata_path
):
    """
    Raise when trying to merge modalities with the same name.
    """
    input_sample1_path, input_sample2_path = split_small_mudata_path
    input_sample2_edited_path = random_h5mu_path()
    output_path = random_h5mu_path()
    data_sample2 = read_h5mu(input_sample2_path)
    data_sample2 = MuData({"mod1": data_sample2.mod["mod2"]})
    data_sample2.write(input_sample2_edited_path, compression="gzip")

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_sample1_path,
                "--input",
                input_sample2_edited_path,
                "--output",
                output_path,
            ]
        )
    assert re.search(
        r"ValueError: Modality 'mod1' was found in more than 1 sample\.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
