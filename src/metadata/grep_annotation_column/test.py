import sys
import pytest
import pandas as pd
import numpy as np
import scipy
import math
from anndata import AnnData
from mudata import MuData, read_h5mu
from subprocess import CalledProcessError
from openpipelinetest_utils.asserters import assert_annotation_objects_equal
from openpipelinetest_utils.utils import remove_annotation_column
from itertools import permutations, islice


## VIASH START
meta = {
    "executable": "./target/executable/metadata/grep_annotation_column/grep_annotation_column",
    "resources_dir": "./resources_test/concat_test_data/",
    "cpus": 2,
    "config": "./src/metadata/grep_annotation_column/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def generate_h5mu():
    # generate data
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)
    tmp_mudata = MuData({"mod1": ad1, "mod2": ad2})
    return tmp_mudata


@pytest.fixture(
    params=[
        np.float32,
        np.float64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.int8,
        np.int16,
        np.int64,
    ]
)
def very_sparse_mudata(request):
    # NOTE: np.float16 is not a supported type in scipy!
    # See https://github.com/scipy/scipy/issues/20200#issuecomment-1982170609
    # and https://github.com/scipy/scipy/issues/20200
    rng = np.random.default_rng()
    shape = (10000, 200)
    random_counts = scipy.sparse.random(
        *shape,
        density=0.00001,
        format="csr",
        dtype=request.param,
        random_state=rng,
        data_rvs=lambda length: np.array([1] * length),
    )
    permutation_length_obs = int(math.log(shape[0]) / math.log(26)) + 1
    obs_index_perms = permutations(
        "abcdefghijklmnopqrstuvwxyz", r=permutation_length_obs
    )
    obs_index = pd.Index(["".join(x) for x in islice(obs_index_perms, shape[0])])
    obs = pd.DataFrame(index=obs_index)
    permutation_length_var = int(math.log(shape[1]) / math.log(26)) + 1
    var_index_perms = permutations(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ", r=permutation_length_var
    )
    var_index = pd.Index(["".join(x) for x in islice(var_index_perms, shape[1])])
    var = pd.DataFrame(index=var_index)
    mod = AnnData(X=random_counts, var=var, obs=obs)
    return MuData({"mod1": mod})


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column_default(
    run_component,
    generate_h5mu,
    random_h5mu_path,
    write_mudata_to_file,
    compression_format,
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--regex_pattern",
            "^var1",
            "--output_match_column",
            "test",
            "--output_compression",
            compression_format,
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [True, False, False]

    assert_annotation_objects_equal(
        input_path,
        remove_annotation_column(output_data, "test", "var", "mod1"),
        check_data=True,
    )


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column(
    run_component,
    generate_h5mu,
    random_h5mu_path,
    write_mudata_to_file,
    compression_format,
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--input_column",
            "Feat",
            "--regex_pattern",
            "^a",
            "--output_match_column",
            "test",
            "--output_compression",
            compression_format,
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [True, False, False]

    assert_annotation_objects_equal(
        input_path,
        remove_annotation_column(output_data, "test", "var", "mod1"),
        check_data=True,
    )


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column_fraction_column(
    run_component,
    generate_h5mu,
    random_h5mu_path,
    write_mudata_to_file,
    compression_format,
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--input_column",
            "Feat",
            "--regex_pattern",
            "^a",
            "--output_match_column",
            "test",
            "--output_fraction_column",
            "test_output_fraction",
            "--output_compression",
            compression_format,
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [True, False, False]
    assert output_data.mod["mod1"].obs["test_output_fraction"].to_list() == [
        1 / 6,
        4 / 15,
    ]
    output_object_without_obs_column = remove_annotation_column(
        output_data, ["test_output_fraction"], "obs", "mod1"
    )
    output_object_without_wo_output = remove_annotation_column(
        output_object_without_obs_column, ["test"], "var", "mod1"
    )
    assert_annotation_objects_equal(
        input_path, output_object_without_wo_output, check_data=True
    )


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_fraction_column_nothing_matches(
    run_component,
    generate_h5mu,
    random_h5mu_path,
    write_mudata_to_file,
    compression_format,
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--input_column",
            "Feat",
            "--regex_pattern",
            "^doesnotmatch$",
            "--output_match_column",
            "test",
            "--output_fraction_column",
            "test_output_fraction",
            "--output_compression",
            compression_format,
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [False, False, False]
    assert output_data.mod["mod1"].obs["test_output_fraction"].to_list() == [0.0, 0.0]
    output_object_without_obs_column = remove_annotation_column(
        output_data, ["test_output_fraction"], "obs", "mod1"
    )
    output_object_without_wo_output = remove_annotation_column(
        output_object_without_obs_column, ["test"], "var", "mod1"
    )
    assert_annotation_objects_equal(
        input_path, output_object_without_wo_output, check_data=True
    )


def test_fraction_column_with_no_counts(
    run_component, generate_h5mu, random_h5mu_path, write_mudata_to_file
):
    output_path = random_h5mu_path()
    generate_h5mu.mod["mod1"].X[0] = 0
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--input_column",
            "Feat",
            "--regex_pattern",
            "^a",
            "--output_match_column",
            "test",
            "--output_fraction_column",
            "test_output_fraction",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [True, False, False]
    assert output_data.mod["mod1"].obs["test_output_fraction"].to_list() == [
        0.0,
        4 / 15,
    ]
    output_object_without_obs_column = remove_annotation_column(
        output_data, ["test_output_fraction"], "obs", "mod1"
    )
    output_object_without_wo_output = remove_annotation_column(
        output_object_without_obs_column, ["test"], "var", "mod1"
    )
    assert_annotation_objects_equal(
        input_path, output_object_without_wo_output, check_data=True
    )


def test_fraction_column_very_sparse(
    run_component, very_sparse_mudata, random_h5mu_path, write_mudata_to_file
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(very_sparse_mudata)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--regex_pattern",
            "^a",
            "--output_match_column",
            "test",
            "--output_fraction_column",
            "test_output_fraction",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    output_object_without_obs_column = remove_annotation_column(
        output_data, ["test_output_fraction"], "obs", "mod1"
    )
    output_object_without_wo_output = remove_annotation_column(
        output_object_without_obs_column, ["test"], "var", "mod1"
    )
    assert_annotation_objects_equal(
        input_path, output_object_without_wo_output, check_data=True
    )


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_fraction_column_input_layer(
    run_component,
    generate_h5mu,
    random_h5mu_path,
    write_mudata_to_file,
    compression_format,
):
    output_path = random_h5mu_path()
    generate_h5mu.mod["mod1"].layers["test_data"] = generate_h5mu.mod["mod1"].X.copy()
    generate_h5mu.mod["mod1"].X = None
    input_path = write_mudata_to_file(generate_h5mu)

    # run component
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "mod1",
            "--matrix",
            "var",
            "--input_column",
            "Feat",
            "--regex_pattern",
            "^a",
            "--output_match_column",
            "test",
            "--output_fraction_column",
            "test_output_fraction",
            "--input_layer",
            "test_data",
            "--output_compression",
            compression_format,
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod["mod1"].var.columns.to_list()
    assert output_data.mod["mod1"].var["test"].to_list() == [True, False, False]
    output_data.mod["mod1"].obs["test_output_fraction"].to_list() == [1 / 6, 4 / 15]
    output_object_without_obs_column = remove_annotation_column(
        output_data, ["test_output_fraction"], "obs", "mod1"
    )
    output_object_without_wo_output = remove_annotation_column(
        output_object_without_obs_column, ["test"], "var", "mod1"
    )
    assert_annotation_objects_equal(
        input_path, output_object_without_wo_output, check_data=True
    )


def test_missing_column(
    run_component, generate_h5mu, random_h5mu_path, write_mudata_to_file
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--modality",
                "mod1",
                "--matrix",
                "var",
                "--input_column",
                "filliberke",
                "--regex_pattern",
                "^a",
                "--output_match_column",
                "test",
                "--output_compression",
                "gzip",
            ]
        )
    assert (
        "ValueError: Column filliberke could not be found for modality mod1"
        in err.value.stdout.decode("utf-8")
    )


def test_invalid_regex_pattern(
    run_component, generate_h5mu, random_h5mu_path, write_mudata_to_file
):
    output_path = random_h5mu_path()
    input_path = write_mudata_to_file(generate_h5mu)

    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--modality",
                "mod1",
                "--matrix",
                "var",
                "--input_column",
                "Feat",
                "--regex_pattern",
                "(a",
                "--output_match_column",
                "test",
                "--output_compression",
                "gzip",
            ]
        )
    assert (
        "ValueError: (a is not a valid regular expression pattern."
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-s"]))
