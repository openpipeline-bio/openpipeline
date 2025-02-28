import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/metadata/add_id/add_id",
    "resources_dir": "./resources_test/concat_test_data/",
    "cpus": 2,
    "config": "./src/metadata/add_id/config.vsh.yaml",
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


@pytest.mark.parametrize("output_compression", ["gzip", "lzf", None])
def test_add_id(
    run_component, small_mudata, small_mudata_path, random_h5mu_path, output_compression
):
    output_path = random_h5mu_path()

    args = [
        "--input",
        str(small_mudata_path),
        "--output",
        str(output_path),
        "--input_id",
        "test_id",
    ]

    small_mudata.obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.mod["mod1"].obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.mod["mod2"].obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.strings_to_categoricals()
    small_mudata.update()

    if output_compression:
        args.extend(["--output_compression", output_compression])
    # run component
    run_component(args)
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)

    assert_annotation_objects_equal(output_data, small_mudata)


def test_add_id_obs_output(
    run_component, small_mudata, small_mudata_path, random_h5mu_path
):
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            str(small_mudata_path),
            "--output",
            str(output_path),
            "--input_id",
            "test_id",
            "--obs_output",
            "test_key",
        ]
    )

    small_mudata.obs["test_key"] = ["test_id", "test_id"]
    small_mudata.mod["mod1"].obs["test_key"] = ["test_id", "test_id"]
    small_mudata.mod["mod2"].obs["test_key"] = ["test_id", "test_id"]
    small_mudata.strings_to_categoricals()
    small_mudata.update()

    assert output_path.is_file()
    output_data = read_h5mu(output_path)

    assert_annotation_objects_equal(output_data, small_mudata)


def test_add_id_observations_unique(
    run_component, small_mudata, small_mudata_path, random_h5mu_path
):
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            str(small_mudata_path),
            "--output",
            str(output_path),
            "--input_id",
            "test_id",
            "--make_observation_keys_unique",
        ]
    )

    small_mudata.obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.mod["mod1"].obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.mod["mod2"].obs["sample_id"] = ["test_id", "test_id"]
    small_mudata.strings_to_categoricals()
    small_mudata.obs.index = pd.Index(["test_id_obs1", "test_id_obs2"])
    small_mudata.mod["mod1"].obs.index = pd.Index(["test_id_obs1", "test_id_obs2"])
    small_mudata.mod["mod2"].obs.index = pd.Index(["test_id_obs1", "test_id_obs2"])
    small_mudata.update()

    assert output_path.is_file()
    output_data = read_h5mu(output_path)

    assert_annotation_objects_equal(output_data, small_mudata)


def test_add_id_overwrites_output_column(
    run_component, small_mudata, small_mudata_path, random_h5mu_path
):
    small_mudata.obs["already_exists"] = "alread_exists"
    for _, modality in small_mudata.mod.items():
        modality.obs["already_exists"] = "alread_exists"
    output_path = random_h5mu_path()

    # run component
    run_component(
        [
            "--input",
            str(small_mudata_path),
            "--output",
            str(output_path),
            "--input_id",
            "test_id",
            "--obs_output",
            "already_exists",
        ]
    )

    small_mudata.obs["already_exists"] = ["test_id", "test_id"]
    small_mudata.mod["mod1"].obs["already_exists"] = ["test_id", "test_id"]
    small_mudata.mod["mod2"].obs["already_exists"] = ["test_id", "test_id"]
    small_mudata.update()
    small_mudata.strings_to_categoricals()

    assert output_path.is_file()
    output_data = read_h5mu(output_path)

    assert_annotation_objects_equal(output_data, small_mudata)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
