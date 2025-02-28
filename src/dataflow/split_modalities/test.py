import sys
import pytest
import mudata as md
import anndata as ad
import pandas as pd
import re
from openpipelinetest_utils.asserters import assert_annotation_objects_equal
from textwrap import dedent

## VIASH START
meta = {
    "name": "./target/native/dataflow/split_modalities/split_modalities",
    "resources_dir": "./resources_test/",
    "config": "./src/dataflow/split_modalities/config.vsh.yaml",
    "executable": "./target/docker/dataflow/split_modalities/split_modalities",
}
## VIASH END


@pytest.fixture
def input_modality_1():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    return ad1


@pytest.fixture
def input_modality_2():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    return ad2


@pytest.fixture
def input_h5mu(input_modality_1, input_modality_2):
    tmp_mudata = md.MuData({"mod1": input_modality_1, "mod2": input_modality_2})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


@pytest.mark.parametrize("compression", ["gzip", None])
def test_split(
    run_component,
    random_path,
    input_h5mu,
    input_h5mu_path,
    input_modality_1,
    input_modality_2,
    compression,
):
    output_dir = random_path()
    output_types = random_path(extension="csv")
    args = [
        "--input",
        input_h5mu_path,
        "--output",
        str(output_dir),
        "--output_types",
        str(output_types),
    ]
    if compression:
        args += ["--output_compression", compression]
    run_component(args)
    assert output_types.is_file()
    assert output_dir.is_dir()

    # check output dir
    dir_content = [
        h5mu_file
        for h5mu_file in output_dir.iterdir()
        if h5mu_file.suffix == ".h5mu" and h5mu_file != input_h5mu_path
    ]
    mod1_file = output_dir / f"{input_h5mu_path.stem}_mod1.h5mu"
    mod2_file = output_dir / f"{input_h5mu_path.stem}_mod2.h5mu"
    assert set(dir_content) == set([mod1_file, mod2_file])
    mod1 = md.read_h5mu(mod1_file)
    mod2 = md.read_h5mu(mod2_file)
    assert mod1.n_mod == 1
    assert mod2.n_mod == 1

    assert_annotation_objects_equal(mod1.mod["mod1"], input_modality_1)
    assert_annotation_objects_equal(mod2.mod["mod2"], input_modality_2)

    assert mod1.n_obs == input_h5mu.n_obs
    assert mod2.n_obs == input_h5mu.n_obs

    # When a var_key is only present for one modality, it is prefixed by the name of the
    # modality followed by a colon and the name of the key (in the global .var).
    replace_regex = r"(^mod1:|^mod2:)"
    expected_var_keys = {
        re.sub(replace_regex, "", col_name) for col_name in input_h5mu.var_keys()
    }
    assert set(mod1.var_keys()) | set(mod2.var_keys()) == expected_var_keys

    assert set(mod1.var_keys()) == set(input_h5mu.mod["mod1"].var.columns)
    assert set(mod2.var_keys()) == set(input_h5mu.mod["mod2"].var.columns)

    expected_csv_output = dedent(
        f"""\
        name,filename
        mod1,{mod1_file.name}
        mod2,{mod2_file.name}
        """
    )
    with open(output_types, "r") as open_csv_file:
        result = open_csv_file.read()
        assert result == expected_csv_output


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
