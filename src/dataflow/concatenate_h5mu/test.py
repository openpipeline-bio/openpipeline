import mudata as md
import anndata as ad
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import pytest
import re
import sys
from openpipelinetest_utils.utils import remove_annotation_column
from operator import attrgetter

## VIASH START
meta = {
    "executable": "./target/docker/dataflow/concatenate_h5mu/concatenate_h5mu",
    "resources_dir": "./resources_test/concat_test_data/",
    "cpus": 2,
    "config": "./src/dataflow/concatenate_h5mu/config.vsh.yaml",
}
## VIASH END

meta["cpus"] = 1 if not meta["cpus"] else meta["cpus"]


@pytest.fixture
def sample_1_modality_1():
    """
    >>> ad1.obs
          Obs1 Shared_obs
     obs1    A          B
     obs2    C          D

    >>> ad1.var
                         Feat1  Shared_feat
    var1                     a            b
    var2                     c            d
    overlapping_var_mod1     e            f

    >>> ad1.X
    array([[1, 2, 3],
           [4, 5, 6]])
    """

    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]],
        index=["obs1", "obs2"],
        columns=["var1", "var2", "overlapping_var_mod1"],
    )
    obs = pd.DataFrame(
        [["A", "B"], ["C", "D"]], index=df.index, columns=["Obs1", "Shared_obs"]
    )
    var = pd.DataFrame(
        [["a", "b"], ["c", "d"], ["e", "f"]],
        index=df.columns,
        columns=["Feat1", "Shared_feat"],
    )
    varm = np.random.rand(df.columns.size, 5)
    ad1 = ad.AnnData(
        df,
        obs=obs,
        var=var,
        varm={"random_vals_mod1": varm},
        uns={
            "uns_unique_to_sample1": pd.DataFrame(
                ["foo"], index=["bar"], columns=["col1"]
            ),
            "overlapping_uns_key": pd.DataFrame(
                ["jing"], index=["jang"], columns=["col2"]
            ),
        },
    )
    # ad1 = ad.AnnData(df, obs=obs, var=var)
    return ad1


@pytest.fixture
def sample_1_input_modality_2():
    """
    >>> ad2.X
    array([[ 7,  8],
           [ 9, 10],
           [11, 12]])

    >>> ad2.obs
         Obs2 Obs3 Shared_obs
    obs3    E    F          G
    obs4    H    I          J
    obs5    K    L          M

    >>> ad2.var
         Feat2  Shared_feat
    var4     d            e
    var5     f            g

    """
    df = pd.DataFrame(
        [[7, 8], [9, 10], [11, 12]],
        index=["obs3", "obs4", "obs5"],
        columns=["var3", "var4"],
    )
    obs = pd.DataFrame(
        [["E", "F", "G"], ["H", "I", "J"], ["K", "L", "M"]],
        index=df.index,
        columns=["Obs2", "Obs3", "Shared_obs"],
    )
    var = pd.DataFrame(
        [["d", "e"], ["f", "g"]], index=df.columns, columns=["Feat2", "Shared_feat"]
    )
    ad2 = ad.AnnData(df, obs=obs, var=var)
    return ad2


@pytest.fixture
def sample_1_h5mu(sample_1_modality_1, sample_1_input_modality_2):
    tmp_mudata = md.MuData(
        {"mod1": sample_1_modality_1, "mod2": sample_1_input_modality_2}
    )
    return tmp_mudata


@pytest.fixture
def sample_2_modality_1():
    """
    >>> ad3.X
    array([[13, 14],
           [15, 16],
           [17, 18]])

    >>> ad3.var
                         Feat3 Shared_feat
    var5                     h           i
    overlapping_var_mod1     j           k

    >>> ad3.obs
         Obs4 Obs5 Shared_obs
    obs6    O    P          Q
    obs7    R    S          T
    obs8    U    V          W
    """
    df = pd.DataFrame(
        [[13, 14], [15, 16], [17, 18]],
        index=["obs6", "obs7", "obs8"],
        columns=["var5", "overlapping_var_mod1"],
    )
    obs = pd.DataFrame(
        [["O", "P", "Q"], ["R", "S", "T"], ["U", "V", "W"]],
        index=df.index,
        columns=["Obs4", "Obs5", "Shared_obs"],
    )
    var = pd.DataFrame(
        [["h", "i"], ["j", "k"]], index=df.columns, columns=["Feat3", "Shared_feat"]
    )
    ad3 = ad.AnnData(
        df,
        obs=obs,
        var=var,
        uns={
            "uns_unique_to_sample2": pd.DataFrame(
                ["baz"], index=["qux"], columns=["col3"]
            ),
            "overlapping_uns_key": pd.DataFrame(
                ["ping"], index=["pong"], columns=["col4"]
            ),
        },
    )
    return ad3


@pytest.fixture
def sample_2_modality_2():
    """
    >>> ad4.X
    array([[19, 20, 21],
           [22, 23, 24]])

    >>> ad4.obs
         Obs6 Shared_obs
    obs8    X          Y
    obs9    Z         AA

    >>> ad4.var
         Feat4 Shared_feat
    var6     l           m
    var7     n           o
    var8     p           q
    """
    df = pd.DataFrame(
        [[19, 20, 21], [22, 23, 24]],
        index=["obs8", "obs9"],
        columns=["var6", "var7", "var8"],
    )
    obs = pd.DataFrame(
        [["X", "Y"], ["Z", "AA"]], index=df.index, columns=["Obs6", "Shared_obs"]
    )
    var = pd.DataFrame(
        [["l", "m"], ["n", "o"], ["p", "q"]],
        index=df.columns,
        columns=["Feat4", "Shared_feat"],
    )
    varm = np.random.rand(df.columns.size, 3)
    ad4 = ad.AnnData(df, obs=obs, var=var, varm={"random_vals_mod2": varm})
    # ad4 = ad.AnnData(df, obs=obs, var=var)
    return ad4


@pytest.fixture
def sample_2_h5mu(sample_2_modality_1, sample_2_modality_2):
    tmp_mudata = md.MuData({"mod1": sample_2_modality_1, "mod2": sample_2_modality_2})
    return tmp_mudata


@pytest.fixture
def sample_3_modality_1():
    """
    >>> ad3.X
    array([[25]])

    >>> ad3.obs
          Obs7
    obs10   AB

    >>> ad3.var
         Feat4
    var9     r

    """
    df = pd.DataFrame([[25]], index=["obs10"], columns=["var9"])
    obs = pd.DataFrame([["AB"]], index=df.index, columns=["Obs7"])
    var = pd.DataFrame([["r"]], index=df.columns, columns=["Feat4"])
    ad5 = ad.AnnData(df, obs=obs, var=var)
    return ad5


@pytest.fixture
def sample_3_modality_3():
    """
    >>> ad6.X
    array([[ 26,  32,  33, 453],
        [ 34,  35,  36, 543]])

    >>> ad6.var
          Feat5 Feat6 Feat7 Feat8
    var10     s     t     u     v
    var11     w     x     y     z
    var12    aa    ab    ac    ad
    var13    ae    af    ag    ah

    >>> ad6.obs
          Obs8 Obs9 obs10 obs11
    obs11   AC   AD    AE    AF
    obs12   AG   AH    AI    AJ
    """
    df = pd.DataFrame(
        [[26, 32, 33, 453], [34, 35, 36, 543]],
        index=["obs11", "obs12"],
        columns=["var10", "var11", "var12", "var13"],
    )
    obs = pd.DataFrame(
        [["AC", "AD", "AE", "AF"], ["AG", "AH", "AI", "AJ"]],
        index=df.index,
        columns=["Obs8", "Obs9", "obs10", "obs11"],
    )
    var = pd.DataFrame(
        [
            ["s", "t", "u", "v"],
            ["w", "x", "y", "z"],
            ["aa", "ab", "ac", "ad"],
            ["ae", "af", "ag", "ah"],
        ],
        index=df.columns,
        columns=["Feat5", "Feat6", "Feat7", "Feat8"],
    )
    ad6 = ad.AnnData(df, obs=obs, var=var)
    return ad6


@pytest.fixture
def sample_3_h5mu(sample_3_modality_1, sample_3_modality_3):
    tmp_mudata = md.MuData({"mod1": sample_3_modality_1, "mod3": sample_3_modality_3})
    return tmp_mudata


@pytest.fixture
def wrap_anndata_to_mudata():
    def wrapper(anndata_obj, mod_name="mod"):
        return md.MuData({mod_name: anndata_obj})

    return wrapper


@pytest.fixture
def change_column_contents():
    def wrapper(mudata_obj, annotation_frame_name, column_name, values_per_modality):
        mudata_obj.update()
        get_frame = attrgetter(annotation_frame_name)
        modality_columns = []
        for mod_name, col_value in values_per_modality.items():
            modality = mudata_obj.mod[mod_name]
            annotation_frame = get_frame(modality)
            annotation_frame[column_name] = col_value
            modality_columns.append(annotation_frame[column_name])
        mudata_obj.update()
        global_annotation_frame = get_frame(mudata_obj)
        if column_name in global_annotation_frame.columns:
            updated_global_column = pd.concat(modality_columns, copy=True, join="inner")
            no_duplicates = updated_global_column.reset_index().drop_duplicates(
                subset=["index"]
            )
            no_duplicates = no_duplicates.set_index("index")
            global_annotation_frame[column_name] = no_duplicates
        setattr(
            mudata_obj,
            annotation_frame_name,
            global_annotation_frame.convert_dtypes(
                infer_objects=True,
                convert_integer=True,
                convert_string=False,
                convert_boolean=True,
                convert_floating=False,
            ),
        )

    return wrapper


def test_concatenate_samples_with_same_observation_ids_raises(
    run_component,
    wrap_anndata_to_mudata,
    write_mudata_to_file,
    sample_1_modality_1,
    sample_2_modality_1,
    random_h5mu_path,
):
    """
    Test how concat handles overlapping observation IDs.
    This should raise.
    """
    # introduce an overlapping observation
    input_1_mudata = wrap_anndata_to_mudata(sample_1_modality_1)
    old_obs_names = sample_2_modality_1.obs_names
    new_obs_names = old_obs_names.where(
        old_obs_names.isin([old_obs_names[0]]), sample_1_modality_1.obs.index[0]
    )
    sample_2_modality_1.obs_names = new_obs_names
    input_2_mudata = wrap_anndata_to_mudata(sample_2_modality_1)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_id",
                "foo;bar",
                "--input",
                write_mudata_to_file(input_1_mudata),
                "--input",
                write_mudata_to_file(input_2_mudata),
                "--output",
                random_h5mu_path(),
                "--other_axis_mode",
                "move",
                "--output_compression",
                "gzip",
            ]
        )
    assert (
        "ValueError: Observations are not unique across samples."
        in err.value.stdout.decode("utf-8")
    )


def test_concat_different_var_columns_per_sample(
    run_component, sample_1_h5mu, sample_2_h5mu, random_h5mu_path, write_mudata_to_file
):
    """
    Test what happens when concatenating samples with differing auxiliary
    (like in .var) columns (present in 1 sample, absent in other).
    When concatenating the samples, all columns should be present in the
    resulting object, filling the values from samples with the missing
    column with NA.

    Looking at Shared_feat here:

                   mod1      mod2
    sample 1    present   present
    sample 2          x         x
    """
    output_path = random_h5mu_path()
    # Before removing the 'Shared_feat' column from one of the samples,
    # check if they are present in both
    assert "Shared_feat" in sample_1_h5mu.var_keys()
    assert "Shared_feat" in sample_2_h5mu.var_keys()

    sample_2_h5mu = remove_annotation_column(sample_2_h5mu, ["Shared_feat"], axis="var")
    assert "Shared_feat" in sample_1_h5mu.var_keys()
    assert "Shared_feat" not in sample_2_h5mu.var_keys()

    # 'Shared_feat' column is not missing from sample2, which is what this test is about
    input_sample1_path = write_mudata_to_file(sample_1_h5mu)
    input_sample2_path = write_mudata_to_file(sample_2_h5mu)

    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            input_sample1_path,
            "--input",
            input_sample2_path,
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )

    assert Path(output_path).is_file()
    concatenated_data = md.read(output_path)

    data_sample1 = md.read(input_sample1_path)
    data_sample2 = md.read(input_sample2_path)

    assert (
        concatenated_data.n_vars
        == data_sample1.var.index.union(data_sample2.var.index).size
    )

    for mod_name in ("mod1", "mod2"):
        # Check if all features are present
        concatenated_mod = concatenated_data.mod[mod_name]
        sample1_original_mod = data_sample1.mod[mod_name]
        sample2_original_mod = data_sample2.mod[mod_name]

        original_var_keys = set(
            sample1_original_mod.var_keys()
            + sample2_original_mod.var_keys()
            + list(sample1_original_mod.varm.keys())
            + list(sample2_original_mod.varm.keys())
        )

        assert original_var_keys == set(concatenated_mod.varm.keys()) | set(
            concatenated_mod.var.columns.tolist()
        )

    # Values from sample2 (which are also not in sample1) should have NA
    non_shared_features = data_sample2.var_names.difference(data_sample1.var_names)
    assert concatenated_data.var["Shared_feat"].loc[non_shared_features].isna().all()

    # Values from sample1 should not have NA, and should be equal to the original values
    var_values = concatenated_data.var["Shared_feat"].loc[data_sample1.var_names]
    data_sample1.var["Shared_feat"].equals(var_values)


def test_concat_different_columns_per_modality(
    run_component, sample_1_h5mu, sample_2_h5mu, write_mudata_to_file, random_h5mu_path
):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that is missing in one modality compared to the other, but the the column
    is missing from the same modalities in both samples.

    Looking at Shared_feat here:

                   mod1      mod2
    sample 1          x   present
    sample 2          x   present
    """
    sample_2_h5mu = remove_annotation_column(
        sample_2_h5mu, ["Shared_feat"], axis="var", modality_name="mod1"
    )
    sample_1_h5mu = remove_annotation_column(
        sample_1_h5mu, ["Shared_feat"], axis="var", modality_name="mod1"
    )

    input_sample1_path = write_mudata_to_file(sample_1_h5mu)
    input_sample2_path = write_mudata_to_file(sample_2_h5mu)

    output_path = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            input_sample1_path,
            "--input",
            input_sample2_path,
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )

    assert Path(output_path).is_file() is True
    concatenated_data = md.read(output_path)

    data_sample1 = md.read(str(input_sample1_path))
    data_sample2 = md.read(str(input_sample2_path))

    # Check if all features are present
    assert (
        concatenated_data.n_vars
        == data_sample1.var.index.union(data_sample2.var.index).size
    )

    for mod_name in ("mod1", "mod2"):
        concatenated_mod = concatenated_data.mod[mod_name]
        data_sample1_mod = data_sample1.mod[mod_name]
        data_sample2_mod = data_sample2.mod[mod_name]
        original_var_keys = set(
            data_sample1_mod.var_keys()
            + data_sample2_mod.var_keys()
            + list(data_sample2_mod.varm.keys())
            + list(data_sample1_mod.varm.keys())
        )

        assert original_var_keys == set(concatenated_mod.varm.keys()) | set(
            concatenated_mod.var.columns.tolist()
        )

    # Check if the shared column stays removed from modality
    assert "Shared_feat" not in concatenated_data.mod["mod1"].var.columns

    # Values from modality 1 have NA
    mod_1_features = data_sample1["mod1"].var_names.union(
        data_sample2["mod1"].var_names
    )
    assert concatenated_data.var.loc[mod_1_features, "mod2:Shared_feat"].isna().all()

    # Values from modalitu should not have NA, and should be equal to the original values
    mod2_data = pd.concat(
        [
            data_sample2["mod2"].var["Shared_feat"],
            data_sample1["mod2"].var["Shared_feat"],
        ]
    )
    mod2_features = mod2_data.index
    assert (
        concatenated_data.var.loc[mod2_features, "mod2:Shared_feat"]
        .astype(str)
        .equals(mod2_data)
    )


def test_concat_different_columns_per_modality_and_per_sample(
    run_component, sample_1_h5mu, sample_2_h5mu, write_mudata_to_file, random_h5mu_path
):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that differ between the modalities and also between samples


    Looking at 'Feat4' from sample 2 here:
                   mod1      mod2
    sample 1          x         x
    sample 2          x   present
    """

    input_sample1_path = write_mudata_to_file(sample_1_h5mu)
    input_sample2_path = write_mudata_to_file(sample_2_h5mu)
    output_path = random_h5mu_path()

    run_component(
        [
            "--input_id",
            "mouse;human",
            "--input",
            input_sample1_path,
            "--input",
            input_sample2_path,
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )

    assert Path(output_path).is_file()
    concatenated_data = md.read(output_path)

    data_sample1 = md.read(input_sample1_path)
    data_sample2 = md.read(input_sample2_path)

    # Check if all features are present
    assert (
        concatenated_data.n_vars
        == data_sample1.var_names.union(data_sample2.var_names).size
    )

    # Check if all features are present
    for mod_name in ("mod1", "mod2"):
        concatenated_mod = concatenated_data.mod[mod_name]
        data_sample1_mod = data_sample1.mod[mod_name]
        data_sample2_mod = data_sample2.mod[mod_name]
        original_var_keys = set(
            data_sample1_mod.var_keys()
            + data_sample2_mod.var_keys()
            + list(data_sample2_mod.varm.keys())
            + list(data_sample1_mod.varm.keys())
        )

        assert original_var_keys == set(
            column_name.removeprefix("conflict_")
            for column_name in concatenated_mod.varm.keys()
        ) | set(concatenated_mod.var.columns.tolist())

    assert "Shared_feat" in concatenated_data.mod["mod2"].var.columns

    # Values from modality 1 have NA
    mod_1_features = data_sample1["mod1"].var_names.union(
        data_sample2["mod1"].var_names
    )
    assert concatenated_data.var.loc[mod_1_features, "mod2:Feat4"].isna().all()

    # Values from modality 2 should not have NA if they originate from sample2
    # These values should be equal to the original values
    mod2_data = data_sample2["mod2"].var["Feat4"].rename("mod2:Feat4")
    mod2_features = mod2_data.index
    assert (
        concatenated_data.var.loc[mod2_features, "mod2:Feat4"]
        .astype(str)
        .equals(mod2_data)
    )

    # Values from modality2 should have NA if they originate from sample1 (and only from sample1)
    non_shared_features = data_sample1.var_names.difference(data_sample2.var_names)
    assert concatenated_data.var.loc[non_shared_features, "mod2:Feat4"].isna().all()


@pytest.mark.parametrize(
    "test_value,test_value_dtype,expected",
    [
        ("bar", "str", "bar"),
        (True, pd.BooleanDtype(), True),
        (1, pd.Int16Dtype(), 1),
        (0.1, float, 0.1),
        (0.1, np.float64, 0.1),
        (np.nan, np.float64, pd.NA),
    ],
)
def test_concat_remove_na(
    run_component,
    sample_1_h5mu,
    sample_2_h5mu,
    write_mudata_to_file,
    random_h5mu_path,
    test_value,
    test_value_dtype,
    expected,
    change_column_contents,
):
    """
    Test concatenation of samples where the column from one sample contains NA values
    NA values should be removed from the concatenated result

                      mod1    mod2
    sample 1            NA      NA
    sample 2    test_value      NA
    """
    change_column_contents(
        sample_1_h5mu, "var", "Shared_feat", {"mod1": np.nan, "mod2": np.nan}
    )
    change_column_contents(
        sample_2_h5mu, "var", "Shared_feat", {"mod1": test_value, "mod2": np.nan}
    )
    sample_2_h5mu.var["Shared_feat"] = sample_2_h5mu.var["Shared_feat"].astype(
        test_value_dtype
    )
    output_path = random_h5mu_path()

    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )

    assert Path(output_path).is_file()
    concatenated_data = md.read(output_path)

    # Values from modality 2 have NA
    mod_2_features = sample_1_h5mu["mod2"].var_names.union(
        sample_2_h5mu["mod2"].var_names
    )
    assert concatenated_data.var.loc[mod_2_features, "Shared_feat"].isna().all()

    # Values from modality 1 should not have NA if they originate from sample 1
    # These values should be equal to the original values
    assert sample_1_h5mu["mod1"].var["Shared_feat"].isna().all()

    # Values from modality 1 should hold a value if they originate from sample 2
    mod1_features = sample_2_h5mu["mod1"].var_names.difference(sample_1_h5mu.var_names)
    if not pd.isna(expected):
        assert (
            concatenated_data.var.loc[mod1_features, "Shared_feat"] == expected
        ).all()
    else:
        assert concatenated_data.var.loc[mod1_features, "Shared_feat"].isna().all()

    # The 'Shared_feat' column for mod1 contains an overlapping feature.
    # For sample 1, it is NA, for sample 2 is is filled with test value.
    # The concat component should choose the test-value over NA
    shared_features = sample_2_h5mu.var_names.intersection(sample_1_h5mu.var_names)
    if not pd.isna(expected):
        assert (
            concatenated_data.var.loc[shared_features, "Shared_feat"] == expected
        ).all()
    else:
        assert concatenated_data.var.loc[shared_features, "Shared_feat"].isna().all()


def test_concat_invalid_h5_error_includes_path(
    run_component, tmp_path, sample_1_h5mu, write_mudata_to_file
):
    empty_file = tmp_path / "empty.h5mu"
    empty_file.touch()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_id",
                "mouse;empty",
                "--input",
                write_mudata_to_file(sample_1_h5mu),
                "--input",
                empty_file,
                "--output",
                "concat.h5mu",
                "--other_axis_mode",
                "move",
            ]
        )
    assert re.search(
        rf"OSError: Failed to load .*{str(empty_file)}\. Is it a valid h5 file?",
        err.value.stdout.decode("utf-8"),
    )


@pytest.mark.parametrize(
    "test_value_1,value_1_dtype,test_value_2,value_2_dtype,expected",
    [
        (1, float, "1", str, pd.CategoricalDtype(categories=["1.0", "1"])),
        (1, np.float64, "1", str, pd.CategoricalDtype(categories=["1.0", "1"])),
        (1, pd.Int16Dtype(), 2.0, pd.Int16Dtype(), pd.Int64Dtype()),
        (True, bool, False, bool, pd.BooleanDtype()),
        (True, pd.BooleanDtype(), False, bool, pd.BooleanDtype()),
        ("foo", str, "bar", str, pd.CategoricalDtype(categories=["bar", "foo"])),
    ],
)
def test_concat_dtypes_per_modality(
    run_component,
    write_mudata_to_file,
    change_column_contents,
    sample_1_h5mu,
    sample_2_h5mu,
    test_value_1,
    value_1_dtype,
    test_value_2,
    value_2_dtype,
    expected,
    random_h5mu_path,
):
    """
    Test joining column with different dtypes to make sure that they are writable.
    The default path is to convert all non-na values to strings and wrap the column into a categorical dtype.
    Here, we test on the level of a single modality only. Because the mod1 modality for both sample 1 and
    sample 2 contain a column 'test_col' and there is an overlapping feature name (overlapping_var_mod1),
    there is a conflict for this var column in mod 1 for this column. Upon concatenation, the column is moved
    to .varm, but for mod1 only. The column is concatenated for mod2 as planned. Here we check if the results
    for the test column in mod2 is still writable.
    """
    change_column_contents(
        sample_1_h5mu, "var", "test_col", {"mod1": test_value_1, "mod2": test_value_1}
    )
    sample_1_h5mu.var["test_col"] = sample_1_h5mu.var["test_col"].astype(value_1_dtype)
    change_column_contents(
        sample_2_h5mu, "var", "test_col", {"mod1": test_value_2, "mod2": test_value_2}
    )
    sample_2_h5mu.var["test_col"] = sample_2_h5mu.var["test_col"].astype(value_2_dtype)

    output_file = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_file,
            "--other_axis_mode",
            "move",
        ]
    )
    concatenated_data = md.read(output_file)
    assert concatenated_data["mod2"].var["test_col"].dtype == expected


@pytest.mark.parametrize(
    "test_value,value_dtype,expected",
    [
        (1, float, pd.Int64Dtype()),
        (1, np.float64, pd.Int64Dtype()),
        (1, pd.Int16Dtype(), pd.Int16Dtype()),
        (True, bool, pd.BooleanDtype()),
        (True, pd.BooleanDtype(), pd.BooleanDtype()),
        ("foo", str, pd.CategoricalDtype(categories=["foo"])),
    ],
)
def test_concat_dtypes_per_modality_multidim(
    run_component,
    write_mudata_to_file,
    sample_1_h5mu,
    sample_2_h5mu,
    test_value,
    value_dtype,
    expected,
    random_h5mu_path,
):
    """
    Test if the result of concatenation is still writable when the input already contain
    data in .varm and this data is kept. Because we are joining observations, the dtype of this
    data may change and the result might not be writable anymore
    """

    sample_1_h5mu["mod1"].varm["test_df"] = pd.DataFrame(
        index=sample_1_h5mu["mod1"].var_names
    )
    sample_1_h5mu["mod1"].varm["test_df"]["test_col"] = test_value
    sample_1_h5mu["mod1"].varm["test_df"]["test_col"] = (
        sample_1_h5mu["mod1"].varm["test_df"]["test_col"].astype(value_dtype)
    )

    output_file = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_file,
            "--other_axis_mode",
            "move",
        ]
    )
    concatenated_data = md.read(output_file)
    assert concatenated_data["mod1"].varm["test_df"]["test_col"].dtype == expected


@pytest.mark.parametrize(
    "test_value_1,test_value_2,expected",
    [(1, "1", pd.CategoricalDtype(categories=["1.0", "1"]))],
)
def test_concat_dtypes_global(
    run_component,
    write_mudata_to_file,
    change_column_contents,
    sample_1_h5mu,
    sample_2_h5mu,
    test_value_1,
    test_value_2,
    expected,
    random_h5mu_path,
):
    """
    Test joining column with different dtypes to make sure that they are writable.
    The default path is to convert all non-na values to strings and wrap the column into a categorical dtype.
    Here, we test on the level of a column that is added to a global annotation matrix.
    """
    change_column_contents(
        sample_1_h5mu, "var", "test_col", {"mod1": test_value_1, "mod2": test_value_1}
    )
    change_column_contents(
        sample_2_h5mu, "var", "test_col", {"mod1": test_value_2, "mod2": test_value_2}
    )
    sample1_mod1_names = sample_2_h5mu["mod1"].var_names
    # Here, we avoid a conflict between sample 1 and sample 2 by making sure there is no overlap in features
    # between sample 1 and sample 2 (no shared var_names). If this change would not be done, a different
    # value for sample 1 and sample 2 would be found by the concat component for the var feature
    # 'overlapping_var_mod1' for modality 'mod1'. The concat component would move the column for mod1 to
    # .varm because of this conflict, and in the global .var column of the concatenated object, only
    # a 'mod2:test_col' column would be present. But here, we want to test the column that is populated by
    # both 'mod1' and 'mod2'
    assert "overlapping_var_mod1" in sample1_mod1_names
    new_names = sample1_mod1_names.where(
        ~sample1_mod1_names.isin(["overlapping_var_mod1"]), "non_overlapping"
    )
    sample_2_h5mu["mod1"].var_names = new_names
    sample_2_h5mu.update()
    output_file = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_file,
            "--other_axis_mode",
            "move",
        ]
    )
    concatenated_data = md.read(output_file)
    assert concatenated_data.var["test_col"].dtype == expected


def test_non_overlapping_modalities(
    run_component, sample_2_h5mu, sample_3_h5mu, random_h5mu_path, write_mudata_to_file
):
    """
    Test that the component does not fail when the modalities are not shared between samples.
    """
    output_path = random_h5mu_path()
    input_file_2 = write_mudata_to_file(sample_2_h5mu)
    input_file_3 = write_mudata_to_file(sample_3_h5mu)

    run_component(
        [
            "--input_id",
            "sample2;sample3",
            "--input",
            input_file_2,
            "--input",
            input_file_3,
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )
    output_data = md.read(output_path)
    assert set(output_data.mod.keys()) == {"mod1", "mod2", "mod3"}


def test_resolve_annotation_conflict_missing_column(
    run_component,
    sample_1_h5mu,
    sample_2_h5mu,
    sample_3_h5mu,
    write_mudata_to_file,
    random_h5mu_path,
):
    """
    Test using mode 'move' and resolving a conflict in metadata between the samples,
    but the metadata column is missing in one of the samples.
    """
    output_path = random_h5mu_path()
    input_file_1 = write_mudata_to_file(sample_1_h5mu)
    input_file_2 = write_mudata_to_file(sample_2_h5mu)
    input_file_3 = write_mudata_to_file(sample_3_h5mu)

    run_component(
        [
            "--input_id",
            "sample1;sample2;sample3",
            "--input",
            input_file_1,
            "--input",
            input_file_2,
            "--input",
            input_file_3,
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )

    concatenated_data = md.read(output_path)
    # 'Shared_feat' is defined for mod1 in sample 1 and 2 and there is a conflict
    assert "conflict_Shared_feat" in concatenated_data["mod1"].varm
    # 'Shared_feat' is defined for mod2 in sample 1 and 2 and there is no conflict
    assert "Shared_feat" in concatenated_data["mod2"].var.columns
    # 'Shared_feat' is not defined in any of the samples samples for modality 3
    assert "Shared_feat" not in concatenated_data["mod3"].var.columns
    assert "Shared_feat" not in concatenated_data["mod3"].varm


def test_mode_move(
    run_component, sample_1_h5mu, sample_2_h5mu, random_h5mu_path, write_mudata_to_file
):
    """
    Test that in case of a conflict, the conflicting columns are move to the multidimensional annotation slot
    (.varm and .obsm). The key of the datafame in the slot should start with 'conflict_' followed by the name
    of the column and the columns of the dataframe should contain the sample names.
    """
    output_path = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )
    assert output_path.is_file()
    concatenated_data = md.read(output_path)

    # Check if observations from all of the samples are present
    assert concatenated_data.n_obs == sample_1_h5mu.n_obs + sample_2_h5mu.n_obs

    # Check if all modalities are present
    sample1_mods, sample2_mods = (
        set(sample_1_h5mu.mod.keys()),
        set(sample_2_h5mu.mod.keys()),
    )
    concatentated_mods = set(concatenated_data.mod.keys())
    assert (sample1_mods | sample2_mods) == concatentated_mods

    varm_check = {
        "mod1": ({"conflict_Shared_feat": ("sample1", "sample2")}),
        "mod2": {},
    }

    # Check if all features are present
    for mod_name in ("mod1", "mod2"):
        concatenated_mod = concatenated_data.mod[mod_name]
        sample_1_mod = sample_1_h5mu.mod[mod_name]
        sample_2_mod = sample_2_h5mu.mod[mod_name]
        original_varm_keys = set(
            list(sample_1_mod.varm.keys()) + list(sample_2_mod.varm.keys())
        )
        original_var_keys = (
            set(sample_1_mod.var_keys() + sample_2_mod.var_keys()) | original_varm_keys
        )

        assert original_var_keys == set(
            column_name.removeprefix("conflict_")
            for column_name in concatenated_mod.varm.keys()
        ) | set(concatenated_mod.var.columns.tolist())

        varm_expected = varm_check[mod_name]
        assert set(concatenated_mod.varm.keys()) == set(
            varm_expected.keys() | original_varm_keys
        )
        for varm_key, expected_columns in varm_expected.items():
            assert tuple(concatenated_mod.varm[varm_key].columns) == expected_columns
        if not varm_expected:
            assert set(concatenated_mod.varm.keys()) == original_varm_keys
        assert concatenated_mod.obsm == {}


# Execute this test multiple times, anndata.concat sometimes returns the observations in a different order
@pytest.mark.parametrize("_", range(10))
def test_concat_var_obs_names_order(
    run_component,
    sample_1_h5mu,
    sample_2_h5mu,
    write_mudata_to_file,
    random_h5mu_path,
    _,
):
    """
    Test that the var_names and obs_names are still linked to the correct count data.
    """
    output_path = random_h5mu_path()
    sample_1_h5mu["mod1"].obs["sample_id"] = "sample1"
    sample_1_h5mu["mod2"].obs["sample_id"] = "sample1"
    sample_2_h5mu["mod1"].obs["sample_id"] = "sample2"
    sample_2_h5mu["mod2"].obs["sample_id"] = "sample2"
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
        ]
    )
    assert output_path.is_file()
    for sample_name, sample_h5mu in {
        "sample1": sample_1_h5mu,
        "sample2": sample_2_h5mu,
    }.items():
        for mod_name in ["mod1", "mod2"]:
            data_sample = sample_h5mu[mod_name].copy()
            processed_data_ad = md.read_h5ad(output_path, mod=mod_name)
            processed_data_ad = processed_data_ad[
                processed_data_ad.obs["sample_id"] == sample_name
            ]
            processed_data_ad = processed_data_ad[:, data_sample.var_names]
            processed_data = pd.DataFrame(
                processed_data_ad.X,
                index=processed_data_ad.obs_names,
                columns=processed_data_ad.var_names,
            )
            data_sample = pd.DataFrame(
                data_sample.X,
                index=data_sample.obs_names,
                columns=data_sample.var_names,
            ).reindex_like(processed_data)
            pd.testing.assert_frame_equal(
                processed_data, data_sample, check_dtype=False
            )


def test_keep_uns(
    run_component, sample_1_h5mu, sample_2_h5mu, write_mudata_to_file, random_h5mu_path
):
    sample_1_h5mu.uns["global_uns_sample1"] = "dolor"
    sample_1_h5mu.uns["overlapping_global"] = "sed"
    sample_2_h5mu.uns["global_uns_sample2"] = "amet"
    sample_2_h5mu.uns["overlapping_global"] = "elit"
    output_path = random_h5mu_path()
    run_component(
        [
            "--input_id",
            "sample1;sample2",
            "--input",
            write_mudata_to_file(sample_1_h5mu),
            "--input",
            write_mudata_to_file(sample_2_h5mu),
            "--output",
            output_path,
            "--other_axis_mode",
            "move",
            "--uns_merge_mode",
            "make_unique",
        ]
    )
    assert output_path.is_file()
    concatenated_data = md.read(output_path)
    mod1 = concatenated_data.mod["mod1"]
    mod2 = concatenated_data.mod["mod2"]
    assert set(concatenated_data.uns.keys()) == set(
        [
            "global_uns_sample1",
            "global_uns_sample2",
            "sample1_overlapping_global",
            "sample2_overlapping_global",
        ]
    )
    assert set(mod1.uns.keys()) == set(
        [
            "sample1_overlapping_uns_key",
            "uns_unique_to_sample1",
            "sample2_overlapping_uns_key",
            "uns_unique_to_sample2",
        ]
    )
    assert set(mod2.uns.keys()) == set()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
