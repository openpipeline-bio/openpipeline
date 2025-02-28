import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/metadata/join_uns_to_obs/join_uns_to_obs",
    "config": "./src/metadata/join_uns_to_obs/config.vsh.yml",
}
## VIASH END


@pytest.fixture
def ad_w_uns():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame(
        [["A", "sample1"], ["B", "sample2"]],
        index=df.index,
        columns=["Obs", "sample_id"],
    )
    var = pd.DataFrame(
        [["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
        index=df.columns,
        columns=["Feat", "sample_id_var"],
    )
    obsm = pd.DataFrame(
        [["X", "W"]], index=pd.Index([0]), columns=["uns_col1", "uns_col2"]
    )
    return AnnData(df, obs=obs, var=var, uns={"obsm1": obsm})


@pytest.fixture
def ad_wo_uns():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    var = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    return AnnData(df, obs=obs, var=var)


@pytest.fixture
def sample_h5mu(ad_w_uns, ad_wo_uns):
    mudata = MuData({"mod1": ad_w_uns, "mod2": ad_wo_uns})
    return mudata


def test_join_uns_to_obs(
    run_component, random_h5mu_path, write_mudata_to_file, sample_h5mu
):
    input_file = write_mudata_to_file(sample_h5mu)
    output_file = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_file),
            "--modality",
            "mod1",
            "--uns_key",
            "obsm1",
            "--output",
            str(output_file),
        ]
    )

    expected_obs = pd.DataFrame(
        {
            "Obs": ["A", "B"],
            "sample_id": ["sample1", "sample2"],
            "uns_col1": ["X", "X"],
            "uns_col2": ["W", "W"],
        },
        index=pd.Index(["obs1", "obs2"]),
    ).astype(
        {
            "Obs": "object",
            "sample_id": "object",
            "uns_col1": "category",
            "uns_col2": "category",
        }
    )

    assert output_file.is_file()
    output_mudata = read_h5mu(output_file)
    assert "obsm1" in output_mudata.mod["mod1"].uns

    sample_h5mu.mod["mod1"].obs = expected_obs
    assert_annotation_objects_equal(sample_h5mu, output_mudata)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
