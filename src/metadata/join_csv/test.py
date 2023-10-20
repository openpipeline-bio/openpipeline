import sys
import pytest
from tempfile import NamedTemporaryFile
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
import subprocess


## VIASH START
meta = {
    'executable': './target/docker/integrate/add_metadata/add_metadata',
}
## VIASH END

@pytest.fixture
def sample_h5mu(tmp_path):
    output = tmp_path / "output.h5mu"
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A", "sample1"], ["B", "sample2"]], index=df.index, columns=["Obs", "sample_id"])
    var = pd.DataFrame([["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
                        index=df.columns, columns=["Feat", "sample_id_var"])
    ad1 = AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)
    mudata = MuData({'mod1': ad1, 'mod2': ad2})

    mudata.write_h5mu(output)
    return output, mudata

def test_add_metadata_var(run_component, tmp_path, sample_h5mu):
    input_h5mu, input_mudata = sample_h5mu
    input_csv = tmp_path / "input.csv"
    output_h5mu = tmp_path / "output.h5mu"

    # create csv
    csv = pd.DataFrame({"id": ["sample1", "sample2"], "foo": ["v", "w"], "bar": ["x", "y"]})
    csv.to_csv(str(input_csv), index=False)

    run_component([
        "--input", str(input_h5mu),
        "--input_csv", str(input_csv),
        "--output", str(output_h5mu),
        "--modality", "mod1",
        "--var_key", "sample_id_var",
        "--csv_key", "id",
        "--compression_output", "gzip"
    ])
    result = read_h5mu(output_h5mu)
    expected_var = pd.DataFrame(
        {
            "Feat": ["a", "b", "c"],
            "sample_id_var": ["sample1", "sample2", "sample1"],
            "foo": ["v", "w", "v"],
            "bar": ["x", "y", "x"]
        },
        index=pd.Index(['var1', 'var2', 'var3']),
    ).astype(
        {
            "Feat": "object",
            "sample_id_var": "category",
            "foo": "category",
            "bar": "category"
        }
    )
    pd.testing.assert_frame_equal(result.mod['mod1'].var, expected_var)
    pd.testing.assert_frame_equal(result.mod['mod1'].obs, input_mudata.mod['mod1'].obs)
    pd.testing.assert_frame_equal(result.mod['mod2'].obs, input_mudata.mod['mod2'].obs)
    pd.testing.assert_frame_equal(result.mod['mod2'].var, input_mudata.mod['mod2'].var)

def test_add_metadata_matrix_sample_column(run_component, tmp_path, sample_h5mu):
    input_h5mu, input_mudata = sample_h5mu
    input_csv = tmp_path / "input.csv"
    output_h5mu = tmp_path / "output.h5mu"

    # create csv
    csv = pd.DataFrame({"id": ["sample1", "sample2"], "foo": ["v", "w"], "bar": ["x", "y"]})
    csv.to_csv(str(input_csv), index=False)

    run_component([
        "--input", str(input_h5mu),
        "--input_csv", str(input_csv),
        "--output", str(output_h5mu),
        "--modality", "mod1",
        "--obs_key", "sample_id",
        "--csv_key", "id",
    ])

    result = read_h5mu(output_h5mu)
    expected_obs = pd.DataFrame(
        {
            "Obs": ["A", "B"],
            "sample_id": ["sample1", "sample2"],
            "foo": ["v", "w"],
            "bar": ["x", "y"]
        },
        index=pd.Index(['obs1', 'obs2']),
    ).astype(
        {
            "Obs": "object",
            "foo": "object",
            "bar": "object"
        }
    )
    pd.testing.assert_frame_equal(result.mod['mod1'].obs, expected_obs)
    pd.testing.assert_frame_equal(result.mod['mod1'].var, input_mudata.mod['mod1'].var)
    pd.testing.assert_frame_equal(result.mod['mod2'].obs, input_mudata.mod['mod2'].obs)
    pd.testing.assert_frame_equal(result.mod['mod2'].var, input_mudata.mod['mod2'].var)

def test_add_not_all_samples_in_csv_raises(run_component, tmp_path, sample_h5mu):
    input_h5mu, input_mudata = sample_h5mu
    input_csv = tmp_path / "input.csv"
    output_h5mu = tmp_path / "output.h5mu"

    csv = pd.DataFrame({"id": ["sample1", "lorem"], "foo": ["v", "w"], "bar": ["x", "y"]})
    csv.to_csv(str(input_csv), index=False)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", str(input_h5mu),
            "--input_csv", str(input_csv),
            "--output", str(output_h5mu),
            "--modality", "mod1",
            "--obs_key", "sample_id",
            "--csv_key", "id",
        ])
        assert not output_h5mu.is_file(), "No output should be created."
        assert "Not all sample IDs selected from obs (using the column selected " \
            "with --var_key or --obs_key) were found in the csv file." in err.stdout.decode('utf-8')

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))