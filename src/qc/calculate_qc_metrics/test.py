import sys
import pytest
from pathlib import Path
import mudata as md
import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import scipy
import uuid
from pandas.testing import assert_series_equal
from itertools import product
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/qc/calculate_qc_metrics/calculate_qc_metrics",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "./src/qc/calculate_qc_metrics/config.vsh.yaml",
    "cpus": 2,
}
## VIASH END


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_mudata_random():
    rng = np.random.default_rng(seed=1)
    random_counts = scipy.sparse.random(
        50000, 100, density=0.8, format="csr", dtype=np.uint32, random_state=rng
    )
    good_dtype = random_counts.astype(np.float32)
    del random_counts
    mod1 = ad.AnnData(
        X=good_dtype,
        obs=pd.DataFrame(index=pd.RangeIndex(50000)),
        var=pd.DataFrame(index=pd.RangeIndex(100)),
    )
    return md.MuData({"mod1": mod1})


@pytest.fixture
def input_mudata(input_path):
    mudata = md.read_h5mu(input_path)
    # create a less sparse matrix to increase the variability in qc statistics
    rng = np.random.default_rng()
    random_counts = scipy.sparse.random(
        *mudata["rna"].X.shape, density=0.8, format="csr", random_state=rng
    )
    mudata["rna"].X = random_counts
    return mudata


@pytest.fixture
def input_mudata_path(tmp_path, input_mudata):
    output_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
    input_mudata.write(output_path)
    return output_path


@pytest.fixture(
    params=product([True, False, np.nan, "random"], ["bool", pd.BooleanDtype()])
)
def mudata_with_boolean_column(tmp_path, input_mudata, request):
    requested_value, requested_type = request.param
    input_var = input_mudata.mod["rna"].var
    input_var["custom"] = requested_value
    if requested_value == "random":
        input_var["custom"] = np.random.choice(
            [True, False], len(input_var), p=[0.20, 0.80]
        )
    input_var["custom"] = input_var["custom"].astype(requested_type)
    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_mudata.write(new_input_path)
    return new_input_path


def test_add_qc(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--modality",
            "rna",
            "--top_n_vars",
            "10;20;90",
            "--output_compression",
            "gzip",
        ]
    )

    assert output_path.is_file()
    data_with_qc = md.read(output_path)
    expected_obs_columns = [
        f"pct_of_counts_in_top_{top_n_var}_vars" for top_n_var in ("10", "20", "90")
    ] + ["total_counts", "num_nonzero_vars"]
    expected_var_columns = [
        "pct_dropout",
        "num_nonzero_obs",
        "obs_mean",
        "total_counts",
    ]
    data_with_qc.mod["rna"].var = data_with_qc.mod["rna"].var.drop(
        columns=expected_var_columns, errors="raise"
    )
    data_with_qc.var = data_with_qc.var.drop(
        columns=list(map("rna:{}".format, expected_var_columns)), errors="errors"
    )
    data_with_qc.mod["rna"].obs = data_with_qc.mod["rna"].obs.drop(
        columns=expected_obs_columns, errors="raise"
    )
    data_with_qc.obs = data_with_qc.obs.drop(
        columns=list(map("rna:{}".format, expected_obs_columns)), errors="errors"
    )
    assert_annotation_objects_equal(input_path, data_with_qc)


@pytest.mark.parametrize(
    "optional_parameter,annotation_matrix,arg_value",
    [
        ("--output_obs_num_nonzero_vars", "obs", "lorem"),
        ("--output_obs_total_counts_vars", "obs", "ipsum"),
        ("--output_var_num_nonzero_obs", "var", "dolor"),
        ("--output_var_total_counts_obs", "var", "amet"),
        ("--output_var_obs_mean", "var", "sit"),
        ("--output_var_pct_dropout", "var", "elit"),
    ],
)
def test_qc_metrics_set_output_column(
    run_component,
    input_mudata_path,
    optional_parameter,
    annotation_matrix,
    arg_value,
    random_h5mu_path,
):
    output_h5mu = random_h5mu_path()
    args = [
        "--input",
        input_mudata_path,
        "--output",
        output_h5mu,
        "--modality",
        "rna",
        "--output_compression",
        "gzip",
        optional_parameter,
        arg_value,
    ]

    run_component(args)
    default_obs_columns = {
        "output_obs_num_nonzero_vars": "num_nonzero_vars",
        "output_obs_total_counts_vars": "total_counts",
    }
    default_var_columns = {
        "output_var_num_nonzero_obs": "num_nonzero_obs",
        "output_var_total_counts_obs": "total_counts",
        "output_var_obs_mean": "obs_mean",
        "output_var_pct_dropout": "pct_dropout",
    }
    defaults = {"var": default_var_columns, "obs": default_obs_columns}
    defaults[annotation_matrix].update({optional_parameter.strip("-"): arg_value})
    assert output_h5mu.is_file()
    data_with_qc = md.read(output_h5mu)
    for attribute_name in ("var", "obs"):
        setattr(
            data_with_qc.mod["rna"],
            attribute_name,
            getattr(data_with_qc.mod["rna"], attribute_name).drop(
                columns=list(defaults[attribute_name].values()), errors="raise"
            ),
        )
        setattr(
            data_with_qc,
            attribute_name,
            getattr(data_with_qc, attribute_name).drop(
                columns=list(
                    map("rna:{}".format, list(defaults[attribute_name].values()))
                ),
                errors="raise",
            ),
        )
    assert_annotation_objects_equal(input_mudata_path, data_with_qc)


@pytest.mark.parametrize(
    "optional_parameter,annotation_matrix,expected_missing,",
    [
        ("--var_qc_metrics", "obs", "total_counts_.*|pct_*"),
        ("--top_n_vars", "obs", "pct_of_counts_in_top_.*"),
        ("--output_obs_num_nonzero_vars", "obs", "num_nonzero_vars"),
        ("--output_obs_total_counts_vars", "obs", "total_counts"),
        ("--output_var_num_nonzero_obs", "var", "num_nonzero_obs"),
        ("--output_var_total_counts_obs", "var", "total_counts"),
        ("--output_var_obs_mean", "var", "obs_mean"),
        ("--output_var_pct_dropout", "var", "pct_dropout"),
    ],
)
def test_qc_metrics_optional(
    run_component,
    input_mudata_path,
    optional_parameter,
    annotation_matrix,
    expected_missing,
):
    args = [
        "--input",
        input_mudata_path,
        "--output",
        "foo.h5mu",
        "--modality",
        "rna",
        "--output_compression",
        "gzip",
    ]
    if optional_parameter not in ["--var_qc_metrics", "--top_n_vars"]:
        args.extend([optional_parameter, ""])

    run_component(args)
    assert Path("foo.h5mu").is_file()
    data_with_qc = md.read("foo.h5mu")
    matrix = getattr(data_with_qc.mod["rna"], annotation_matrix)
    assert matrix.filter(regex=expected_missing, axis=1).empty


def test_calculcate_qc_var_qc_metrics(
    run_component, mudata_with_boolean_column, random_h5mu_path
):
    output_path = random_h5mu_path()
    input_data = md.read_h5mu(mudata_with_boolean_column)
    args = [
        "--input",
        str(mudata_with_boolean_column),
        "--output",
        str(output_path),
        "--modality",
        "rna",
        "--top_n_vars",
        "10;20;90",
        "--var_qc_metrics",
        "custom",
    ]
    if input_data.mod["rna"].var["custom"].isna().any():
        args.extend(["--var_qc_metrics_fill_na_value", "True"])

    run_component(args)
    assert output_path.is_file()
    data_with_qc = md.read(output_path)["rna"]
    for qc_metric in ("pct_custom", "total_counts_custom"):
        assert qc_metric in data_with_qc.obs
    # Do a percentage calculation based on indexes
    # and compare it to the calculations from the component
    custom_column_true = data_with_qc.var["custom"].fillna(True)
    gene_counts = sc.get.obs_df(data_with_qc, keys=data_with_qc.var_names.to_list())
    gene_counts_custom = gene_counts.loc[:, custom_column_true]
    sum_custom_column = gene_counts_custom.sum(axis=1)
    sum_all = gene_counts.sum(axis=1)
    percentage = (sum_custom_column / sum_all) * 100
    pd.testing.assert_series_equal(
        percentage,
        data_with_qc.obs["pct_custom"],
        check_exact=False,  # Comparing floats
        check_names=False,
    )
    assert (data_with_qc.obs["pct_custom"] <= 100).all()


def test_compare_scanpy(
    run_component, mudata_with_boolean_column, input_mudata, random_h5mu_path
):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            str(mudata_with_boolean_column),
            "--output",
            str(output_path),
            "--modality",
            "rna",
            "--top_n_vars",
            "10;20;90",
            "--var_qc_metrics",
            "custom",
            "--var_qc_metrics_fill_na_value",
            "False",
        ]
    )
    assert output_path.is_file()

    component_data = md.read(output_path)
    rna_mod = component_data.mod["rna"]

    # Replicate --var_qc_metrics_fill_na_value False
    # Scanpy also does not work with pd.BooleanDtype()
    # So cast from 'boolean' to 'bool'
    input_mudata.mod["rna"].var["custom"] = (
        input_mudata.mod["rna"].var["custom"].fillna(False).astype("bool")
    )
    sc.pp.calculate_qc_metrics(
        input_mudata.mod["rna"],
        expr_type="counts",
        var_type="genes",
        qc_vars=["custom"],
        percent_top=[10, 20, 90],
        use_raw=False,
        inplace=True,
        log1p=False,
    )
    scanpy_var = input_mudata.mod["rna"].var
    component_var = rna_mod.var

    vars_to_compare = {
        "pct_dropout": "pct_dropout_by_counts",
        "num_nonzero_obs": "n_cells_by_counts",
        "obs_mean": "mean_counts",
        "total_counts": "total_counts",
    }
    for from_var, to_var in vars_to_compare.items():
        assert_series_equal(
            component_var[from_var],
            scanpy_var[to_var],
            check_names=False,
            check_dtype=False,
        )

    scanpy_obs = input_mudata.mod["rna"].obs
    component_obs = rna_mod.obs
    obs_to_compare = {
        "num_nonzero_vars": "n_genes_by_counts",
        "pct_custom": "pct_counts_custom",
        "total_counts_custom": "total_counts_custom",
        "total_counts": "total_counts",
    }
    obs_to_compare |= {
        f"pct_of_counts_in_top_{i}_vars": f"pct_counts_in_top_{i}_genes"
        for i in (10, 20, 90)
    }
    for from_obs, to_obs in obs_to_compare.items():
        assert_series_equal(
            component_obs[from_obs],
            scanpy_obs[to_obs],
            check_names=False,
            check_dtype=False,
        )


def test_total_counts_less_precision_dtype(
    run_component, input_mudata_random, random_h5mu_path
):
    input_path = random_h5mu_path()
    input_mudata_random.write(input_path)
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--modality",
            "mod1",
        ]
    )
    output_data = md.read_h5mu(output_path)
    matrix_good_type = input_mudata_random.mod["mod1"].X
    var_names = input_mudata_random.var_names
    obs_names = input_mudata_random.obs_names
    del input_mudata_random
    input_df = pd.DataFrame(
        matrix_good_type.todense(), columns=var_names, index=obs_names
    )
    total_sums_manual = input_df.to_numpy().sum(axis=0, dtype=np.float128)
    total_counts = output_data.mod["mod1"].var["total_counts"]
    np.testing.assert_allclose(total_sums_manual, total_counts.to_numpy())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
