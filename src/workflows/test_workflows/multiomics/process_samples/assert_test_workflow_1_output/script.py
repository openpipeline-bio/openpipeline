import sys
import mudata as mu
from pathlib import Path
import shutil
import pandas as pd
import os
import pytest


##VIASH START
par = {
    "input": "work/1f/992ff0941265de4e4d49e5a287f41b/_viash_par/input_1/merged.umap.output.h5mu",
    "orig_input": [
        "work/1f/992ff0941265de4e4d49e5a287f41b/_viash_par/orig_input_1/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
        "work/1f/992ff0941265de4e4d49e5a287f41b/_viash_par/orig_input_2/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    ],
}

meta = {"resources_dir": "src/base"}
sys.path.append("src/utils")
sys.path.append("src/base")
##VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def sample_1_h5mu():
    return mu.read_h5mu(par["orig_input"][0])


@pytest.fixture
def sample_2_h5mu():
    return mu.read_h5mu(par["orig_input"][1])


@pytest.fixture
def output_h5mu():
    return mu.read_h5mu(par["input"])


def test_modalities_present(output_h5mu):
    assert all(k in output_h5mu.obsm.keys() for k in ("rna", "rna"))


def test_rna_embeddings_present(output_h5mu):
    assert all(k in output_h5mu["rna"].obsm.keys() for k in ("X_pca", "X_umap"))
    assert "pca_loadings" in output_h5mu["rna"].varm.keys()
    assert "pca_variance" in output_h5mu["rna"].uns.keys()


def test_rna_layers(output_h5mu):
    assert "log_normalized" in output_h5mu["rna"].layers.keys()


def test_rna_raw_expressions(output_h5mu, sample_1_h5mu, sample_2_h5mu):
    rna_output = output_h5mu["rna"]
    sample_1_output = rna_output[
        rna_output.obs["sample_id"] == "pbmc", sample_1_h5mu["rna"].var_names
    ]
    obs_index = pd.Index(
        sample_1_output.obs_names.to_series().str.extract("pbmc_(.+)", expand=False)
    )
    sample_1_output = pd.DataFrame(
        sample_1_output.X.toarray(), index=obs_index, columns=sample_1_output.var_names
    )

    sample_2_output = rna_output[
        rna_output.obs["sample_id"] == "pbmc_with_more_args",
        sample_2_h5mu["rna"].var_names,
    ]
    obs_index = pd.Index(
        sample_2_output.obs_names.to_series().str.extract(
            "pbmc_with_more_args_(.+)", expand=False
        )
    )
    sample_2_output = pd.DataFrame(
        sample_2_output.X.toarray(), index=obs_index, columns=sample_2_output.var_names
    )
    for input_counts, output_counts in [
        (sample_1_h5mu["rna"], sample_1_output),
        (sample_2_h5mu["rna"], sample_2_output),
    ]:
        input_counts_df = pd.DataFrame(
            input_counts.X.toarray(),
            index=input_counts.obs_names,
            columns=input_counts.var_names,
        )
        pd.testing.assert_frame_equal(input_counts_df, output_counts)


def test_rna_uns_slots(output_h5mu):
    assert all(
        k in output_h5mu["rna"].uns.keys()
        for k in (
            "neighbors",
            "pbmc_metrics_cellranger",
            "pbmc_with_more_args_metrics_cellranger",
            "pca_variance",
        )
    )


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    shutil.copyfile(
        os.path.join(meta["resources_dir"], "openpipelinetestutils", "conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))
