import subprocess
import sys
import os
import pytest
import re
import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal
import scvi

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"


@pytest.fixture
def create_scvi_model(random_path, tmp_path):
    def wrapper(input_file, reference_file):
        input_data = mu.read_h5mu(input_file)
        input_modality = input_data.mod["rna"]
        reference_data = mu.read_h5mu(reference_file)
        reference_modality = reference_data.mod["rna"]

        reference_data.var["gene_symbol"] = list(reference_data.var.index)
        reference_data.var.index = [
            re.sub("\\.[0-9]+$", "", s) for s in reference_data.var["ensemblid"]
        ]
        reference_modality.var["gene_symbol"] = list(reference_modality.var.index)
        reference_modality.var.index = [
            re.sub("\\.[0-9]+$", "", s) for s in reference_modality.var["ensemblid"]
        ]

        common_ens_ids = list(
            set(reference_modality.var.index).intersection(
                set(input_modality.var.index)
            )
        )

        reference = reference_modality[:, common_ens_ids].copy()
        query = input_modality[:, common_ens_ids].copy()

        scvi.model.SCVI.setup_anndata(reference, labels_key="cell_ontology_class")

        scvi_model = scvi.model.SCVI(
            reference,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=1,
        )
        scvi_model.train(max_epochs=10)

        input_data.mod["rna"] = query
        # reference_data.mod["rna"] = reference

        input_data_file = random_path(extension="h5mu")
        # reference_file = random_path(extension="h5mu")
        scvi_model_file = tmp_path

        input_data.write_h5mu(input_data_file)
        # reference_data.write_h5mu(reference_file)
        scvi_model.save(scvi_model_file, save_anndata=True, overwrite=True)

        return scvi_model_file, input_data_file

    return wrapper


@pytest.fixture
def create_scanvi_model(create_scvi_model, tmp_path):
    def scanvi_wrapper():
        scvi_model_file, input_data_file = create_scvi_model(input_file, reference_file)

        scvi_model = scvi.model.SCVI.load(scvi_model_file)
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category="Unkown",
            labels_key="cell_ontology_class",
        )
        scanvi_model.train(max_epochs=10)

        scanvi_model_file = tmp_path
        scanvi_model.save(scanvi_model_file, save_anndata=True, overwrite=True)

        return scanvi_model_file, input_data_file

    return scanvi_wrapper


def test_simple_execution(run_component, random_h5mu_path, create_scvi_model):
    scvi_model_file, input_file_scvi = create_scvi_model(input_file, reference_file)
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file_scvi,
            "--scvi_reference_model",
            scvi_model_file,
            "--reference_max_epochs",
            "10",
            "--query_max_epochs",
            "10",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file_scvi)
    output_mudata = mu.read_h5mu(output_file)

    assert (
        input_mudata.mod["rna"].n_obs == output_mudata.mod["rna"].n_obs
    ), "Number of observations changed"
    assert (
        input_mudata.mod["rna"].n_vars == output_mudata.mod["rna"].n_vars
    ), "Number of variables changed"
    assert (
        "scanvi_embedding" in output_mudata.mod["rna"].obsm.keys()
    ), "Latent representation not added"
    assert "scanvi_pred" in output_mudata.mod["rna"].obs.keys(), "Predictions not added"
    assert (
        "scanvi_probability" in output_mudata.mod["rna"].obs.keys()
    ), "Probabilities not added"

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])


def test_multiple_arguments(
    run_component, random_h5mu_path, create_scvi_model, tmp_path
):
    scvi_model_file, input_file_scvi = create_scvi_model(input_file, reference_file)
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file_scvi,
            "--scvi_reference_model",
            scvi_model_file,
            "--output",
            output_file,
            "--reference_max_epochs",
            "10",
            "--reference_reduce_lr_on_plateau",
            "True",
            "--reference_lr_patience",
            "5",
            "--reference_lr_factor",
            "0.5",
            "--reference_train_size",
            "0.8",
            "--reference_early_stopping",
            "True",
            "--reference_early_stopping_patience",
            "5",
            "--reference_early_stopping_min_delta",
            "0.01",
            "--query_max_epochs",
            "10",
            "--query_reduce_lr_on_plateau",
            "True",
            "--query_lr_patience",
            "5",
            "--query_lr_factor",
            "0.5",
            "--query_train_size",
            "0.8",
            "--query_early_stopping",
            "True",
            "--query_early_stopping_patience",
            "5",
            "--query_early_stopping_min_delta",
            "0.01",
            "--output_obs_predictions",
            "scanvi_pred",
            "--output_obs_probabilities",
            "scanvi_probabilitity",
            "--output_compression",
            "gzip",
            "--output_model",
            tmp_path,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"
    assert os.path.exists(tmp_path / "model.pt"), "Model file does not exist"

    input_mudata = mu.read_h5mu(input_file_scvi)
    output_mudata = mu.read_h5mu(output_file)

    assert (
        input_mudata.mod["rna"].n_obs == output_mudata.mod["rna"].n_obs
    ), "Number of observations changed"
    assert (
        input_mudata.mod["rna"].n_vars == output_mudata.mod["rna"].n_vars
    ), "Number of variables changed"
    assert (
        "scanvi_embedding" in output_mudata.mod["rna"].obsm.keys()
    ), "Latent representation not added"
    assert "scanvi_pred" in output_mudata.mod["rna"].obs.keys(), "Predictions not added"
    assert (
        "scanvi_probability" in output_mudata.mod["rna"].obs.keys()
    ), "Probabilities not added"

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])


def test_pretrained_scanvi(run_component, random_h5mu_path, create_scanvi_model):
    scanvi_model_file, input_file_scanvi = create_scanvi_model()
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file_scanvi,
            "--scanvi_reference_model",
            scanvi_model_file,
            "--reference_obs_label",
            "cell_ontology_class",
            "--reference_max_epochs",
            "10",
            "--query_max_epochs",
            "10",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file_scanvi)
    output_mudata = mu.read_h5mu(output_file)

    assert (
        input_mudata.mod["rna"].n_obs == output_mudata.mod["rna"].n_obs
    ), "Number of observations changed"
    assert (
        input_mudata.mod["rna"].n_vars == output_mudata.mod["rna"].n_vars
    ), "Number of variables changed"
    assert (
        "scanvi_embedding" in output_mudata.mod["rna"].obsm.keys()
    ), "Latent representation not added"
    assert "scanvi_pred" in output_mudata.mod["rna"].obs.keys(), "Predictions not added"
    assert (
        "scanvi_probability" in output_mudata.mod["rna"].obs.keys()
    ), "Probabilities not added"

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])


def test_raises(
    run_component, random_h5mu_path, create_scvi_model, create_scanvi_model
):
    scvi_model_file, input_file_scvi = create_scvi_model(input_file, reference_file)
    scanvi_model_file, input_file_scanvi = create_scanvi_model()
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file_scanvi,
                "--scanvi_reference_model",
                scanvi_model_file,
                "--scvi_reference_model",
                scvi_model_file,
                "--reference_obs_label",
                "cell_ontology_class",
                "--reference_max_epochs",
                "10",
                "--query_max_epochs",
                "10",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"ValueError: Make sure to provide either an '--scvi_reference_model' or a '--scanvi_reference_model', but not both.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
