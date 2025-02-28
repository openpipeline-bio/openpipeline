import sys
import os
import pytest
import subprocess
import re
import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
cl_nlp_emb_file = (
    f"{meta['resources_dir']}/annotation_test_data/ontology/cl.ontology.nlp.emb"
)
cl_ontology_file = f"{meta['resources_dir']}/annotation_test_data/ontology/cl.ontology"
cl_obo_file = f"{meta['resources_dir']}/annotation_test_data/ontology/cl.obo"
model_file = (
    f"{meta['resources_dir']}/annotation_test_data/onclass_model/example_file_model"
)


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--reference",
            reference_file,
            "--reference_obs_target",
            "cell_ontology_class",
            "--cl_nlp_emb_file",
            cl_nlp_emb_file,
            "--cl_ontology_file",
            cl_ontology_file,
            "--cl_obo_file",
            cl_obo_file,
            "--max_iter",
            "10",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == ["onclass_pred", "onclass_prob"]

    obs_values = output_mudata.mod["rna"].obs["onclass_prob"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at cell_ontology_class_prob has values outside the range [0, 1]"


def test_custom_obs(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--reference",
            reference_file,
            "--reference_obs_target",
            "cell_ontology_class",
            "--output_obs_predictions",
            "dummy_pred_1",
            "--output_obs_probability",
            "dummy_prob_1",
            "--cl_nlp_emb_file",
            cl_nlp_emb_file,
            "--cl_ontology_file",
            cl_ontology_file,
            "--cl_obo_file",
            cl_obo_file,
            "--max_iter",
            "10",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert set(output_mudata.mod["rna"].obs.keys()) == {"dummy_pred_1", "dummy_prob_1"}

    obs_values = output_mudata.mod["rna"].obs["dummy_prob_1"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at dummy_prob_1 has values outside the range [0, 1]"


def test_no_model_no_reference_error(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--input_var_gene_names",
                "gene_symbol",
                "--output",
                output_file,
                "--cl_nlp_emb_file",
                cl_nlp_emb_file,
                "--cl_ontology_file",
                cl_ontology_file,
                "--cl_obo_file",
                cl_obo_file,
                "--reference_obs_target",
                "cell_ontology_class",
            ]
        )
    assert re.search(
        r"ValueError: Make sure to provide either 'model' or 'reference', but not both.",
        err.value.stdout.decode("utf-8"),
    )


def test_pretrained_model(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--cl_nlp_emb_file",
            cl_nlp_emb_file,
            "--cl_ontology_file",
            cl_ontology_file,
            "--cl_obo_file",
            cl_obo_file,
            "--reference_obs_target",
            "cell_ontology_class",
            "--model",
            model_file,
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == ["onclass_pred", "onclass_prob"]

    obs_values = output_mudata.mod["rna"].obs["onclass_prob"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at cell_ontology_class_prob has values outside the range [0, 1]"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
