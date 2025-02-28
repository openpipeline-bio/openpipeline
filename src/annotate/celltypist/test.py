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

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
model_file = (
    f"{meta['resources_dir']}/annotation_test_data/celltypist_model_Immune_All_Low.pkl"
)
celltypist_input_file = (
    f"{meta['resources_dir']}/annotation_test_data/demo_2000_cells.h5mu"
)


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_layer",
            "log_normalized",
            "--reference",
            reference_file,
            "--reference_layer",
            "log_normalized",
            "--reference_obs_target",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert {"celltypist_pred", "celltypist_probability"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs"

    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at celltypist_probability has values outside the range [0, 1]"


def test_set_params(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_layer",
            "log_normalized",
            "--reference",
            reference_file,
            "--reference_layer",
            "log_normalized",
            "--reference_obs_target",
            "cell_ontology_class",
            "--reference_var_gene_names",
            "ensemblid",
            "--feature_selection",
            "True",
            "--majority_voting",
            "True",
            "--C",
            "0.5",
            "--max_iter",
            "100",
            "--use_SGD",
            "--min_prop",
            "0.1",
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert {"celltypist_pred", "celltypist_probability"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs"

    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at celltypist_probability has values outside the range [0, 1]"


def test_with_model(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            celltypist_input_file,
            "--model",
            model_file,
            "--reference_obs_targets",
            "cell_type",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    output_mudata = mu.read_h5mu(output_file)

    assert {"celltypist_pred", "celltypist_probability"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs"

    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(
        0 <= value <= 1 for value in obs_values
    ), ".obs at celltypist_probability has values outside the range [0, 1]"


def test_fail_invalid_input_expression(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    # fails because input data are not lognormalized
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--reference",
                reference_file,
                "--reference_var_gene_names",
                "ensemblid",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell",
        err.value.stdout.decode("utf-8"),
    )

    # fails because reference data are not lognormalized
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--layer",
                "log_normalized",
                "--reference",
                reference_file,
                "--reference_var_gene_names",
                "ensemblid",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
