from pathlib import Path
import mudata as mu
import numpy as np
import logging
import sys
import pytest
import uuid

## VIASH START
meta = {
    'functionality_name': './target/native/filter/do_filter/do_filter',
    'resources_dir': 'resources_test/',
     'executable': './target/docker/filter/do_filter/do_filter',
     'config': './src/filter/do_filter/config.vsh.yaml'
}
## VIASH END
resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_path = f"{resources_dir}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

@pytest.fixture
def random_h5mu_path(tmp_path):
    unique_filename = f"{str(uuid.uuid4())}.h5mu"
    temp_file = tmp_path / unique_filename
    return temp_file

@pytest.fixture
def write_to_temp_file(tmp_path):
    def write_h5mu_wrapper(mudata_object):
        unique_filename = f"{str(uuid.uuid4())}.h5mu"
        temp_file = tmp_path / unique_filename
        mudata_object.write(temp_file)
        return temp_file
    return write_h5mu_wrapper

@pytest.fixture()
def input_data():
    mu_in = mu.read_h5mu(input_path)
    return mu_in

@pytest.fixture()
def original_n_obs(input_data):
    return input_data.mod['rna'].n_obs

@pytest.fixture()
def original_n_vars(input_data):
    return input_data.mod['rna'].n_vars

@pytest.fixture()
def test_data_filter_nothing(input_data, write_to_temp_file):
    rna_mod = input_data.mod['rna']
    rna_mod.obs["filter_none"] = np.repeat(True, rna_mod.n_obs)
    return write_to_temp_file(input_data)

@pytest.fixture()
def test_data_filter_with_random(input_data, write_to_temp_file):
    rna_mod = input_data.mod['rna']
    rna_mod.obs["filter_with_random"] = np.random.choice([False, True], size=rna_mod.n_obs)
    rna_mod.var["filter_with_random"] = np.random.choice([False, True], size=rna_mod.n_vars)
    return write_to_temp_file(input_data)

def test_filtering_a_little_bit(run_component,
                                test_data_filter_with_random,
                                random_h5mu_path,
                                original_n_obs,
                                original_n_vars):
    component_output = run_component([
        "--input", test_data_filter_with_random,
        "--output", random_h5mu_path,
        "--obs_filter", "filter_with_random",
        "--var_filter", "filter_with_random"
        ])
    assert random_h5mu_path.is_file(), "Output file not found"
    mu_out = mu.read_h5mu(random_h5mu_path)
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs < original_n_obs, "Some RNA obs should have been filtered"
    assert new_vars < original_n_vars,"Some RNA vars should have been filtered"
    assert b"Filtering modality 'rna' observations by .obs['filter_with_random']" in component_output
    assert b"Filtering modality 'rna' variables by .var['filter_with_random']" in component_output

def test_filter_nothing(run_component,
                        test_data_filter_nothing,
                        random_h5mu_path,
                        original_n_obs,
                        original_n_vars):
    run_component([
        "--input", test_data_filter_nothing,
        "--output", random_h5mu_path,
        "--obs_filter", "filter_none"])
    assert random_h5mu_path.is_file(), "Output file not found"
    mu_out = mu.read_h5mu(random_h5mu_path)
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs == original_n_obs, "No RNA obs should have been filtered"
    assert new_vars == original_n_vars, "No RNA vars should have been filtered"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))
