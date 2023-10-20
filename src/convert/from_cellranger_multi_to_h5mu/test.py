import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    'executable': './target/docker/convert/from_cellranger_multi_to_h5mu/from_cellranger_multi_to_h5mu',
    'resources_dir': 'resources_test/',
    'config': 'src/convert/from_cellranger_multi_to_h5mu/config.vsh.yaml'
}
## VIASH END

input_anticmv = f"{meta['resources_dir']}/10x_5k_anticmv/processed/10x_5k_anticmv.cellranger_multi.output.output"
input_lung_crispr = f"{meta['resources_dir']}/10x_5k_lung_crispr/processed/10x_5k_lung_crispr.cellranger_multi.output.output"

def test_cellranger_multi_basic(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_anticmv,
        "--output", str(output_path),
        "--output_compression", "gzip"
    ])
    assert output_path.is_file()

    # check output
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ['rna', 'prot', 'vdj_t']
    assert list(converted_data.uns.keys()) == ['metrics_cellranger']
    expected_metrics = ['Category', 'Library Type', 'Grouped By', 'Group Name', 'Metric Name', 'Metric Value']
    assert converted_data.uns['metrics_cellranger'].columns.to_list() == expected_metrics
    
def test_cellranger_multi_to_h5mu_crispr(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_lung_crispr,
        "--output", str(output_path),
        "--output_compression", "gzip"])
    assert output_path.is_file()

    # check output
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ['rna', 'gdo']
    assert list(converted_data.uns.keys()) == ['metrics_cellranger']
    assert 'perturbation_efficiencies_by_feature' in converted_data.mod['gdo'].uns
    assert 'perturbation_efficiencies_by_target' in converted_data.mod['gdo'].uns
    assert 'feature_reference' not in converted_data.mod['rna'].uns
    assert 'feature_reference' in converted_data.mod['gdo'].uns

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))