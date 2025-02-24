import sys
from pathlib import Path
import os
import pytest

import mudata as md
import muon as mu
import numpy as np
import pandas as pd

## VIASH START
meta = {
    'executable': './target/docker/qc/calculate_atac_qc_metrics/calculate_atac_qc_metrics',
    'resources_dir': "./resources_test/cellranger_atac_tiny_bcl/counts/",
    'config': './src/qc/calculate_atac_qc_metrics/config.vsh.yaml',
    'cpus': 2
}
## VIASH END

@pytest.fixture
def tiny_atac_mudata(tmp_path):
    resources_dir = Path(meta["resources_dir"])
    mdata = mu.read_10x_h5(resources_dir / "counts" / "filtered_peak_bc_matrix.h5")
    mu.atac.tl.locate_fragments(mdata, fragments=str(resources_dir / "counts" / "fragments.tsv.gz"))
    assert "files" in mdata.mod["atac"].uns.keys()
    assert "fragments" in mdata.mod["atac"].uns["files"].keys()

    # Read features annotation and save it to uns
    peak_annotation = pd.read_csv(resources_dir / "counts" / "peak_annotation.tsv", sep="\t")
    peak_annotation["gene"] = peak_annotation["gene"].astype(str)  # Fixes saving error
    mu.atac.tl.add_peak_annotation(mdata.mod["atac"], peak_annotation)
    
    # Simulate clustering to not install leiden dependencies
    mdata.mod["atac"].obs["leiden"] = pd.Categorical(np.random.choice(np.arange(5), size=mdata.n_obs))

    mdata_path = tmp_path / "tiny_atac.h5mu"
    mdata.write(mdata_path)

    return mdata_path

@pytest.mark.parametrize("mudata", ["tiny_atac_mudata"])
def test_marker_peaks(run_component, request, mudata, tmp_path):
    input_path = request.getfixturevalue(mudata)
    output_path = tmp_path / "foo.h5mu"

    args = [
        "--input", str(input_path),
        "--output", str(output_path),
        "--output_marker_peaks", str(tmp_path / "marker_peaks.tsv"),
        "--modality", "atac",
        "--groupby", "leiden"
    ]

    run_component(args)
    assert output_path.is_file()
    assert "marker_peaks.tsv" in os.listdir(output_path.parent)

    data_with_markers = md.read(output_path)

    assert "rank_genes_groups" in data_with_markers.mod["atac"].uns

    marker_peaks = pd.read_csv(tmp_path / "marker_peaks.tsv", sep="\t")
    assert marker_peaks.shape[0] > 0
    assert "0_pvals_adj" in marker_peaks.columns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))