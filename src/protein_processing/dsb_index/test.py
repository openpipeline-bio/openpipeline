import subprocess
from os import path
import muon as mu
import logging
from sys import stdout
import pandas as pd

## VIASH START
meta = {
    'functionality_name': 'dsb_index',
    'executable': './target/dsb_index',
    'resources_dir': 'resources_test/'
    
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5"
cell_index = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix/barcodes.tsv.gz"
output_cell_idx = "dsb_output/cell_idx.csv"
output_empty_idx = "dsb_output/empty_idx.csv"
cmd_pars = [
    meta["executable"],
    "--data_raw", input,
    "--output", "dsb_output",
    "--cell_index", cell_index,
    "--empty_counts_range", "1.5:2.8"
]
try:
    subprocess.check_output(cmd_pars)
except subprocess.CalledProcessError as e:
    print(e.stdout.decode("utf-8"))
    raise e

logger.info("> Check if output was created.")
assert path.exists(output_cell_idx), "No output for cell index was created."
assert path.exists(output_empty_idx), "No output for empty index was created."


logger.info("> Check whether output cell index has the samw shape as the cell_index input.")
cell_index = pd.read_csv(cell_index, header=None).iloc[:, 0].tolist()
output_cell_index = pd.read_csv(output_cell_idx, header=None).iloc[:, 0].tolist()
assert len(output_cell_index) == len(cell_index), 'Output cell index has the samw shape as the cell_index input.'

