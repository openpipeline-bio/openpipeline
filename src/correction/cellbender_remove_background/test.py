import subprocess
from os import path
import mudata as mu

## VIASH START
meta = {
  'executable': 'target/docker/correction/cellbender_remove_background/cellbender_remove_background',
  'resources_dir': 'resources_test/pbmc_1k_protein_v3'
}
## VIASH END

file_input = meta["resources_dir"] + "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
file_output = "output.h5mu"

print("> Check whether cellbender works when it should be working")

# run cellbender
cmd_pars = [
    meta["executable"],
    "--input", file_input,
    "--output", file_output,
    "--epochs", "5",
    "--output_compression", "gzip"
]
# todo: if cuda is available, add --cuda
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5mu"), "No output was created."

data = mu.read_h5mu("output.h5mu")

# check whether gex was found
assert data.mod["rna"].var["feature_types"].unique() == [
    "Gene Expression"
], "Output X should only contain Gene Expression vars."

# check whether ab counts were found
assert "prot" in data.mod, 'Output should contain data.mod["rna"].'