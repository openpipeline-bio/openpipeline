import subprocess
from os import path
import mudata as mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5mu"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--obsm_name", "fooleiden",
    "--resolution", "1",
    "--resolution", "0.25",
    "--output", output,
    "--output_compression", "gzip"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5mu"), "No output was created."

# read it with scanpy
data = mu.read_h5mu("output.h5mu")

# check whether leiden.custom.resolution was found
assert "1.0" in data.mod["rna"].obsm["fooleiden"].columns, 'Output should contain resolution 1.0.'
assert "0.25" in data.mod["rna"].obsm["fooleiden"].columns, 'Output should contain resolution 0.25.'
