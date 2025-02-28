import csv
import os
import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal


##VIASH START
par = {
    "input": "output_test/split_modalities/foo_types.csv",
    "mod_dir": "output_test/split_modalities/h5mu",
    "orig_input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3.h5mu",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END


print("Loading data", flush=True)
with open(par["input"], "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    data = list(reader)

input_mu = mu.read_h5mu(par["orig_input"])

num_mod = len(data) - 1
num_files = len(os.listdir(par["mod_dir"]))

# Check if the number of files is equal to the number of lines in the csv
assert num_mod == num_files, f"Expected {num_mod} files, but found {num_files}."
assert (
    input_mu.n_mod == num_mod
), f"Expected {num_mod} modalities in {par['orig_input']} got {input_mu.n_mod} modalities."

rna_mod = mu.read_h5mu(os.path.join(par["mod_dir"], data[1][1]))
prot_mod = mu.read_h5mu(os.path.join(par["mod_dir"], data[2][1]))

# Check if the files exist and if the modality name is in the file name
for i, row in enumerate(data):
    if i == 0:
        continue
    # Check if the files exist and if the modality name is in the file name
    assert row[0] in row[1], f"Expected {row[0]} to be in {row[1]}."
    mod_fp = os.path.join(par["mod_dir"], row[1])
    assert os.path.exists(mod_fp), f"Expected {row[1]} to exist."
    # Check modality is correct in the h5mu file
    mod_mu = mu.read_h5mu(mod_fp)
    assert mod_mu.n_mod == 1, f"Expected 1 modality in {row[1]}."
    assert row[0] in mod_mu.mod.keys(), f"Expected {row[0]} to be the mod in {row[1]}."
    assert (
        row[0] in input_mu.mod.keys()
    ), f"Expected {row[0]} to be a mod in {par['orig_input']}."

# Check if extracted modalities are equal to the original modalities
assert_annotation_objects_equal(rna_mod.mod["rna"], input_mu.mod["rna"])
assert_annotation_objects_equal(prot_mod.mod["prot"], input_mu.mod["prot"])

print("Test successful!", flush=True)
