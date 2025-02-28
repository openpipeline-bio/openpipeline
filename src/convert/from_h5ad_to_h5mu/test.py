import sys
import pytest
import mudata as mu
from openpipelinetest_utils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/executable/convert/from_h5ad_to_h5mu/from_h5ad_to_h5mu",
    "config": "./src/convert/from_h5ad_to_h5mu/config.vsh.yaml",
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def test_run(run_component, random_h5mu_path, random_path):
    mdata = mu.read_h5mu(input)
    tmp_rna = random_path(extension="h5ad")
    tmp_prot = random_path(extension="h5ad")
    mdata.mod["rna"].write_h5ad(tmp_rna)
    mdata.mod["prot"].write_h5ad(tmp_prot)

    tmp_output = random_h5mu_path()

    cmd_pars = [
        "--modality",
        "rna",
        "--input",
        tmp_rna,
        "--modality",
        "prot",
        "--input",
        tmp_prot,
        "--output",
        tmp_output,
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert tmp_output.is_file(), "No output was created."

    mdata2 = mu.read_h5mu(tmp_output)

    assert list(mdata2.mod.keys()) == ["rna", "prot"]

    assert_annotation_objects_equal(mdata2.mod["rna"], tmp_rna)
    assert_annotation_objects_equal(mdata2.mod["prot"], tmp_prot)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
