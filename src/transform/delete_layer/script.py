import sys
from mudata import read_h5ad, write_h5ad
import shutil
from pathlib import Path

## VIASH START
from mudata import read_h5mu
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": ['log_normalized'],
    "missing_ok": False,
    "output_compression": "lzf"
}
meta = {
    "functionality_name": "delete_layer",
    "resources_dir": "resources_test"
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

from compress_h5mu import compress_h5mu

def main():
    input_file, output_file, mod_name = Path(par["input"]), Path(par["output"]), par['modality']

    logger.info('Reading input file %s, modality %s.', input_file, mod_name)
    mod = read_h5ad(input_file, mod=mod_name)
    for layer in par['layer']:
        if layer not in mod.layers:
            if par['missing_ok']:
                continue
            raise ValueError(f"Layer '{layer}' is not present in modality {mod_name}.")
        logger.info('Deleting layer %s from modality %s.', layer, mod_name)
        del mod.layers[layer]

    logger.info('Writing output to %s.', par['output'])
    output_file_uncompressed = output_file.with_name(output_file.stem + "_uncompressed.h5mu") \
        if par["output_compression"] else output_file
    shutil.copyfile(par['input'], output_file_uncompressed)
    write_h5ad(filename=output_file_uncompressed, mod=mod_name, data=mod)
    if par["output_compression"]:
        compress_h5mu(output_file_uncompressed, output_file, compression=par["output_compression"])
        output_file_uncompressed.unlink()

    logger.info('Finished.')

if __name__ == "__main__":
    main()