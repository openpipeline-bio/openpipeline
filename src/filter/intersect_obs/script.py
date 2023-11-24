import mudata as mu
import anndata as ad
import sys
from pathlib import Path
import shutil

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modalities": ["rna", "prot"],
    "output": "output.h5mu"
}
meta = {
    
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

# START TEMPORARY WORKAROUND compress_h5mu
# reason: resources aren't available when using Nextflow fusion

# from compress_h5mu import compress_h5mu
from h5py import Group, Dataset
from h5py import File as H5File
from typing import Union, Literal
from functools import partial
from pathlib import Path

def compress_h5mu(input_path: Union[str, Path], 
                output_path: Union[str, Path], 
                compression: Union[Literal['gzip'], Literal['lzf']]):
    input_path, output_path = str(input_path), str(output_path)

    def copy_attributes(in_object, out_object):
        for key, value in in_object.attrs.items():
            out_object.attrs[key] = value

    def visit_path(output_h5: H5File,
                   compression: Union[Literal['gzip'], Literal['lzf']], 
                   name: str, object: Union[Group, Dataset]):
            if isinstance(object, Group):
                new_group = output_h5.create_group(name)
                copy_attributes(object, new_group)
            elif isinstance(object, Dataset):
                # Compression only works for non-scalar Dataset objects
                # Scalar objects dont have a shape defined
                if not object.compression and object.shape not in [None, ()]: 
                    new_dataset = output_h5.create_dataset(name, data=object, compression=compression)
                    copy_attributes(object, new_dataset)
                else:
                    output_h5.copy(object, name)
            else:
                raise NotImplementedError(f"Could not copy element {name}, "
                                          f"type has not been implemented yet: {type(object)}")

    with H5File(input_path, 'r') as input_h5, H5File(output_path, 'w', userblock_size=512) as output_h5:
        copy_attributes(input_h5, output_h5)
        input_h5.visititems(partial(visit_path, output_h5, compression))

    with open(input_path, "rb") as input_bytes:
        # Mudata puts metadata like this in the first 512 bytes:
        # MuData (format-version=0.1.0;creator=muon;creator-version=0.2.0)
        # See mudata/_core/io.py, read_h5mu() function
        starting_metadata = input_bytes.read(100)
        # The metadata is padded with extra null bytes up until 512 bytes
        truncate_location = starting_metadata.find(b"\x00")
        starting_metadata = starting_metadata[:truncate_location]
    with open(output_path, "br+") as f:
        nbytes = f.write(starting_metadata)
        f.write(b"\0" * (512 - nbytes))
# END TEMPORARY WORKAROUND compress_h5mu

logger = setup_logger()

def main():
    modality_names = par['modalities']

    if len(modality_names) < 2:
        raise ValueError("Please provide two more more modalities.")
    
    obs_names = {}
    for mod_name in par['modalities']:
        try:
            modality = mu.read_h5ad(filename=par['input'], mod=mod_name)
        except KeyError:
            raise ValueError(f"Modality {mod_name} does not exist for file {par['input']}.")

        obs_names[mod_name] = modality.obs_names.copy()
        del modality
    
    intersected_index = None
    for mod_name, mod_index in obs_names.items():
        if intersected_index is None:
            intersected_index = mod_index
            continue
        intersected_index = intersected_index.intersection(mod_index)
    

    output_file = Path(par['output'])
    output_file_uncompressed = output_file.with_name(output_file.stem + "_uncompressed.h5mu")
    output_file_uncompressed.touch()

    mdata = mu.MuData({modality: ad.AnnData() for modality in modality_names})
    mdata.write(output_file_uncompressed, compression=par['output_compression'])
    
    for mod_name in modality_names:
        modality = mu.read_h5ad(filename=par['input'], mod=mod_name)
        intersected_modality = modality[intersected_index]
        mu.write_h5ad(output_file_uncompressed, data=intersected_modality, mod=mod_name)

    if par['output_compression']:
        compress_h5mu(output_file_uncompressed, output_file, compression=par['output_compression'])
        output_file_uncompressed.unlink()
    else:
        shutil.move(output_file_uncompressed, output_file)
    
if __name__ == "__main__":
    main()