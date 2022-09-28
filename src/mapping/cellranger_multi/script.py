from __future__ import annotations

import os
import re
import subprocess
import tempfile
import logging
from sys import stdout
import pandas as pd
from typing import Optional, Any, Union
import tarfile
import pathlib

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
# The following code has been auto-generated by Viash.
par = {
  'output': '/path/to/output',
  'input': ['resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz'],
  'gex_reference': 'resources_test/reference_gencodev41_chr1//reference_cellranger.tar.gz',
  'vdj_reference': 'resources_test/10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz',
  'feature_reference': 'resources_test/10x_5k_anticmv/raw/feature_reference.csv',
  'library_id': ['5k_human_antiCMV_T_TBNK_connect_GEX_1_subset',
                 '5k_human_antiCMV_T_TBNK_connect_AB_subset',
                 '5k_human_antiCMV_T_TBNK_connect_VDJ_subset'],
  'library_type': ['Gene Expression', 'Antibody Capture', 'VDJ'],
  'library_lanes': None,
  'library_subsample': None,
  'gex_expect_cells': None,
  'gex_chemistry': 'auto',
  'gex_secondary_analysis': True,
  'gex_generate_bam': True,
  'gex_include_introns': True,
  'cell_multiplex_sample_id': None,
  'cell_multiplex_oligo_ids': None,
  'cell_multiplex_description': None,
  'dryrun': False
}
meta = {
  'n_proc': None,
  'memory_b': None,
  'memory_kb': None,
  'memory_mb': None,
  'memory_gb': None,
  'memory_tb': None,
  'memory_pb': None,
  'temp_dir': '/tmp'
}
## VIASH END

fastq_regex = r'([A-Za-z0-9\-_\.]+)_S(\d+)_L(\d+)_R(\d+)_(\d+)\.fastq\.gz'
# assert re.match(fastq_regex, "5k_human_GEX_1_subset_S1_L001_R1_001.fastq.gz") is not None


REFERENCES = ("gex_reference", "feature_reference", "vdj_reference")
REFERENCE_CONFIG_KEYS = {
    "gex_reference": "gene-expression",
    "feature_reference": "feature",
    "vdj_reference": "vdj"
}

LIBRARY_PARAMS = ("library_id", "library_type", "library_subsample", "library_lanes")
LIBRARY_CONFIG_KEYS = {'library_id': 'fastq_id',
                       'library_type': 'feature_types',
                       'library_subsample': 'subsample_rate',
                       'library_lanes': 'lanes'}


SAMPLE_PARAMS = ("cell_multiplex_sample_id", "cell_multiplex_oligo_ids", "cell_multiplex_description")
SAMPLE_PARAMS_CONFIG_KEYS = {'cell_multiplex_sample_id': 'sample_id', 
                             'cell_multiplex_oligo_ids': 'cmo_ids',
                             'cell_multiplex_description': 'description'}


def lengths_gt1(dic: dict[str, Optional[list[Any]]]) -> dict[str, int]:
    return {key: len(li) for key, li in dic.items() 
            if li is not None and len(li) > 1}
  
def strip_margin(text: str) -> str:
    return re.sub('(\n?)[ \t]*\|', '\\1', text)


def subset_dict(dictionary: dict[str, str], 
                keys: Union[dict[str, str], list[str]]) -> dict[str, str]:
  if isinstance(keys, (list, tuple)):
    keys = {key: key for key in keys}
  return {dest_key: dictionary[orig_key] 
          for orig_key, dest_key in keys.items() 
          if dictionary[orig_key] is not None}

def check_subset_dict_equal_length(group_name: str, 
                                   dictionary: dict[str, list[str]]) -> None:
    lens = lengths_gt1(dictionary)
    assert len(set(lens.values())) <= 1, f"The number of values passed to {group_name} "\
                                         f"arguments must be 0, 1 or all the same. Offenders: {lens}"

def process_params(par: dict[str, Any]) -> str:
    # if par_input is a directory, look for fastq files
    if len(par["input"]) == 1 and os.path.isdir(par["input"][0]):
        logger.info("Detected '--input' as a directory, "
                    "traversing to see if we can detect any FASTQ files.")
        par["input"] = [os.path.join(dp, f) 
                        for dp, _, filenames in os.walk(par["input"][0])
                        for f in filenames if re.match(fastq_regex, f) ]

    # check input fastq files
    for input in par["input"]:
        assert re.match(fastq_regex, os.path.basename(input)) is not None, \
               f"File name of --input '{input}' should match regex {fastq_regex}."
    
    # check lengths of libraries metadata 
    library_dict = subset_dict(par, LIBRARY_PARAMS)
    check_subset_dict_equal_length("Library", library_dict)
    # storing for later use
    par["libraries"] = library_dict

    cmo_dict = subset_dict(par, SAMPLE_PARAMS)
    check_subset_dict_equal_length("Cell multiplexing", cmo_dict)
    # storing for later use
    par["cmo"] = cmo_dict

    # use absolute paths
    par["input"] = [ os.path.abspath(f) for f in par["input"] ]
    for file_path in REFERENCES + ('output', ):
        if par[file_path]:
            logger.info('Making path %s absolute', par[file_path])
            par[file_path] = os.path.abspath(par[file_path])
    return par

def generate_dict_category(name: str, args: dict[str, str]) -> list[str]:
    title = [ f'[{name}]' ]
    values = [ f'{key},{val}' for key, val in args.items() if val is not None ]
    if len(values) > 0:
        return title + values + [""]
    else:
        return []

def generate_csv_category(name: str, args: dict[str, str]) -> list[str]:
    title = [ f'[{name}]' ]
    if len(args) > 0:
        values = [ pd.DataFrame(args).to_csv(index=False) ]
        return title + values + [""]
    else:
        return []

def generate_config(par: dict[str, Any], fastq_dir: str) -> str:
    serialized_refs = []
    for reference in REFERENCES:
        pars = subset_dict(par, {reference: 'ref'})
        serialized = generate_dict_category(REFERENCE_CONFIG_KEYS[reference], pars)
        serialized_refs.extend(serialized)

    # process libraries parameters
    library_pars = subset_dict(par, LIBRARY_CONFIG_KEYS)
    library_pars['fastqs'] = fastq_dir
    libraries_strs = generate_csv_category("libraries", library_pars)

    # process samples parameters
    cmo_pars = subset_dict(par, SAMPLE_PARAMS_CONFIG_KEYS)
    cmo_strs = generate_csv_category("samples", cmo_pars)
    
    # combine content
    content_list = serialized_refs + libraries_strs + cmo_strs
    return '\n'.join(content_list)

def main(par: dict[str, Any], meta: dict[str, Any]):
    logger.info("  Processing params")
    par = process_params(par)
    logger.info(par)

    # TODO: throw error or else Cell Ranger will
    # # Create output dir if not exists
    with tempfile.TemporaryDirectory(prefix="cellranger_multi-",
                                     dir=meta["temp_dir"]) as temp_dir:
        for reference_par_name in REFERENCES:
            reference = par[reference_par_name]
            logger.info('Looking at %s to check if it needs decompressing', reference)
            if tarfile.is_tarfile(reference):
                extaction_dir_name, _ = os.path.splitext(reference)
                extaction_dir_name, _ = os.path.splitext(extaction_dir_name)

                unpacked_directory = os.path.join(temp_dir, os.path.basename(extaction_dir_name))
                logger.info('Extracting %s to %s', reference, unpacked_directory)

                with tarfile.open(reference, 'r') as open_tar:
                    rootDirs = [ rootDir for rootDir in open_tar.getnames() if '/' not in rootDir ]
                    members = open_tar.getmembers()
                    # if there is only one rootDir (and there are files in that directory)
                    # strip that directory name from the destination folder
                    if len(rootDirs) == 1 and len(members) > 1 and rootDirs[0] != '.':
                        for mem in members:
                            mem.path = '/'.join(pathlib.Path(mem.path).parts[1:])
                    open_tar.extractall(unpacked_directory, members=[mem for mem in members if len(mem.path) > 0])
                par[reference_par_name] = unpacked_directory

        # Creating symlinks of fastq files to tempdir
        input_symlinks_dir =  os.path.join(temp_dir, "input_symlinks")
        os.mkdir(input_symlinks_dir)
        for fastq in par['input']:
            destination = os.path.join(input_symlinks_dir, os.path.basename(fastq))
            os.symlink(fastq, destination)

        logger.info("  Creating config file")
        config_content = generate_config(par, input_symlinks_dir)

        logger.info("  Creating Cell Ranger argument")
        temp_id="run"
        proc_pars=["--disable-ui", "--id", temp_id]

        if meta["n_proc"]:
            proc_pars.append(f"--localcores={meta['n_proc']}")

        if meta["memory_gb"]:
            proc_pars.append(f"--localmem={int(meta['memory_gb']) - 2}")

        ## Run pipeline
        if par["dryrun"]:
            cmd = ["cellranger multi"] + proc_pars + ["--csv=config.csv"]
            logger.info("> " + ' '.join(cmd))
            logger.info("Contents of 'config.csv':")
            logger.info(config_content)
        else:
            # write config file
            config_file = os.path.join(temp_dir, "config.csv")
            with open(config_file, "w") as f:
                f.write(config_content)
            proc_pars.append(f"--csv={config_file}")

        # run process
        cmd = ["cellranger", "multi"] + proc_pars
        logger.info("> " + ' '.join(cmd))
        _ = subprocess.check_call(
            cmd,
            cwd=temp_dir
        )

        # look for output dir file
        tmp_output_dir = os.path.join(temp_dir, temp_id, "outs")
        expected_files = {
            "multi": "directory", 
            "per_sample_outs": "directory", 
            "config.csv": "file",
        }
        for file, type_ in expected_files.items():
            path = os.path.join(tmp_output_dir, file)
            if not os.path.exists(path):
                raise ValueError(f"Could not find expected {type_} '{path}'")

if __name__ == "__main__":
    main(par, meta)

