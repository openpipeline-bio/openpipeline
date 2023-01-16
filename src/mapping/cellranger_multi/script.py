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
from pathlib import Path
import shutil
from itertools import chain

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
# The following code has been auto-generated by Viash.
par = {
  'output': './cellranger_test_output',
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
  'gex_secondary_analysis': False,
  'gex_generate_bam': False,
  'gex_include_introns': False,
  'cell_multiplex_sample_id': None,
  'cell_multiplex_oligo_ids': None,
  'cell_multiplex_description': None,
  'dryrun': False
}
meta = {
  'cpus': 10,
  'memory_b': None,
  'memory_kb': None,
  'memory_mb': None,
  'memory_gb': 15,
  'memory_tb': None,
  'memory_pb': None,
  'temp_dir': '/tmp'
}
## VIASH END

fastq_regex = r'([A-Za-z0-9\-_\.]+)_S(\d+)_L(\d+)_[RI](\d+)_(\d+)\.fastq\.gz'
# assert re.match(fastq_regex, "5k_human_GEX_1_subset_S1_L001_R1_001.fastq.gz") is not None

# Invert some parameters. Keep the original ones in the config for compatibility
inverted_params = {
    "gex_generate_no_bam": "gex_generate_bam",
    "gex_no_secondary_analysis": "gex_secondary_analysis"
}
for inverted_param, param in inverted_params.items():
    par[inverted_param] = not par[param] if par[param] is not None else None
    del par[param]

GEX_CONFIG_KEYS = {
    "gex_reference": "reference",
    "gex_expect_cells": "expect-cells",
    "gex_chemistry": "chemistry",
    "gex_no_secondary_analysis": "no-secondary",
    "gex_generate_no_bam": "no-bam",
    "gex_include_introns": "include-introns"
}
FEATURE_CONFIG_KEYS = {"feature_reference": "reference"}
VDJ_CONFIG_KEYS = {"vdj_reference": "reference"}

REFERENCE_SECTIONS = {
    "gene-expression": (GEX_CONFIG_KEYS, "index"),
    "feature": (FEATURE_CONFIG_KEYS, "index"),
    "vdj": (VDJ_CONFIG_KEYS, "index")
}

LIBRARY_CONFIG_KEYS = {'library_id': 'fastq_id',
                       'library_type': 'feature_types',
                       'library_subsample': 'subsample_rate',
                       'library_lanes': 'lanes'}
SAMPLE_PARAMS_CONFIG_KEYS = {'cell_multiplex_sample_id': 'sample_id',
                             'cell_multiplex_oligo_ids': 'cmo_ids',
                             'cell_multiplex_description': 'description'}


# These are derived from the dictionaries above
REFERENCES = tuple(reference_param for reference_param, cellranger_param
                   in chain(GEX_CONFIG_KEYS.items(), FEATURE_CONFIG_KEYS.items(), VDJ_CONFIG_KEYS.items())
                   if cellranger_param == "reference")
LIBRARY_PARAMS = tuple(LIBRARY_CONFIG_KEYS.keys())
SAMPLE_PARAMS = tuple(SAMPLE_PARAMS_CONFIG_KEYS.keys())


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
    par["input"] = [Path(fastq) for fastq in par["input"]]
    if len(par["input"]) == 1 and par["input"][0].is_dir():
        logger.info("Detected '--input' as a directory, "
                    "traversing to see if we can detect any FASTQ files.")
        par["input"] = [input_path for input_path in par["input"].rglob('*')
                        if re.match(fastq_regex, input_path.name) ]

    # check input fastq files
    for input_path in par["input"]:
        assert re.match(fastq_regex, input_path.name) is not None, \
               f"File name of --input '{input_path}' should match regex {fastq_regex}."

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
    par["input"] = [input_path.resolve() for input_path in par["input"]]
    for file_path in REFERENCES + ('output', ):
        if par[file_path]:
            logger.info('Making path %s absolute', par[file_path])
            par[file_path] = Path(par[file_path]).resolve()
    return par

def generate_dict_category(name: str, args: dict[str, str]) -> list[str]:
    title = [ f'[{name}]' ]
    values = [ f'{key},{val}' for key, val in args.items() if val is not None ]
    if len(values) > 0:
        return title + values + [""]
    else:
        return []

def generate_csv_category(name: str, args: dict[str, str], orient: str) -> list[str]:
    assert orient in ("index", "columns")
    if not args:
        return []
    title = [ f'[{name}]' ]
    # Which index to include in csv section is based on orientation
    to_csv_args = {"index": (orient=="index"), "header": (orient=="columns")}
    values = [pd.DataFrame.from_dict(args, orient=orient).to_csv(**to_csv_args).strip()]
    return title + values + [""]


def generate_config(par: dict[str, Any], fastq_dir: str) -> str:
    content_list = []
    par["fastqs"] = fastq_dir
    libraries = dict(LIBRARY_CONFIG_KEYS, **{"fastqs": "fastqs"})
    #TODO: use the union (|) operator when python is updated to 3.9
    all_sections = dict(REFERENCE_SECTIONS, 
                        **{"libraries": (libraries, "columns")},
                        **{"samples": (SAMPLE_PARAMS_CONFIG_KEYS, "columns")})
    for section_name, (section_params, orientation) in all_sections.items():
        reference_pars = subset_dict(par, section_params)
        content_list += generate_csv_category(section_name, reference_pars, orient=orientation)

    return '\n'.join(content_list)

def main(par: dict[str, Any], meta: dict[str, Any]):
    logger.info("  Processing params")
    par = process_params(par)
    logger.info(par)

    # TODO: throw error or else Cell Ranger will
    with tempfile.TemporaryDirectory(prefix="cellranger_multi-",
                                     dir=meta["temp_dir"]) as temp_dir:
        temp_dir_path = Path(temp_dir)
        for reference_par_name in REFERENCES:
            reference = par[reference_par_name]
            logger.info('Looking at %s to check if it needs decompressing', reference)
            if Path(reference).is_file() and tarfile.is_tarfile(reference):
                extaction_dir_name = Path(reference.stem).stem # Remove two extensions (if they exist)
                unpacked_directory = temp_dir_path / extaction_dir_name
                logger.info('Extracting %s to %s', reference, unpacked_directory)

                with tarfile.open(reference, 'r') as open_tar:
                    members = open_tar.getmembers()
                    root_dirs = [member for member in members if member.isdir()
                                 and member.name != '.' and '/' not in member.name]
                    # if there is only one root_dir (and there are files in that directory)
                    # strip that directory name from the destination folder
                    if len(root_dirs) == 1:
                        for mem in members:
                            mem.path = Path(*Path(mem.path).parts[1:])
                    members_to_move = [mem for mem in members if mem.path != Path('.')]
                    open_tar.extractall(unpacked_directory, members=members_to_move)
                par[reference_par_name] = unpacked_directory

        # Creating symlinks of fastq files to tempdir
        input_symlinks_dir = temp_dir_path / "input_symlinks"
        input_symlinks_dir.mkdir()
        for fastq in par['input']:
            destination = input_symlinks_dir / fastq.name
            destination.symlink_to(fastq)

        logger.info("  Creating config file")
        config_content = generate_config(par, input_symlinks_dir)

        logger.info("  Creating Cell Ranger argument")
        temp_id="run"
        proc_pars=["--disable-ui", "--id", temp_id]

        command_line_parameters = {
            "--localcores": meta['cpus'],
            "--localmem": int(meta['memory_gb']) - 2 if meta['memory_gb'] else None,
        }
        for param, param_value in command_line_parameters.items():
            if param_value:
                proc_pars.append(f"{param}={param_value}")

        ## Run pipeline
        if par["dryrun"]:
            par['output'].mkdir(parents=True, exist_ok=True)

            # write config file
            config_file = par['output'] / "config.csv"
            with open(config_file, "w") as f:
                f.write(config_content)
            proc_pars.append(f"--csv={config_file}")

            # display command that would've been used
            cmd = ["cellranger multi"] + proc_pars + ["--csv=config.csv"]
            logger.info("> " + ' '.join(cmd))
        else:
            # write config file to execution directory
            config_file = temp_dir_path / "config.csv"
            with open(config_file, "w") as f:
                f.write(config_content)
            proc_pars.append(f"--csv={config_file}")

            # Already copy config file to output directory
            par['output'].mkdir(parents=True, exist_ok=True)
            with (par['output'] / "config.csv").open('w') as open_config:
                open_config.write(config_content)

            # run process
            cmd = ["cellranger", "multi"] + proc_pars
            logger.info("> " + ' '.join(cmd))
            try:
                process_output = subprocess.run(
                    cmd,
                    cwd=temp_dir,
                    check=True,
                    capture_output=True
                )
            except subprocess.CalledProcessError as e:
                print(e.output.decode('utf-8'), flush=True)
                raise e
            else:
                # Write stdout output to output folder
                with (par["output"] / "cellranger_multi.log").open('w') as open_log:
                    open_log.write(process_output.stdout.decode('utf-8'))
                print(process_output.stdout.decode('utf-8'), flush=True)

            # look for output dir file
            tmp_output_dir = temp_dir_path / temp_id / "outs"
            expected_files = {
                Path("multi"): Path.is_dir,
                Path("per_sample_outs"): Path.is_dir,
                Path("config.csv"): Path.is_file,
            }
            for file_path, type_func in expected_files.items():
                output_path = tmp_output_dir / file_path
                if not type_func(output_path):
                    raise ValueError(f"Could not find expected '{output_path}'")

            for output_path in tmp_output_dir.rglob('*'):
                if output_path.name != "config.csv": # Already created
                    shutil.move(str(output_path), par['output'])

if __name__ == "__main__":
    main(par, meta)