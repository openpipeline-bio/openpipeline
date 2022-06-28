import os
import re
import subprocess
import tempfile

## VIASH START
par = {
  'input': ['resources_test/bd_rhapsody_wta_test/raw/sample_R1_.fastq.gz', 'resources_test/bd_rhapsody_wta_test/raw/sample_R2_.fastq.gz'],
  'output': 'output_dir',
  'subsample': None,
  'reference_genome': 'resources_test/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz',
  'transcriptome_annotation': 'resources_test/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf',
  'exact_cell_count': None,
  'disable_putative_calling': False,
  'parallel': True,
  'timestamps': False,
  'abseq_reference': [],
  'supplemental_reference': []
}
meta = {
  'resources_dir': 'src/mapping/bd_rhapsody_wta',
  'temp_dir': os.getenv("VIASH_TEMP")
}
## VIASH END

def strip_margin(text):
  return re.sub('\n[ \t]*\|', '\n', text)

# if par_input is a directory, look for fastq files
if len(par["input"]) == 1 and os.path.isdir(par["input"][0]):
  par["input"] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk(par["input"]) for f in filenames if re.match(r'.*\.fastq.gz', f) ]

# use absolute paths
par["input"] = [ os.path.abspath(f) for f in par["input"] ]
par["reference_genome"] = os.path.abspath(par["reference_genome"])
par["transcriptome_annotation"] = os.path.abspath(par["transcriptome_annotation"])
par["output"] = os.path.abspath(par["output"])
if par["abseq_reference"]:
  par["abseq_reference"] = [ os.path.abspath(f) for f in par["abseq_reference"] ]
if par["supplemental_reference"]:
  par["supplemental_reference"] = [ os.path.abspath(f) for f in par["supplemental_reference"] ]

# create output dir if not exists
if not os.path.exists(par["output"]):
  os.makedirs(par["output"])

# Create params file
config_file = os.path.join(par["output"], "config.yml")
endl = "\n"

content_list = [f"""#!/usr/bin/env cwl-runner

cwl:tool: rhapsody

# This is a YML file used to specify the inputs for a BD Genomics WTA Rhapsody Analysis pipeline run. See the
# BD Genomics Analysis Setup User Guide (Doc ID: 47383) for more details.

## Reads (required) - Path to your read files in the FASTQ.GZ format. You may specify as many R1/R2 read pairs as you want.
Reads:
"""]

for file in par["input"]:
  content_list.append(strip_margin(f"""\
    | - class: File
    |   location: "{file}"
    |"""))

content_list.append(strip_margin(f"""\
  |
  |## Reference_Genome (required) - Path to STAR index for tar.gz format. See Doc ID: 47383 for instructions to obtain pre-built STAR index file.
  |Reference_Genome:
  |   class: File
  |   location: "{par["reference_genome"]}"
  |
  |## Transcriptome_Annotation (required) - Path to GTF annotation file
  |Transcriptome_Annotation:
  |   class: File
  |   location: "{par["transcriptome_annotation"]}"
  |"""))

if par["abseq_reference"]:
  content_list.append(strip_margin(f"""\
    |
    |## AbSeq_Reference (optional) - Path to the AbSeq reference file in FASTA format.  Only needed if BD AbSeq Ab-Oligos are used.
    |## For putative cell calling using an AbSeq dataset, please provide an AbSeq reference fasta file as the AbSeq_Reference.
    |AbSeq_Reference:
    |"""))
  for file in par["abseq_reference"]:
    content_list.append(strip_margin(f"""\
      | - class: File
      |   location: {file}
      |"""))

if par["supplemental_reference"]:
  content_list.append(strip_margin(f"""\
    |
    |## Supplemental_Reference (optional) - Path to the supplemental reference file in FASTA format.  Only needed if there are additional transgene sequences used in the experiment.
    |Supplemental_Reference:
    |"""))
  for file in par["supplemental_reference"]:
    content_list.append(strip_margin(f"""\
      | - class: File
      |   location: {file}
      |"""))

## Putative Cell Calling Settings
content_list.append(strip_margin(f"""\
  |
  |####################################
  |## Putative Cell Calling Settings ##
  |####################################
  |"""))

if par["exact_cell_count"]:
  content_list.append(strip_margin(f"""\
    |## Exact cell count (optional) - Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count
    |Exact_Cell_Count: {par["exact_cell_count"]}
    |"""))

if par["disable_putative_calling"]:
  content_list.append(strip_margin(f"""\
    |## Disable Refined Putative Cell Calling (optional) - Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.
    |## The values can be true or false. By default, the refined algorithm is used.
    |Basic_Algo_Only: {str(par["disable_putative_calling"]).lower()}
    |"""))

## Subsample Settings
content_list.append(strip_margin(f"""\
  |
  |########################
  |## Subsample Settings ##
  |########################
  |"""
))

if par["subsample"]:
  content_list.append(strip_margin(f"""\
    |## Subsample (optional) - A number >1 or fraction (0 < n < 1) to indicate the number or percentage of reads to subsample.
    |Subsample: {par["subsample"]}
    |"""))

if par["subsample_seed"]:
  content_list.append(strip_margin(f"""\
    |## Subsample seed (optional) - A seed for replicating a previous subsampled run.
    |Subsample_seed: {par["subsample_seed"]}
    |"""))


## Multiplex options
content_list.append(strip_margin(f"""\
  |
  |#######################
  |## Multiplex options ##
  |#######################
  |"""
))

if par["sample_tags_version"]:
  content_list.append(strip_margin(f"""\
    |## Sample Tags Version (optional) - Specify if multiplexed run: human, hs, mouse or mm
    |Sample_Tags_Version: {par["sample_tags_version"]}
    |"""))

if par["tag_names"]:
  content_list.append(strip_margin(f"""\
    |## Tag_Names (optional) - Specify the tag number followed by '-' and the desired sample name to appear in Sample_Tag_Metrics.csv
    |# Do not use the special characters: &, (), [], {{}},  <>, ?, |
    |Tag_Names: [{', '.join(par["tag_names"])}]
    |"""))

## Write config to file
config_content = ''.join(content_list)

with open(config_file, "w") as f:
  f.write(config_content)

## Process parameters
proc_pars = ["--no-container"]
# proc_pars = []

if par["parallel"]:
  proc_pars.append("--parallel")

if par["timestamps"]:
  proc_pars.append("--timestamps")

## Run pipeline
cwl_file=os.path.abspath(os.path.join(meta["resources_dir"], "rhapsody_wta_1.9.1_nodocker.cwl"))

# sed -i 's#"ramMin": [^,]*,#"ramMin": 2000,#' "$meta_resources_dir/rhapsody_wta_1.9.1_nodocker.cwl"
# sed -i 's#"coresMin": [^,]*,#"coresMin": 1,#' "$meta_resources_dir/rhapsody_wta_1.9.1_nodocker.cwl"

with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]) as temp_dir:
  cmd = ["cwl-runner"] + proc_pars + [cwl_file, os.path.basename(config_file)]

  env = dict(os.environ)
  env["TMPDIR"] = temp_dir

  print("> " + ' '.join(cmd))

  p = subprocess.Popen(
    cmd,
    cwd=os.path.dirname(config_file),
    env=env
  )
  p.wait()

  if p.returncode != 0:
    raise Exception(f"cwl-runner finished with exit code {p.returncode}") 
