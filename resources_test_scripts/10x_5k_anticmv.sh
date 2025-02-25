#!/bin/bash

set -eo pipefail


# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_anticmv
OUT=resources_test/$ID

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# dataset page:
# https://www.10xgenomics.com/resources/datasets/integrated-gex-totalseqc-and-tcr-analysis-of-connect-generated-library-from-5k-cmv-t-cells-2-standard

# check whether reference is available
reference_dir="resources_test/reference_gencodev41_chr1/"
genome_tar="$reference_dir/reference_cellranger.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
    exit 1
fi

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/5k_human_antiCMV_T_TBNK_connect_Multiplex"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_fastqs.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit head -n 200000 "$input" | gzip > "$output"
  fi
}

orig_sample_id="5k_human_antiCMV_T_TBNK_connect"

seqkit_head "$tar_dir/gex_1/${orig_sample_id}_GEX_1_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_GEX_1_subset_S1_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/gex_1/${orig_sample_id}_GEX_1_S1_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_GEX_1_subset_S1_L001_R2_001.fastq.gz"

seqkit_head "$tar_dir/ab/${orig_sample_id}_AB_S2_L004_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_AB_subset_S2_L004_R1_001.fastq.gz"
seqkit_head "$tar_dir/ab/${orig_sample_id}_AB_S2_L004_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_AB_subset_S2_L004_R2_001.fastq.gz"

seqkit_head "$tar_dir/vdj/${orig_sample_id}_VDJ_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_VDJ_subset_S1_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/vdj/${orig_sample_id}_VDJ_S1_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_VDJ_subset_S1_L001_R2_001.fastq.gz"

# download immune panel fasta if needed
feature_reference="$raw_dir/feature_reference.csv"
if [[ ! -f "$feature_reference" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv" -O "$feature_reference"
fi

# download vdj reference if needed
vdj_ref="$raw_dir/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"
if [[ ! -f "$vdj_ref" ]]; then
  wget "https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz" -O "$vdj_ref"
fi


# Run mapping pipeline
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - "${orig_sample_id}_GEX_1_subset"
    - "${orig_sample_id}_AB_subset"
    - "${orig_sample_id}_VDJ_subset"
  library_type:
    - "Gene Expression"
    - "Antibody Capture"
    - "VDJ"

gex_reference: "$genome_tar"
vdj_reference: "$vdj_ref"
feature_reference: "$feature_reference"
publish_dir: "$OUT/processed"
HERE

nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config \
  -c src/workflows/utils/errorstrat_ignore.config

# Convert to h5mu
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/processed/10x_5k_anticmv.cellranger_multi.output.output"
publish_dir: "$OUT/"
output: "$orig_sample_id.h5mu"
HERE

nextflow \
  run . \
  -main-script target/nextflow/convert/from_cellranger_multi_to_h5mu/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config

# run qc workflow
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/$orig_sample_id.h5mu"
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
publish_dir: "$OUT/"
output: "${orig_sample_id}_qc.h5mu"
HERE

nextflow \
  run . \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config


# Run full pipeline
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/${orig_sample_id}_qc.h5mu"
publish_dir: "$OUT/"
output: "${orig_sample_id}_mms.h5mu"
HERE

nextflow \
  run . \
  -main-script src/workflows/multiomics/full_pipeline/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config

# create fastqc directory
fastqc_dir="$OUT/fastqc"
mkdir -p "$fastqc_dir"

./target/docker/qc/fastqc/fastqc \
  --input "$raw_dir" \
  --mode "dir" \
  --output "$fastqc_dir"