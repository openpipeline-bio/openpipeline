#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT=resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# download metrics summary
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv \
  --output "${OUT}_metrics_summary.csv"

# download counts h5 file
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output "${OUT}_filtered_feature_bc_matrix.h5"

# download counts matrix tar gz file
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz \
  --output "${OUT}_filtered_feature_bc_matrix.tar.gz"

# extract matrix tar gz
mkdir -p "${OUT}_filtered_feature_bc_matrix"
tar -xvf "${OUT}_filtered_feature_bc_matrix.tar.gz" \
  -C "${OUT}_filtered_feature_bc_matrix" \
  --strip-components 1
rm "${OUT}_filtered_feature_bc_matrix.tar.gz"

# convert 10x h5 to h5mu
target/docker/convert/from_10xh5_to_h5mu/from_10xh5_to_h5mu \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --output "${OUT}_filtered_feature_bc_matrix.h5mu"

# run single sample
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/process_rna/singlesample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_uss \
  --input "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --output "`basename $OUT`_uss.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run multisample
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/process_rna/multisample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_ums \
  --input "${OUT}_uss.h5mu" \
  --output "`basename $OUT`_ums.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run integration
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/integration/multimodal_integration/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_mms \
  --input "${OUT}_ums.h5mu" \
  --output "`basename $OUT`_mms.h5mu" \
  --publishDir `dirname $OUT` \
  -resume