#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=pbmc_1k_protein_v3
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# dataset page:
# https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-gene-expression-and-cell-surface-protein-3-standard-3-0-0

# download metrics summary
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv \
  -O "${OUT}_metrics_summary.csv"

# download counts h5 file
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  -O "${OUT}_filtered_feature_bc_matrix.h5"

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  -O "${OUT}_raw_feature_bc_matrix.h5"

# download counts matrix tar gz file
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz \
  -O "${OUT}_filtered_feature_bc_matrix.tar.gz"

# extract matrix tar gz
mkdir -p "${OUT}_filtered_feature_bc_matrix"
tar -xvf "${OUT}_filtered_feature_bc_matrix.tar.gz" \
  -C "${OUT}_filtered_feature_bc_matrix" \
  --strip-components 1
rm "${OUT}_filtered_feature_bc_matrix.tar.gz"

# convert 10x h5 to h5mu
target/docker/convert/from_10xh5_to_h5mu/from_10xh5_to_h5mu \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --input_metrics_summary "${OUT}_metrics_summary.csv" \
  --output "${OUT}_filtered_feature_bc_matrix.h5mu"

# run single sample
NXF_VER=21.10.6 nextflow \
  run . \
  -main-script src/workflows/multiomics/rna_singlesample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_uss \
  --input "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --output "`basename $OUT`_uss.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run multisample
NXF_VER=21.10.6 nextflow \
  run . \
  -main-script src/workflows/multiomics/rna_multisample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_ums \
  --input "${OUT}_uss.h5mu" \
  --output "`basename $OUT`_ums.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run integration
NXF_VER=21.10.6 nextflow \
  run . \
  -main-script src/workflows/multiomics/integration/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_mms \
  --input "${OUT}_ums.h5mu" \
  --output "`basename $OUT`_mms.h5mu" \
  --publishDir `dirname $OUT` \
  --obs_covariates sample_id \
  -resume

python <<HEREDOC
import mudata as mu
mudata = mu.read_h5mu("${DIR}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu")
mudata.mod["rna"].write_h5ad("${DIR}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5ad")
HEREDOC