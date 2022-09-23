#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# The output folder
OUT="resources_test/concat_test_data/"

# create it if it doesn't exist already
[ -d "$OUT" ] || mkdir -p "$OUT"

echo "> Downloading files"
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-arc/1.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5 \
  --output "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5"

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-arc/1.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.h5 \
  --output "${OUT}/human_brain_3k_filtered_feature_bc_matrix.h5"

echo "> Converting to h5mu"
bin/viash run src/convert/from_10xh5_to_h5mu/config.vsh.yaml -- \
  --input "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5" \
  --output "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5mu"

bin/viash run src/convert/from_10xh5_to_h5mu/config.vsh.yaml -- \
  --input "$OUT/human_brain_3k_filtered_feature_bc_matrix.h5" \
  --output "$OUT/human_brain_3k_filtered_feature_bc_matrix.h5mu"

echo "> Subsetting datasets"
bin/viash run src/filter/subset_h5mu/config.vsh.yaml -p docker -- \
  --input "$OUT/human_brain_3k_filtered_feature_bc_matrix.h5mu" \
  --output "$OUT/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu" \
  --number_of_observations 2000

bin/viash run src/filter/subset_h5mu/config.vsh.yaml -p docker -- \
  --input "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5mu" \
  --output "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu" \
  --number_of_observations 2000

echo "Removing temp files"
rm "$OUT/human_brain_3k_filtered_feature_bc_matrix.h5" \
  "$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5" \
  "${OUT}/human_brain_3k_filtered_feature_bc_matrix.h5mu" \
  "${OUT}/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5mu"

echo "> Running concat component"
bin/viash run src/integrate/concat/config.vsh.yaml -- \
  --input "$OUT/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu,$OUT/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu" \
  --sample_names "mouse,human" \
  --output "$OUT/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"
