#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --cellData "$meta_resources_dir/vireo_test_data/cells.cellSNP.vcf.gz" \
    --nDonor 4 \
    --output vireo_result/

[[ ! -f vireo_result/GT_donors.vireo.vcf.gz ]] && echo "Output donor genotype file could not be found!" && exit 1
[[ ! -f vireo_result/assignment.tsv ]] && echo "Output donor assignment tsv could not be found!" && exit 1
echo ">>> Test finished successfully"
