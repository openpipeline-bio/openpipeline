#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_sub" ]; then
  extra_params+=( "--sub $par_sub" )
fi
  
if [ ! -z "$par_ems" ]; then
  extra_params+=( "--ems $par_ems" )
fi

if [ ! -z "$par_dbl" ]; then
  extra_params+=( "--dbl $par_dbl" )
fi

if [ ! -z "$par_vcf_known" ] ; then
  extra_params+=( "--vcf $par_vcf_known" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

scSplit_loc = "/usr/local/lib/python3.6/site-packages/scSplit"

python3 ${scSplit_loc}/scSplit count --vcf $par_vcf --bam $par_bam \
        --bar $par_bar --tag $par_tag --ref $par_ref --alt $par_alt \
        --com $par_com --out $par_output
python3 ${scSplit_loc}/scSplit run --ref ${par_output}$ref \
        --alt ${par_output}$alt --out $par_output --num $par_num ${extra_params[@]}

if [ "$par_geno" = true ]; then
  python3 ${scSplit_loc}/scSplit genotype --ref ${par_output}$ref\
        --alt ${par_output}$alt --psc ${par_output}${par_psc} $par_output
fi

