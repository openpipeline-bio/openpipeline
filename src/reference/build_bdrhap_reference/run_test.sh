#!/bin/bash



## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

echo "> Running $meta_functionality_name, writing to $tmpdir."
$meta_executable \
  --genome_fasta "$meta_resources_dir/reference_gencodev41_chr1/reference.fa.gz" \
  --transcriptome_gtf "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" \
  --output "$tmpdir/myreference.tar.gz"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/myreference.tar.gz" ]] && echo "Output tar file could not be found!" && exit 1

echo "> Test succeeded!"