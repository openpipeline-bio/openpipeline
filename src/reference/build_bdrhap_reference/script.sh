#!/bin/bash

set -eo pipefail

## VIASH START
par_genome_fasta="temp/reference.fa.gz"
par_transcriptome_gtf="temp/reference.gtf.gz"
par_output="temp/reference_star.tar.gz"
meta_n_proc=20
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_resources_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

meta_n_proc="${meta_n_proc:-1}"

# process params
extra_params=( )

if [ ! -z "$meta_n_proc" ]; then 
  extra_params+=( "--runThreadN $meta_n_proc" )
fi

echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmpdir/genome.fa"
unpigz -c "$par_transcriptome_gtf" > "$tmpdir/transcriptome.gtf"

echo "> Building star index"
mkdir "$tmpdir/out"
STAR \
  --runMode genomeGenerate \
  --genomeDir "$tmpdir/out" \
  --genomeFastaFiles "$tmpdir/genome.fa" \
  --sjdbGTFfile "$tmpdir/transcriptome.gtf" \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 11 \
  "${extra_params[@]}"

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "$tmpdir/out" .