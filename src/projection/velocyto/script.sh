#!/bin/bash

## VIASH START
par_input="resources_test/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"
par_transcriptome="resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf"
par_output="./foo/output.loom"
par_barcode=""
## VIASH END

extra_params=( )

if [ ! -z "$par_barcode" ]; then 
  extra_params+=( "--bcfile=$par_barcode" )
fi

if [ "$par_without_umi" == "true" ]; then
  extra_params+=( "--without-umi" )
fi

output_dir=`dirname "$par_output"`
sample_id=`basename "$par_output" .loom`

echo "$par_input"
echo "$par_output"
echo "$output_dir"
velocyto run \
  "$par_input" \
  "$par_transcriptome" \
  "${extra_params[@]}" \
  --outputfolder "$output_dir" \
  --sampleid "$sample_id" \
  --samtools-threads 4 \
  --samtools-memory 3500
