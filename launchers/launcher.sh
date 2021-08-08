#!/usr/bin/env bash
cd ../


input_mirna="real-data/chrom22_mut/CHROM22.mirna.fa"
input_utr="real-data/chrom22_ref/CHROM22.utr.fa"
output_directory="real-data/chrom22_mut/results-2"

nextflow run compare-miRNA-pairs.nf \
	--mirna_fa $input_mirna \
	--utr_fa	$input_utr \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
