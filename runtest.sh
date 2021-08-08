mirnabed="test/data/sample.mirna.bed"
utrbed="test/data/sample.utr.bed"
vcf="test/data/sample.vcf.gz"
fasta="test/data/"
output_directory="$(dirname $mirnabed)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run extract-sequences.nf \
	--mirnabed $mirnabed \
  --utrbed $utrbed \
	--vcf $vcf \
	--fasta_dir $fasta \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
