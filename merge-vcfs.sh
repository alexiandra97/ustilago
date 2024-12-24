#!/bin/bash

samples="sample1 sample2 sample3"

# Define Input/Output variables
working_dir=${PWD}
input_dir="${working_dir}/vcfs"
count_dir="${working_dir}/bams/readcounts"
output_dir="${working_dir}/vcfs"
referentni="${working_dir}/ref/UCM520_ref.fa"
gff_file="${working_dir}/ref/genomic.gff"
vcf_files=""

[ -d ${count_dir} ] || mkdir ${count_dir}

for sample in $samples
do
	vcf_files+=" $input_dir/${sample}.norm.vcf.gz"
done

# Merge vcf files
bcftools merge \
	$vcf_files \
	-o ${output_dir}/all_samples_merged.vcf

vcf_count_in=${output_dir}/all_samples_merged_in.vcf
vcf_count_out=${output_dir}/all_samples_merged_out.vcf

cp ${output_dir}/all_samples_merged.vcf $vcf_count_in

for sample in $samples
do
	vcf-readcount-annotator -s ${sample} \
		-o $vcf_count_out \
		$vcf_count_in \
		${count_dir}/${sample}.counts DNA
	mv $vcf_count_out $vcf_count_in
done

bcftools +setGT $vcf_count_in -o $vcf_count_out -- -t q -n X -i 'GT="." && FORMAT/DP>15'

# Export genes from GFF to "BED format"
awk -F"\t" '
	BEGIN {
		OFS=IFS
	}
	$3=="gene" {
		split($9,id,";");
		gene=id[1];
		gsub(/ID=gene-/,"",gene);
		print $1,$4,$5,gene,$6,$7}' \
	$gff_file > ${output_dir}/ustilago-genes.tsv

# Annotation vcf files
bcftools annotate \
	-a ${output_dir}/ustilago-genes.tsv \
	-c CHROM,FROM,TO,GENE \
	-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
	$vcf_count_out \
	> ${output_dir}/${samples_name}_merged_genes.vcf

bgzip -f ${output_dir}/${samples_name}_merged_genes.vcf
bcftools index ${output_dir}/${samples_name}_merged_genes.vcf.gz
