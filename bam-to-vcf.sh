#!/bin/bash

sample=$1

# Define variables
working_dir=$PWD
input_dir="${working_dir}/bams"
count_dir="${working_dir}/bams/readcounts"
output_dir="${working_dir}/vcfs"
ref_genome="${working_dir}/ref/UCM520_ref.fa"
base_quality=20

[ -d ${output_dir} ] || mkdir ${output_dir}
[ -d ${count_dir} ] || mkdir ${count_dir}

freebayes -p 1 \
	-q $base_quality -m 60 --min-coverage 30 \
	-f  $ref_genome \
	${input_dir}/${sample}.bam > ${output_dir}/${sample}.raw.vcf

bam-readcount -w1 -q 60 -b $base_quality --insertion-centric \
	-f $ref_genome ${input_dir}/${sample}.bam \
	> ${count_dir}/${sample}.counts

bcftools filter \
	-e 'FMT/DP < 15' \
	${output_dir}/${sample}.raw.vcf \
| bcftools norm -f $ref_genome \
	--atomize --atom-overlaps . \
	-Oz -o ${output_dir}/${sample}.norm.vcf.gz

bcftools index ${output_dir}/${sample}.norm.vcf.gz
