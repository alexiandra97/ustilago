#!/bin/bash

sample=$1

# Define variables
working_dir=$PWD
input_data=${working_dir}/fastq
ref_genome="${working_dir}/ref/UCM520_ref.fa"
threads_trimming=4
threads_mapping=16
base_quality=20

[ -d clean ] || mkdir clean
[ -d reports ] || mkdir reports
[ -d bams ] || mkdir bams

# read group
RGID="$sample"
RGLB='wgs-lib'
RGPL='ILLUMINA'
RGPU="$sample"

fastp -w $threads_trimming \
	-i $input_data/${sample}_1.fq.gz \
	-I $input_data/${sample}_2.fq.gz \
	-o $working_dir/clean/${sample}_1.fq.gz \
	-O $working_dir/clean/${sample}_2.fq.gz \
	-j $working_dir/reports/${sample}.json \
	-h $working_dir/reports/${sample}.html

bwa mem -t $threads_trimming \
	-R "@RG\tID:${RGID}\tSM:${sample}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
	$ref_genome \
	$working_dir/clean/${sample}_1.fq.gz \
	$working_dir/clean/${sample}_2.fq.gz \
| samtools sort -m 8G \
	-o $working_dir/bams/${sample}_raw.bam \
			-
samtools view -b -F 2048 -f 2 -q ${base_quality} \
	-o $working_dir/bams/${sample}_paired.bam \
	$working_dir/bams/${sample}_raw.bam

gatk MarkDuplicates \
	--MAX_RECORDS_IN_RAM 10000000 \
	--CREATE_INDEX true \
	-I $working_dir/bams/${sample}_paired.bam \
	-O $working_dir/bams/${sample}.bam \
	-M $working_dir/bams/${sample}_metrics.txt
