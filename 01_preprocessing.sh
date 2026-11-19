# QC cleaning with Fastp
for f in ./raw_dat/*_R1_001.fastq.gz; do
	base=$(basename "$f" _R1_001.fastq.gz)
	fastp \
		-i ./raw_dat/${base}_R1_001.fastq.gz -I ./raw_dat/${base}_R2_001.fastq.gz \
		-o ./cleaned_fastq/${base}_R1.clean.fastq.gz -O ./cleaned_fastq/${base}_R2.clean.fastq.gz \
		--detect_adapter_for_pe \
		-h ../res/fastp_reports/${base}.html \
		-j ../res/fastp_reports/${base}.json \
		-q 20 -l 25 \
		--trim_poly_g \
		--thread 16
done

# Mapping with BWA-MEM2
index_file="/home/tongta/shared/ref/GRCh38.p14/bwa_index/GRCh38.p14"
clean_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/dat/cleaned_fastq"
output_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/res/mapped_reads"

for f in "$clean_dir"/*_R1.clean.fastq.gz; do
	base=$(basename "$f" _R1.clean.fastq.gz)
	fw_read="$clean_dir/${base}_R1.clean.fastq.gz"
	rv_read="$clean_dir/${base}_R2.clean.fastq.gz"

	bwa-mem2 mem -M -t 16 \
		-R "@RG\tID:$base\tSM:$base\tPL:ILLUMINA" \
		"$index_file" "$fw_read" "$rv_read" \
		| samtools sort -@ 4 -o "$output_dir/${base}.sorted.bam"

done

# Indexing bam files
for f in *.sorted.bam; do samtools index "$f"; done

# Mark duplicate
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/res"
run_picard() {
	sample=$1
	input_file="${prj_dir}/mapped_reads/${sample}.sorted.bam"
	output_file="${prj_dir}/mark_dup/${sample}.markdup.bam"
	matrix_file="${prj_dir}/mark_dup/${sample}.matrics.txt"

	picard MarkDuplicates \
		I="$input_file" \
		O="$output_file" \
		M="$matrix_file" \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		TMP_DIR="${prj_dir}/tmp"
}

export -f run_picard
export prj_dir

parallel -j 4 run_picard {} ::: $(for f in "${prj_dir}"/mapped_reads/*.sorted.bam; do basename "$f" .sorted.bam; done)

# Filter lowQC mapped reads
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/res/"
run_filter() {
	sample=$1
	input_file="${prj_dir}/mark_dup/${sample}.markdup.bam"
	output_file="${prj_dir}/filtered_bam/${sample}.filtered.bam"

	samtools index "$input_file"
	chroms=$(samtools idxstats "$input_file" | cut -f1 | grep -v -e '^chrM$' -e '^\*$')

	samtools view -@ 4 -b -F 1804 -f 2 -q 30 "$input_file" $chroms > "$output_file"
	samtools index "$output_file"
}

export -f run_filter
export prj_dir

parallel -j 4 run_filter {} ::: $(for f in "${prj_dir}"/mark_dup/*.markdup.bam; do basename "$f" .markdup.bam; done)

# filter blacklist
prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/res"
blacklist_file="/home/tongta/shared/ref/hg38-blacklist.v2.bed"
output_dir="${prj_dir}/final_bam"

for f in "${prj_dir}"/filtered_bam/*.filtered.bam; do
	sample=$(basename "$f" .filtered.bam)
	bedtools intersect -v -abam "$f" -b "$blacklist_file" > "${output_dir}/${sample}.final.bam"
	samtools index "${output_dir}/${sample}.final.bam"
done

# Tn5 shift correction
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

run_shift() {
	sample=$1
	input_file="${prj_dir}/res/final_bam/${sample}.final.bam"
	output_file="${prj_dir}/res/shifted_bam/${sample}.shifted.bam"
	export TMPDIR="${prj_dir}/tmp"

	alignmentSieve --ATACshift \
		-b "$input_file" \
		-o "$output_file" \
		--numberOfProcessors 4 

	samtools sort -@ 4 -T "${prj_dir}/tmp/${sample}.sort_tmp" \
		-o "${output_file%.bam}.sorted.bam" "$output_file"
	samtools index "${output_file%.bam}.sorted.bam"
}

export -f run_shift
export prj_dir

parallel -j 4 run_shift {} ::: $(for f in "${prj_dir}"/res/final_bam/*.final.bam; do basename "$f" .final.bam; done)


# Create BigWig files
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

for bam in "$prj_dir"/res/shifted_bam/*.shifted.sorted.bam; do
	sample=$(basename "$bam" .shifted.sorted.bam)
	bamCoverage -b "$bam" -o "$prj_dir/res/bigwig/${sample}.bw" \
		--normalizeUsing CPM \
		--binSize 10 \
		-p max
done

bigwigAverage -b C260409012-bAT.bw
C260409013-bAT.bw
C260409014-bAT.bw
C260409021-bAT.bw
C260409022-bAT.bw
C260409023-bAT.bw -o ./acute_batch2_avg.bw --binSize 10 -p 24

bigwigAverage -b Acute1.bw
Acute2.bw
Acute3.bw
C260409012-bAT.bw
C260409013-bAT.bw
C260409014-bAT.bw
C260409021-bAT.bw
C260409022-bAT.bw
C260409023-bAT.bw -o ./acute_avg.bw --binSize 10 -p 24

bigwigAverage -b C260409015-bAT.bw
C260409016-bAT.bw
C260409017-bAT.bw
C260409024-bAT.bw
C260409025-bAT.bw
C260409026-bAT.bw -o ./chronic_avg.bw --binSize 10 -p 24

bigwigAverage -b C260409009-bAT.bw
C260409010-bAT.bw
C260409011-bAT.bw
C260409018-bAT.bw
C260409019-bAT.bw
C260409020-bAT.bw
Control1.bw
Control2.bw
Control3.bw -o ./control_acute_avg.bw --binSize 10 -p 24

bigwigAverage -b C260409009-bAT.bw
C260409010-bAT.bw
C260409011-bAT.bw
C260409018-bAT.bw
C260409019-bAT.bw
C260409020-bAT.bw -o ./control_chronic_avg.bw --binSize 10 -p 24
