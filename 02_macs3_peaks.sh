# Peak calling
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

run_macs3() {
	sample=$1
	input_file="${prj_dir}/res/shifted_bam/${sample}.shifted.sorted.bam"
	output_dir="${prj_dir}/res/macs3_peaks"

	macs3 callpeak \
		-t $input_file \
		-f BAMPE \
		-g hs \
		-n "$sample" \
		--outdir "$output_dir" \
		--keep-dup all \
		-p 0.01 \
		--tempdir "${prj_dir}/tmp"
}

export -f run_macs3
export prj_dir

parallel -j 8 run_macs3 {} ::: $(for f in "${prj_dir}"/res/shifted_bam/*.shifted.sorted.bam; do basename "$f" .shifted.sorted.bam; done)

# Peak QC ###########

# Calculate FRiP score using
# Filtering out the blacklist peaks
for f in ./macs3_peaks/*_peaks.narrowPeak; do
	bedtools intersect -v -a "$f" -b ~/shared/ref/hg38-blacklist.v2.bed > "${f$.narrowPeak}.filtered.narrowPeak"
done

# Consensus peaks
cat ./macs3_peaks/*_peaks.filtered.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n \
	| bedtools merge -i - > ./macs3_peaks/consensus_peaks.bed

awk 'BEGIN{OFS=t; print GeneID, Chr, Start, End, Strand}
	{print peakNR, , +1, , .}' ./macs3_peaks/consensus_peaks.bed > ./macs3_peaks/consensus_peaks.saf

# FeatureCounts
featureCounts -a consensus_peaks.saf -F SAF -o counts.txt -p --countReadPairs -T 8 \
	./shifted_bam/*.shifted.bam

# Extract FRiP scores from FeatureCounts
awk -F'\t' '
NR==1 {for(i=2;i<=NF;i++) name[i]=}
==Assigned              {for(i=2;i<=NF;i++) assigned[i]=}
==Unassigned_NoFeatures {for(i=2;i<=NF;i++) nofeat[i]=}
==Unassigned_Ambiguity  {for(i=2;i<=NF;i++) ambig[i]=}
END {
  for(i=2;i<=length(name)+1;i++) {
    total = assigned[i] + nofeat[i] + ambig[i]
    printf %stFRiP=%.3ftAssigned=%dtTotal=%dn, name[i], assigned[i]/total, assigned[i], total
  }
}' counts.txt.summary > FRiP_scores.txt

# ataqv
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

conda activate ataqv_env
run_ataqv() {
	sample=$1
	bamfile="${prj_dir}/res/shifted_bam/${sample}.shifted.sorted.bam"
	input_peak="${prj_dir}/res/macs3_peaks/${sample}_peaks.filtered.narrowPeak"
	tss_file="/home/tongta/shared/ref/gencode/gencode_tss.bed"
	metric_file="${prj_dir}/res/qc/peak_ataqv/${sample}.ataqv.json"
	output_file="${prj_dir}/res/qc/peak_ataqv/${sample}.ataqv.out"

	ataqv --peak-file "$input_peak" \
		--tss-file "$tss_file" \
		--name "$sample" \
		--metrics-file "$metric_file" \
		human "$bamfile" > "$output_file"
}

export -f run_ataqv
export prj_dir

samples=$(for f in "${prj_dir}"/res/shifted_bam/*.shifted.sorted.bam; do
	basename "$f" .shifted.sorted.bam
done)

parallel -j 8 run_ataqv {} ::: $samples

# Summarise QC
mkarv qc_report ./qc/peak_ataqv/*.json

conda deactivate
