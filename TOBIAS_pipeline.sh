# ATACorrect
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

run_tobias_atacorrect() {

	bam_file="$1"
	genome_file="/home/tongta/shared/ref/GRCh38.p14/bwa_index/GRCh38.p14.genome.fa"
	peak_file="${prj_dir}/res/macs3_peaks/consensus_peaks.bed"
	output_dir="${prj_dir}/dat/TOBIAS_atacorrected"

	TOBIAS ATACorrect \
		--bam "$bam_file" \
		--genome "$genome_file" \
		--peaks  "$peak_file" \
		--outdir "$output_dir" \
		--cores 4
}
export -f run_tobias_atacorrect
export prj_dir

parallel -j 6 run_tobias_atacorrect {} ::: "${prj_dir}"/dat/final_bam/*.final.bam


# ScoreBigwig
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

run_tobias_scorebigwig() {

	signal_file="$1"
	sample=$(basename "$signal_file" .final_corrected.bw)
	output_file="${prj_dir}/dat/TOBIAS_scorebigwig/${sample}_footprints.bw"
	peak_regions="${prj_dir}/res/macs3_peaks/consensus_peaks.bed"
	
	TOBIAS ScoreBigwig \
		--signal "$signal_file" \
		--regions "$peak_regions" \
		--output "$output_file" \
		--cores 4
}

export -f run_tobias_scorebigwig
export prj_dir

parallel -j 4 run_tobias_scorebigwig {} ::: "${prj_dir}"/dat/TOBIAS_atacorrected/*.final_corrected.bw


# BINDetect
#!/bin/bash

prj_dir="/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn"

run_tobias_bindetect() {

	case_file="$1"
	ctrl_file="$2"
	case_name="$3"
	ctrl_name="$4"
	motif_file="${prj_dir}/dat/jaspar_human_core.meme"
	peak_file="${prj_dir}/res/macs3_peaks/consensus_peaks.bed"
	genome_file="/home/tongta/shared/ref/GRCh38.p14/bwa_index/GRCh38.p14.genome.fa"
	output_dir="${prj_dir}/dat/TOBIAS_bindetect/${case_name}_vs_${ctrl_name}"

	TOBIAS BINDetect \
		--motifs "$motif_file" \
		--signals "$case_file" "$ctrl_file" \
		--genome "$genome_file" \
		--peaks "$peak_file" \
		--cond_names "$case_name" "$ctrl_name" \
		--outdir "$output_dir" \
		--cores 16
}

run_tobias_bindetect "${prj_dir}/dat/TOBIAS_scorebigwig/acute_avg.bw" "${prj_dir}/dat/TOBIAS_scorebigwig/control_acute_avg.bw" "Acute" "Control"
run_tobias_bindetect "${prj_dir}/dat/TOBIAS_scorebigwig/chronic_avg.bw" "${prj_dir}/dat/TOBIAS_scorebigwig/control_chronic_avg.bw" "Chronic" "Control"

