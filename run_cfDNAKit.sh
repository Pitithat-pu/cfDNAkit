cfDNAkit_dir="/abi/data2/puranach/cfDNAkit/R/"
config_file=${cfDNAkit_dir}/config_files/set_parameters_no_trimming_1MB_PON.R
runtime_variable_filename="runtime_variables.R"
pid_dir="/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/whole_genome_sequencing/results_per_pid/"
analysis_dir="/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/fragment_length_analysis/"
cfDNAkit_outdir_prefix="cfDNAkit"
##The analysis output will be ${analysis_dir}/${sample_id}/${cfDNAkit_outdir_prefix}
# sample_ids=" K34R-UBTXZ9 K34R-WF1SJH K34R-UF2SML K34R-T4SXYQ K34R-S9VXXJ K34R-XDDFED K34R-3TQGFY K34R-C6EYM2 K34R-UG4ZRB K34R-RU832L K34R-LQT3MY K34R-CE5SFA K34R-23MWWD K34R-SFGNT9 K34R-PKN8UU K34R-AABD6B K34R-4KGB9M K34R-EGSBF1 K34R-1PGADQ K34R-7CMR9V K34R-BZMK6K K34R-E75U6S K34R-X3NB4Z K34R-GKD6BS K34R-2VL6V1 K34R-58ULPU K34R-DBD7LR K34R-1DQJQJ K34R-KSGT6G K34R-BFHGCJ K34R-1UWART K34R-8EG32X K34R-WLFZQ5 K34R-WS4PXJ K34R-YUE6ZP K34R-R1U98V K34R-6RLFZV K34R-8AM37B K34R-F7FXM9 K34R-ZPK6WB K34R-8GE8KS K34R-V36XX9 K34R-CNSP79 K34R-314YE3 K34R-S77QY2 K34R-3DLJD9 K34R-Z3WTAN K34R-NTPZA3 K34R-UBTXZ9 K34R-QBKEXL K34R-6GMRNB K34R-NBREXJ K34R-764LG4 K34R-QF63QA K34R-S6EHTR K34R-SHQW63 K34R-GPPYLL K34R-HE38LD K34R-Z6PFY9 K34R-4B8D7J K34R-TMQ2BU K34R-HBAKYV K34R-FJ414V K34R-F84QAF K34R-XVWSLV K34R-KC6K15 K34R-PCC1RG K34R-HYL1T3 K34R-9T1CDU K34R-6WNQBB K34R-6FCKR2 K34R-NJXV7A K34R-35FKAL K34R-DHKDDL K34R-QJRC8U K34R-4UWF2Y " ## exclude K34R-4UWF2Y
sample_ids="K34R-H7ARKZ"
# sample_ids="K34R-4UWF2Y"
# sample_ids="K34R-WYE4XG"
load_R_cmd="module load R/3.5.1"

for sample_id in ${sample_ids};do
	echo ${sample_id};
	alignment_dir="${pid_dir}/${sample_id}/alignment/"
	for bamfile in ${alignment_dir}/*_${sample_id}_merged.mdup.bam; do
		filename=$(basename ${bamfile} .bam);
		sample_output_dir=${analysis_dir}/${sample_id}/${cfDNAkit_outdir_prefix}/${filename};
		runtime_variables_file="${sample_output_dir}/${runtime_variable_filename}"
		echo "Create jobs to submit for sample file ${filename}.bam"
		EXTRACT_BAM_INFO_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/extract_bam_alignment_info.R ${bamfile} ${config_file} ${sample_output_dir}"
		EXTRACT_DENSITY_INFO_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/fragment-length-distribution/extract_density_info.R ${runtime_variables_file} ${config_file}"
		PLOT_DENSITY_INFO_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/fragment-length-distribution/plots/plot_density_info.R ${runtime_variables_file} ${config_file}"
		INITIAL_ZSCORE_SEGMENTATION_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/cnv-size-based-zscore/initial_zscore_segmentation.R ${runtime_variables_file} ${config_file}"
		PSCBS_SEGMENTATION_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/cnv-size-based-zscore/PSCBS_segmentation.R ${runtime_variables_file} ${config_file}"




		EXTRACT_BAM_INFO_job_id=$(echo "${EXTRACT_BAM_INFO_CMD}" | bsub -J ${sample_id}_extract_bam_info -W 10:00 -M 30GB | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
		echo "Submit EXTRACT_BAM_INFO process : job id ${EXTRACT_BAM_INFO_job_id}"
		sleep 1s
		EXTRACT_DENSITY_INFO_job_id=$(echo "${EXTRACT_DENSITY_INFO_CMD}" | bsub -J ${sample_id}_EXTRACT_DENSITY_INFO -W 3:00 -M 10GB -w "done(${EXTRACT_BAM_INFO_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
echo "Submit EXTRACT_DENSITY_INFO process : job id ${EXTRACT_DENSITY_INFO_job_id}"
		PLOT_DENSITY_INFO_job_id=$(echo "${PLOT_DENSITY_INFO_CMD}" | bsub -J ${sample_id}_PLOT_DENSITY_INFO -W 3:00 -M 5GB -w "done(${EXTRACT_DENSITY_INFO_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
		echo "Submit PLOT_DENSITY_INFO process : job id ${PLOT_DENSITY_INFO_job_id}"
		# sleep 1s
		# INITIAL_ZSCORE_SEGMENTATION_job_id=$(echo "${INITIAL_ZSCORE_SEGMENTATION_CMD}" | bsub -J ${sample_id}_INITIAL_ZSCORE_SEGMENTATION -W 5:00 -M 5GB -w "done(${EXTRACT_BAM_INFO_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
		# echo "Submit INITIAL_ZSCORE_SEGMENTATION  process : job id ${INITIAL_ZSCORE_SEGMENTATION_job_id}"
		# sleep 1s
		# PSCBS_SEGMENTATION_job_id=$(echo "${PSCBS_SEGMENTATION_CMD}" | bsub -J ${sample_id}_PSCBS_SEGMENTATION -W 5:00 -M 5GB -w "done(${INITIAL_ZSCORE_SEGMENTATION_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
		# echo "Submit PSCBS_SEGMENTATION  process : job id ${PSCBS_SEGMENTATION_job_id}"
		echo "Done: all jobs submitted"
done; done


