cfDNAkit_dir="/abi/data2/puranach/cfDNAkit/R/"
config_file=${cfDNAkit_dir}/config_files/generate_control/generate_PoN_chrarm.R
runtime_variable_filename="pon_chrarm_runtime_variables.R"
pid_dir="/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/whole_genome_sequencing/results_per_pid/"
analysis_dir="/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/fragment_length_analysis/"
cfDNAkit_outdir_prefix="createPON_chrarm"
sample_ids=" K34R-WYE4XG "
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
		CREATE_PON_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/PON_generation/create_size-based-pon.R ${runtime_variables_file} ${config_file}"


		EXTRACT_BAM_INFO_job_id=$(echo "${EXTRACT_BAM_INFO_CMD}" | bsub -J ${sample_id}_extract_bam_info -W 10:00 -M 30GB | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
		echo "Submit EXTRACT_BAM_INFO process : job id ${EXTRACT_BAM_INFO_job_id}"
		sleep 1s
		CREATE_PON_job_id=$(echo "${CREATE_PON_CMD}" | bsub -J ${sample_id}_CREATE_PON -W 3:00 -M 20GB -w "done(${EXTRACT_BAM_INFO_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
echo "Submit CREATE_PON process : job id ${CREATE_PON_job_id}"
		echo "Done: all jobs submitted"
done; done

