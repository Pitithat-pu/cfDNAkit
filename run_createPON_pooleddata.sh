cfDNAkit_dir="/abi/data2/puranach/cfDNAkit/R/"
config_file=${cfDNAkit_dir}/config_files/generate_control/generate_PoN_pooleddata_100k.R
runtime_variable_filename="pon_runtime_variables.R"
pid_dir="/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/"
cfDNAkit_outdir_prefix="createPON_100k"
load_R_cmd="module load R/3.5.1"
bamfile=${pid_dir}/BH01_rmdup_paired_mapped.bam

filename=$(basename ${bamfile} .bam);
sample_output_dir=${pid_dir}/${cfDNAkit_outdir_prefix}/${filename};
runtime_variables_file="${sample_output_dir}/${runtime_variable_filename}"
echo "Create jobs to submit for sample file ${filename}.bam"
EXTRACT_BAM_INFO_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/extract_bam_alignment_info.R ${bamfile} ${config_file} ${sample_output_dir}"
CREATE_PON_CMD="${load_R_cmd};Rscript ${cfDNAkit_dir}/PON_generation/create_size-based-pon_from_pooleddata.R ${runtime_variables_file} ${config_file}"




EXTRACT_BAM_INFO_job_id=$(echo "${EXTRACT_BAM_INFO_CMD}" | bsub -J ${bamfile}_extract_bam_info -W 10:00 -M 30GB | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
echo "Submit EXTRACT_BAM_INFO process : job id ${EXTRACT_BAM_INFO_job_id}"
sleep 1s
CREATE_PON_job_id=$(echo "${CREATE_PON_CMD}" | bsub -J ${bamfile}_CREATE_PON -W 3:00 -M 20GB -w "done(${EXTRACT_BAM_INFO_job_id})" | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
echo "Submit CREATE_PON process : job id ${CREATE_PON_job_id}"
echo "Done: all jobs submitted"


