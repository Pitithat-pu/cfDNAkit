cfDNAkit_DIR="/home/puranach/Documents/Rscripts/cfDNAkit/cnv-size-based-zscore"
binsize=50
sample_id="LB-I023_004"
sample_bamfile="/icgc/dkfzlsdf/project/inform/liquid_biopsies/sequencing/exon_sequencing/view-by-pid/${sample_id}/tumor01/paired/merged-alignment/tumor01_${sample_id}_merged.mdup.bam"
sample_bamindex="/icgc/dkfzlsdf/project/inform/liquid_biopsies/sequencing/exon_sequencing/view-by-pid/${sample_id}/tumor01/paired/merged-alignment/tumor01_${sample_id}_merged.mdup.bam"
output_folder="/abi/data2/puranach/ct_DNA/inform_liquidBiopsy/test_with_new_deftaf"
delta_f_control_file="/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/deltaf_Agilent6v2_withoutUTR/BH01_rmdup_paired_mapped_deltaf_table_n10_50_mincov5000.csv"

capture_targets_file="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/Agilent6v2withUTRs_plain.bed.gz"
#capture_targets_file="/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/Nimblegen_SeqCap_EZ_Exome_v3_plain.bed.gz"

BAF_file="/abi/data2/puranach/ct_DNA/inform_liquidBiopsy/BAFFiles/tumor01_${sample_id}_merged.mdup_BAF_FILE.tsv.gz"
log_dir="/abi/data2/puranach/ct_DNA/inform_liquidBiopsy/qsublog/"

filename=$(basename ${sample_bamfile%.*})



job_id=$(echo "module load R/3.4.0; Rscript ${cfDNAkit_DIR}/initial_zscore_segmentation.R ${sample_bamfile} ${sample_bamindex} ${output_folder} ${binsize} ${delta_f_control_file} ${capture_targets_file}" | \
bsub -J initial_zscore_segmentation_${filename} -M 20GB -W 120:00 -oo ${log_dir} -eo ${log_dir}  | cut -d '.' -f 1)


job_id=$(echo "module load R/3.4.0; Rscript ${cfDNAkit_DIR}/BAF_segmentation.R ${BAF_file}  ${output_folder}/${filename}_initzscore_${binsize}_segmentation.csv  ${output_folder}" | \
bsub -w \"done(${job_id})\" -J BAF_segmentation_${filename} -M 20GB -W 120:00 -oo ${log_dir} -eo ${log_dir} | cut -d '.' -f 1)


job_id=$(echo "module load R/3.4.0; Rscript ${cfDNAkit_DIR}/segmentbyBAF.R ${output_folder}/${filename}_initzscore_${binsize}_reads_extracted_info.csv ${delta_f_control_file} ${output_folder}/${filename}_initzscore_${binsize}_segmentation_baf_segmentation.csv ${output_folder}" | \
qsub -w \"done(${job_id})\" -J segmentbyBAF_${filename} -M 20GB -W 120:00 -oo ${log_dir} -eo ${log_dir} | cut -d '.' -f 1)


#echo "/home/puranach/Documents/Rscripts/cfDNAkit/BAFzscore_workflow.sh" | qsub -l walltime=120:00:00,mem=20GB -N BAFzscore_workflow -e /home/puranach/Documents/Rscripts/cfDNAkit/qsublog/ -o /home/puranach/Documents/Rscripts/cfDNAkit/qsublog/
