options(digits=10)
package_dir = "/abi/data2/puranach/cfDNAkit/R/"
binsize = 100 #kb
minimum_coverage_per_bin = 30
rerun_readbam = TRUE  ## read bam file again
resource_dir=paste0(package_dir,"/resources/")
main_function_file = paste0(package_dir,"/kit_functions/main_functions.R")
runtime_variables_filename = "runtime_variables_chrarm.R"
output_regions_bed="/abi/data2/puranach/cfDNAkit/R/resources/chromosome_arms.bed"
### Resource
chrLength_file = paste0(resource_dir,"hg19_chrTotalLength.tsv")
control_data_dir = paste0("/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/fragment_length_analysis/K34R-WYE4XG/createPON_",binsize,"k/merged_PON/")

# delta_f_control_file= paste0(control_data_dir,"BH01_rmdup_paired_mapped_deltaf_table_n30_1000.csv")
control_shortread_file = paste0(control_data_dir,"PON_shortread.tsv")
control_longread_file = paste0(control_data_dir,"PON_longread.tsv")
control_zscore_file = paste0(control_data_dir,"PON_zscore.tsv")
control_total_fragment_file = paste0(control_data_dir,"PON_total_fragment.tsv")
control_density_file = paste0(resource_dir,"BH01_rmdup_paired_mapped_Fragment-length_report_50_insert_size_density.csv")
qdnaseq_sliding_windows_RDS = paste0(resource_dir,"/AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds")
duke_blacklist_region = paste0(resource_dir,"/wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz")
dac_blacklist_region = paste0(resource_dir,"/wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz")
centromere_region = paste0(resource_dir,"/hg19_centromere.tsv.gz")


plot_insertsize_density_with_PON = TRUE
calculate_mad = TRUE
calculate_gwzscore = TRUE

what <- c("qname","rname","pos", "isize")
extract_density_outdir_prefix="density_info"
init_zscore_outdir_prefix="initzscore_PON"
coverage_zscore_outdir_prefix="coverage_zscore_PON"
segmentation_outdir_prefix="zscore-segmentation_PON"
maximum_length = 400
minimum_length = 1
# cutoff_zscore = c(1.96,-1.96)
cutoff_zscore = c(1.645,-1.645)
# cutoff_zscore = c(3,-3)
short_length_range = c(100,150)
long_length_range = c(151,250)
plot_upperbound = 10
plot_lowerbound = -10
min_heterozygous_SNPs = 5
smooth_moving_window = 10
