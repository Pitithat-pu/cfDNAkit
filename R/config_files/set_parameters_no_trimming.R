options(digits=10)
binsize = 100 #kb
package_dir = "/abi/data2/puranach/cfDNAkit/R/"
resource_dir=paste0(package_dir,"/resources/")
main_function_file = paste0(package_dir,"/kit_functions/main_functions.R")
runtime_variables_filename = "runtime_variables.R"
### Resource
chrLength_file = "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_chrTotalLength.tsv"
control_data_dir = "/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/deltaf_whole_genome/size-based_zscore_2019Sep19_0954/"
delta_f_control_file= paste0(control_data_dir,"BH01_rmdup_paired_mapped_deltaf_table_n40_100.csv")
control_shortread_file = paste0(control_data_dir,"BH01_rmdup_paired_mapped_deltaf_table_n40_100_shortread.corrected.csv")
control_longread_file = paste0(control_data_dir,"BH01_rmdup_paired_mapped_deltaf_table_n40_100_longread.corrected.csv")
control_zscore_file = paste0(control_data_dir,"BH01_rmdup_paired_mapped_deltaf_table_n40_100_zscore.csv")
control_density_file = "/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/BH01_rmdup_paired_mapped_Fragment-length_report_50_insert_size_density.csv"
qdnaseq_sliding_windows_RDS = paste0(resource_dir,"/AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds")
duke_blacklist_region = paste0(resource_dir,"/wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz")
dac_blacklist_region = paste0(resource_dir,"/wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz")
plot_insertsize_density_with_PON = TRUE

what <- c("qname","rname","pos", "isize")
extract_density_outdir_prefix="density_info"
init_zscore_outdir_prefix="initzscore"
segmentation_outdir_prefix="zscore-segmentation"
maximum_length = 400
minimum_length = 1
minimum_coverage_per_bin = 30
# cutoff_zscore = c(1.96,-1.96)
cutoff_zscore = c(1.645,-1.645)
# cutoff_zscore = c(3,-3)
short_length_range = c(100,150)
long_length_range = c(151,250)
plot_upperbound = 10
plot_lowerbound = -10
min_heterozygous_SNPs = 5

fragment_perwindow_lower_bound_percent = 0
fragment_perwindow_upper_bound_percent  = 100
smooth_moving_window = 10
rerun_readbam = TRUE
