options(digits=10)
binsize = 100 #kb
package_dir = "/abi/data2/puranach/cfDNAkit/R/"
resource_dir=paste0(package_dir,"/resources/")
main_function_file = paste0(package_dir,"/kit_functions/main_functions.R")
runtime_variables_filename = "pon_runtime_variables.R"
### Resource
qdnaseq_sliding_windows_RDS = paste0(resource_dir,"/AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds")
duke_blacklist_region = paste0(resource_dir,"/wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz")
dac_blacklist_region = paste0(resource_dir,"/wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz")
create_pon_outdir_prefix = "create_pon"
rerun_readbam = FALSE

what <- c("qname","rname","pos", "isize")
maximum_length = 400
minimum_length = 1

# cutoff_zscore = c(3,-3)
short_length_range = c(100,150)
long_length_range = c(151,250)

fragment_perwindow_lower_bound_percent = 0
fragment_perwindow_upper_bound_percent  = 100

minimum_coverage_per_bin = 100
sampling_size = 1 ## used when running pooled sample file.
