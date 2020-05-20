options(digits=10)
### run specific configuration, modify as needed
binsize = 1000 #kb
package_dir = "/abi/data2/puranach/cfDNAkit/R/"
resource_dir=paste0(package_dir,"/resources/")
main_function_file = paste0(package_dir,"/kit_functions/main_functions.R")
runtime_variables_filename = paste0("pon_runtime_variables_",binsize,"k.R")
rerun_readbam = FALSE
create_pon_outdir_prefix = "create_pon"
###


### software behaviour config
maximum_length = 1000
minimum_length = 1
short_length_range = c(100,150)
long_length_range = c(151,250)
minimum_coverage_per_bin = 500
sampling_size = 1 ## changed when running a pooled sample bamfile. this tell number of sample to split into
###

### Resource file, no need to change
qdnaseq_sliding_windows_RDS = paste0(resource_dir,"/AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds")
duke_blacklist_region = paste0(resource_dir,"/wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz")
dac_blacklist_region = paste0(resource_dir,"/wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz")
centromere_region = paste0(resource_dir,"/hg19_centromere.tsv.gz")
###

### Rsamtools configs
what <- c("qname","rname","pos", "isize")
###
