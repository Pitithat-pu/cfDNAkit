.file_exists <- function(files_vector){
  all_exist = TRUE
  for (file in files_vector) {
    if(!file.exists(file)){
      paste(file,"do not existed",sep=" ")
      all_exist = FALSE
    }else{
      paste(file,"existed",sep=" ")
    }
  }
  return(all_exist)
}

args <- commandArgs(T)
if (length(args) < 3){
  stop("Usage: extract_bam_alignment_info.R bamfile* config_file* output_folder* capture_targets_file \n
       (optional capture_targets_file in (/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/))")
}
sample_bamfile <- args[1]
sample_bamindex = paste(sample_bamfile,".bai",sep="")
config_file = args[2]
output_folder <- args[3]
if(!.file_exists(c(sample_bamfile,sample_bamindex,config_file))){
  print("Program quit : incomplete input files")
  q('no')
}
print("All input files exist.")
paste0("Load configuration file ",config_file)
source(config_file)
paste0("Done : Load configuration file ",config_file)
paste0("Load main function file ",main_function_file)
if (!.file_exists(c(main_function_file))) {
  paste0("Program quit : missing main function file ",main_function_file)
  q('no')
} else source(main_function_file)
paste0("Done Load main function file ",main_function_file)



## Load smooth function and related parameters
# source(paste0(package_dir,"/cnv-size-based-zscore/test_smoothing_zscore.R"))

if(length(args) == 4 ) {
  capture_targets_file <- args[4]
  if (!main_functions.file_exists(c(capture_targets_file))) {
    q('no')
  }
  capture_targets <- as.data.frame(read.table(gzfile(capture_targets_file),header = FALSE,
                                              sep = "\t"))[,1:3]
  colnames(capture_targets) = c("chromosome","start","end")
  capture_targets_gr=GRanges(seqnames = capture_targets$chromosome,
                             ranges = IRanges(start = as.numeric(capture_targets$start), end=as.numeric(capture_targets$end)))
}

paste0("Create blacklist genomic regions")
blacklist_targets_gr = main_functions.create_blacklist_gr(c(duke_blacklist_region,dac_blacklist_region,centromere_region))
paste0("Done : Create blacklist genomic regions")



output_file_name <- paste(tools::file_path_sans_ext(basename(sample_bamfile)),
                          "_",as.character(binsize),"k",sep="")
paste0("Set output prefix to ",output_file_name)

if(!dir.exists(output_folder)) {
  print(paste("Creating output folder ",output_folder,sep=""))
  dir.create(output_folder,recursive = TRUE)
}

paste0("Copy config file to output directory ",output_folder)
file.copy(config_file, paste(output_folder,"/",basename(config_file),sep=""), overwrite=TRUE)

sliding_windows = main_function.get_sliding_windows(binsize)


readbam_rds_file = paste0(output_folder,"/",output_file_name,"_readbam.rds")
if (main_functions.file_exists(c(readbam_rds_file)) & !rerun_readbam) {
  print(paste0("Reading existing readbam rds file",readbam_rds_file))
  reads_extracted_df = readRDS(readbam_rds_file)
} else {
  print(paste("Reading bam file",sample_bamfile, sep=" "))
  reads_extracted_df = apply(sliding_windows, 1, main_functions.readbam)
  print(paste0("Done : Reading bam file ",sample_bamfile))
  print(paste0("Writing bam rds data to ",readbam_rds_file))
  saveRDS(reads_extracted_df,readbam_rds_file)
  print(paste0("Done : Writing bam matrix data to ",readbam_rds_file))
}

runtime_var_file = paste0(output_folder,"/",runtime_variables_filename)
if (main_functions.file_exists(c(runtime_var_file))) {
  print(paste0("Overwrite existing runtime variables file ",runtime_var_file))
  file.remove(runtime_var_file)
} else print(paste0("Writing runtime variables file ",runtime_var_file))

cat(paste0("readbam_rds_file = \"",readbam_rds_file,"\"\n"),file=runtime_var_file,append = TRUE)
cat(paste0("output_folder = \"",output_folder,"\"\n"),file=runtime_var_file,append = TRUE)
