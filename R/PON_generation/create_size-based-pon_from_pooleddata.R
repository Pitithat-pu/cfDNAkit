# source("/abi/data2/puranach/ct_DNA/cfDNAkit/kit_functions.R")
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
if (length(args) < 2)
  stop("Usage: create_size-based-pon.R *runtime_variable_file(runtime_variables.R) *config_file")
runtime_var_file <- args[1]
config_file = args[2]
if(!.file_exists(c(runtime_var_file,config_file))){
  paste0("Program quit : incomplete input files")
  q('no')
}
paste0("All input files exist.")
paste0("Load configuration file ",config_file)
source(config_file)
paste0("Done : Load configuration file ",config_file)
paste0("Load main function file ",main_function_file)
if (!.file_exists(c(main_function_file))) {
  paste0("Program quit : missing main function file ",main_function_file)
  q('no')
} else source(main_function_file)
paste0("Done Load main function file ",main_function_file)
paste0("Read runtime variable file ",runtime_var_file)
source(runtime_var_file)

paste("Setting Min-Max length of Read Insert ",minimum_length," - ",maximum_length,sep="")
paste("Setting Short Length interval ",short_length_range[1],"-",short_length_range[2],sep="")
paste("Setting Long Length interval ",long_length_range[1],"-",long_length_range[2],sep="")


output_folder_temp = paste0(output_folder,"/",create_pon_outdir_prefix,"_",
                            as.character(binsize),"k")
if(!dir.exists(output_folder_temp)) {
  paste("Creating output folder ",output_folder_temp,sep="")
  dir.create(output_folder_temp,recursive = TRUE)
} else {
  paste("Output folder ",output_folder_temp," is already existing.",sep="")
  output_folder = paste(output_folder_temp,format(Sys.time(), "%Y%b%d_%H%M"),sep="_")
  paste("Creating output folder ",output_folder_temp,sep="")
  dir.create(output_folder_temp,recursive = TRUE)
}
output_folder = output_folder_temp
# else {
#   print(paste("Output directory ",output_folder," is already existing",sep = ""))
#   output_folder = paste(output_folder,format(Sys.time(), "%Y%b%d_%H%M"),sep="")
#   dir.create(output_folder,recursive = TRUE)
#   print(paste("Change output directory to ",output_folder,sep=""))
# }

# if(file.exists(paste(output_folder,"/",output_file_name,".csv",sep = ""))){
#   print(paste("Destination output file ",output_folder,"/",output_file_name,".csv is already exist",sep = ""))
#   output_file_name = paste(output_file_name,format(Sys.time(), "%Y%b%d_%H%M"),sep="")
#   print(paste("Change output filename to ",output_file_name,".csv",sep=""))
# }

.split_reads = function(all_isize){
  mod_isize = length(all_isize)%%sampling_size
  if (mod_isize > 0) {
    split(all_isize,c(rep(1:sampling_size, length(all_isize)/sampling_size),
                      1:mod_isize))
  } else split(all_isize,rep(1:sampling_size, length(all_isize)/sampling_size))
}

.get_isize_info = function(x){
  if(all(is.na(x))){
    return(c("Total Fragments"=NA,"Read Pairs in range"= NA,"Mean"=NA,"short"=NA,"long" = NA,"SLRatio" = NA))
  }
  insert_short_fragments = x[which(x >= short_length_range[1] & x <= short_length_range[2])]
  insert_long_fragments = x[which(x >= long_length_range[1] & x <= long_length_range[2])]
  c("Total Fragments"=if(mean(x) > 0) length(x) else 0,
    "Read Pairs in range"= sum(length(insert_short_fragments), length(insert_long_fragments)),
    "Mean"=mean(x),
    "short"=if(length(insert_short_fragments) > 0) length(insert_short_fragments) else 0,
    "long" = if(length(insert_long_fragments) > 0) length(insert_long_fragments) else 0,
    "SLRatio" = if(length(insert_long_fragments) > 0) {
      length(insert_short_fragments) / length(insert_long_fragments)
    } else {
      ifelse(length(insert_short_fragments) > 0,
             (length(insert_short_fragments)+1) / (length(insert_long_fragments)+1),
             NA)
    }
  )
}

.get_avg_deltaF_random = function(region_name){
  x = reads_extracted_df_control[region_name]
  passed_isize = abs(x[[1]]$isize[which(abs(x[[1]]$isize) <= maximum_length & abs(x[[1]]$isize) >= minimum_length)])
  region_coverage = length(passed_isize)
  chr=strsplit(region_name,"[:-]+")[[1]][1]
  ## trim outliner bins based on given quantile
  # q_lowerbound = chromosomal_coverage_quantile[chr,1]
  # q_upperbound = chromosomal_coverage_quantile[chr,2]
  q_lowerbound = minimum_coverage_per_bin
  q_upperbound = region_coverage
  if (q_lowerbound >  region_coverage | region_coverage > q_upperbound ){
    ## Regions with low coverage than lower bound or higher than upper bound
    sampling_read = split(rep(NA,sampling_size),cut(seq_along(rep(NA,sampling_size)),sampling_size, labels = FALSE))
    ## give list of reads with NA
  } else  {
    ## if generating control from pooled sample of normal
    passed_isize = sample(passed_isize)
    sampling_read = .split_reads(passed_isize)
    ## shuffle the reads inside the bin and split reads into sampling_size samples
  }
  reads_extracted_df_test = as.data.frame(t(sapply(sampling_read,.get_isize_info)))
}

print(paste(sep="","Deriving non-overlapping bin from QDNAseq binsize = ",binsize))
sliding_windows = main_function.get_sliding_windows(binsize)

# readbam_rds_file = paste0(output_folder,"/",output_file_name,"_readbam.rds")
output_file_name <- paste(tools::file_path_sans_ext(basename(readbam_rds_file)),
                          "_pon",sep="")
paste0("Set output prefix to ",output_file_name)

reads_extracted_df_control = readRDS(readbam_rds_file)


print(paste("Splitting reads into ",sampling_size," control samples and getting short read ratio",sep = ""))
all_region_name = names(reads_extracted_df_control)
control_regions_df = as.data.frame(t(sapply(all_region_name,.get_avg_deltaF_random) ))
print(paste("Done : Splitting reads into ",sampling_size," control samples and getting short read ratio",sep = ""))
print("Getting DeltaF")


SLRatio_df = as.data.frame(do.call(rbind,control_regions_df$`SLRatio`))
totalread_df = as.data.frame(do.call(rbind,control_regions_df$`Total Fragments`))
readpair_df = as.data.frame(do.call(rbind,control_regions_df$`Read Pairs in range`))
shortread_df = as.data.frame(do.call(rbind,control_regions_df$short))
longread_df = as.data.frame(do.call(rbind,control_regions_df$long))

print("Correct percent GC-bias.")
# all_region_name_non_NA = rownames(totalread_df[which(complete.cases(totalread_df)),])
all_region_name_non_NA = rownames(totalread_df)[which(complete.cases(totalread_df))]
gc_contents = sliding_windows[all_region_name_non_NA,]$gc/100
totalread.corrected_df = sapply(totalread_df[all_region_name_non_NA,], function(x){
  main_functions.bias_correct(x,gc_contents)
})
rownames(totalread.corrected_df) = all_region_name_non_NA
readpair.corrected_df = sapply(readpair_df[all_region_name_non_NA,], function(x){
  main_functions.bias_correct(x,gc_contents)
})
rownames(readpair.corrected_df) = all_region_name_non_NA
shortread.corrected_df = sapply(shortread_df[all_region_name_non_NA,], function(x){
  main_functions.bias_correct(x,gc_contents)
})
rownames(shortread.corrected_df) = all_region_name_non_NA
longread.corrected_df = sapply(longread_df[all_region_name_non_NA,], function(x){
  main_functions.bias_correct(x,gc_contents)
})
rownames(longread.corrected_df) = all_region_name_non_NA

temp_shortread.corrected_df = shortread_df
temp_shortread.corrected_df[all_region_name_non_NA,] = shortread.corrected_df
shortread.corrected_df=temp_shortread.corrected_df
temp_longread.corrected_df = shortread_df
temp_longread.corrected_df[all_region_name_non_NA,] = longread.corrected_df
longread.corrected_df = temp_longread.corrected_df
temp_totalread.corrected_df = totalread_df
temp_totalread.corrected_df[all_region_name_non_NA,] = totalread.corrected_df
totalread.corrected_df = temp_totalread.corrected_df
temp_readpair.corrected_df = readpair_df
temp_readpair.corrected_df[all_region_name_non_NA,] = readpair.corrected_df
readpair.corrected_df = temp_totalread.corrected_df


SLRatio.corrected_df = shortread.corrected_df / longread.corrected_df
colnames(SLRatio.corrected_df) = rep(1:sampling_size)
print("Done : correct percent GC-bias.")
print(paste("Writing result to ",output_folder,"/",output_file_name,"_shortread.corrected.csv", sep=""))
write.table(shortread.corrected_df, file = paste(output_folder,"/",output_file_name,"_shortread.corrected.csv", sep=""), sep = "\t")
print(paste("Done : Writing result to ",output_folder,"/",output_file_name,"_shortread.corrected.csv", sep=""))
print(paste("Writing result to ",output_folder,"/",output_file_name,"_longread.corrected.csv", sep=""))
write.table(longread.corrected_df, file = paste(output_folder,"/",output_file_name,"_longread.corrected.csv", sep=""), sep = "\t")
print(paste("Done : Writing result to ",output_folder,"/",output_file_name,"_longread.corrected.csv", sep=""))


  # reference_short_read = colSums(na.rm = TRUE,shortread.corrected_df[which(startsWith(rownames(shortread.corrected_df),"X")),])
  # reference_long_read = colSums(na.rm = TRUE,longread.corrected_df[which(startsWith(rownames(longread.corrected_df),"X")),])
reference_short_read = colSums(na.rm = TRUE,shortread.corrected_df)
reference_long_read = colSums(na.rm = TRUE,longread.corrected_df)

P150_reference =  reference_short_read / reference_long_read
delta_f_df.corrected = as.data.frame(t(sapply(all_region_name,function(region_name){
  if (all(is.na(control_regions_df$`Read Pairs in range`[region_name][[1]])==TRUE)) {
    Delta_F = rep(NA,sampling_size)
  } else {
    SLRatio.corrected = SLRatio.corrected_df[region_name,]
    Delta_F = unlist(SLRatio.corrected - P150_reference)
  }

}) ))

reference_long_read = colSums(na.rm = TRUE,longread_df[which(startsWith(rownames(longread_df),"X")),])
reference_short_read = colSums(na.rm = TRUE,shortread_df[which(startsWith(rownames(shortread_df),"X")),])
P150_reference =  reference_short_read / reference_long_read
delta_f_df = as.data.frame(t(sapply(all_region_name,function(region_name){
  if (length(control_regions_df$`Read Pairs in range`[region_name][[1]])==0) {
    Delta_F = rep(NA,sampling_size)
  } else {
    SLRatio = SLRatio_df[region_name,]
    Delta_F = unlist(SLRatio - P150_reference)
  }
}) ))
print("Done : Getting DeltaF")

print("Calculate z-score for each sample.")
all_region_name = rownames(delta_f_df.corrected)
control_zscore_df = as.data.frame(sapply(delta_f_df.corrected, function(sample_delta_f){
  sample_delta_f = do.call(rbind.data.frame, as.list(sample_delta_f))
  rownames(sample_delta_f) = all_region_name
  colnames(sample_delta_f) = c("Delta_F")
  zscore_matrix = main_functions.get_zscore_per_binsize(sample_delta_f,delta_f_df.corrected,binsize)
  zscore_matrix$size_based_zscore
}))
rownames(control_zscore_df) = all_region_name
print(paste("Writing result to ",output_folder,"/",output_file_name,"_zscore.csv", sep=""))
write.table(x = control_zscore_df, file = paste(output_folder,"/",output_file_name,"_zscore.csv", sep=""), sep = "\t")
print(paste("Done : Writing result to ",output_folder,"/",output_file_name,"_zscore.csv", sep=""))

control_zscore_df_non_NA = control_zscore_df[complete.cases(control_zscore_df),]
control_genomewide_zscores = sapply(control_zscore_df_non_NA, function(sample_zscore){
  main_functions.calculate_gw_zscore(sample_zscore,control_zscore_df_non_NA)
})

print(paste("Writing result to ",output_folder,"/",output_file_name,"_genomewide_zscores.csv", sep=""))
write.table(x = control_genomewide_zscores, file = paste(output_folder,"/",output_file_name,"_genomewide_zscores.csv", sep=""), sep = "\t")
print(paste("Done : Writing result to ",output_folder,"/",output_file_name,"_genomewide_zscores.csv", sep=""))

# delta_f_df_correlation = as.data.frame(sapply(all_region_name,function(region_name){
#   cor(use = "pairwise.complete.obs",sapply(delta_f_df[region_name,],function(x) as.numeric(as.character(x)))
#       , sapply(delta_f_df.corrected[region_name,],function(x) as.numeric(as.character(x))))
# }))
# colnames(delta_f_df_correlation) = c("GC_uncorrected_vs_corrected_correlation")
# print(paste("Writing deltaF correlation result to ",output_folder,"/",output_file_name,"_deltaf_cor.csv", sep=""))
# write.table(x = delta_f_df_correlation, file = paste(output_folder,"/",output_file_name,"_deltaf_cor.csv", sep=""), sep = "\t")
# print(paste("Done : Writing delta F correlation result to ",output_folder,"/",output_file_name,"_deltaf_cor.csv", sep=""))

colnames(delta_f_df.corrected) = rep(1:sampling_size)
colnames(delta_f_df) = rep(1:sampling_size)
print(paste("Writing result to ",output_folder,"/",output_file_name,".csv", sep=""))
write.table(x = delta_f_df.corrected, file = paste(output_folder,"/",output_file_name,".csv", sep=""), sep = "\t")
print(paste("Done : Writing result to ",output_folder,"/",output_file_name,".csv", sep=""))


control_regions_df_towrite = as.data.frame(apply(t(control_regions_df), 1, function(x){
  unlist(x,use.names = FALSE)
}))
control_regions_df_towrite$`Total Fragments_corrected` = c(sapply(totalread.corrected_df, function(x){
   x
}))
control_regions_df_towrite$`Read Pairs in range_corrected` = c(sapply(readpair.corrected_df, function(x){
  x
}))
control_regions_df_towrite$`short_corrected` = c(sapply(shortread.corrected_df, function(x){
  x
}))
control_regions_df_towrite$`long_corrected` = c(sapply(longread.corrected_df, function(x){
  x
}) )
control_regions_df_towrite$`SLRatio_corrected` = c(sapply(SLRatio.corrected_df, function(x){
  x
}))
control_regions_df_towrite$chr_region = c(rep(rownames(control_regions_df),sampling_size))
control_regions_df_towrite$sampling_id = rep(1:sampling_size, each=length(rownames(control_regions_df)))

print(paste("Writing regions_df result to ",output_folder,"/",output_file_name,"_regions_df.csv", sep=""))
write.table(x = control_regions_df_towrite, file = paste(output_folder,"/",output_file_name,"_regions_df.csv", sep=""), sep = "\t")
print(paste("Done : Writing regions_df result to ",output_folder,"/",output_file_name,"_regions_df.csv", sep=""))
