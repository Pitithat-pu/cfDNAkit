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
  stop("Usage: initial_zscore_segmentation.R *runtime_variable_file(runtime_variables.R) *config_file")
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


## Load smooth function and related parameters
source(paste0(package_dir,"/cnv-size-based-zscore/test_smoothing_zscore.R"))


paste0("Begin Sized-based Fragment Analysis")
paste("Setting Min-Max length of Read Insert ",minimum_length," - ",maximum_length,sep="")
paste("Setting Short Length interval ",short_length_range[1],"-",short_length_range[2],sep="")
paste("Setting Long Length interval ",long_length_range[1],"-",long_length_range[2],sep="")
paste("Setting z-score cutoff [",cutoff_zscore[2],",",cutoff_zscore[1],"]",sep="")
# paste0("Create blacklist genomic regions")
# blacklist_targets_gr = main_functions.create_blacklist_gr(c(duke_blacklist_region,dac_blacklist_region))
# paste0("Done : Create blacklist genomic regions")



# output_file_name <- paste(tools::file_path_sans_ext(basename(sample_bamfile)),
#                           "_initzscore_",as.character(binsize),sep="")
# paste0("Set output prefix to ",output_file_name)
output_folder_temp = paste0(output_folder,"/",init_zscore_outdir_prefix,"_",
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

# if(!dir.exists(paste0(output_folder,"/",init_zscore_outdir_prefix,"_",
#                       as.character(binsize),"k"))) {
#   paste("Creating output folder ",output_folder,sep="")
#   dir.create(output_folder,recursive = TRUE)
# } else {
#   paste("Output folder ",output_folder," is already existing.",sep="")
#   output_folder = paste(output_folder,format(Sys.time(), "%Y%b%d_%H%M"),sep="_")
#   paste("Creating output folder ",output_folder,sep="")
#   dir.create(output_folder,recursive = TRUE)
# }
# file.copy(config_file, paste(output_folder,"/",basename(config_file),sep=""), overwrite=TRUE)
#
sliding_windows = main_function.get_sliding_windows(binsize)
#
# print(paste("Reading bam file",sample_bamfile, sep=" "))
# reads_extracted_df = apply(sliding_windows, 1, readbam)
# print(paste("Reading bam file",sample_bamfile,"(Done)", sep=" "))
# saveRDS(reads_extracted_df, paste(output_folder,"/",output_file_name,"_",binsize,"k_readbam.rds",sep=""))
output_file_name <- paste(tools::file_path_sans_ext(basename(readbam_rds_file)),
                          "_initzscore",sep="")
paste0("Set output prefix to ",output_file_name)

reads_extracted_df = readRDS(readbam_rds_file)
# chromosomal_coverage_quantile = main_functions.get_chromosomal_coverage_quantile(sliding_windows,reads_extracted_df,
#                                                                   fragment_perwindow_lower_bound_percent,
#                                                                   fragment_perwindow_upper_bound_percent)


paste0("Getting Test Regions info")
reads_extracted_df_info = as.data.frame(t(sapply(reads_extracted_df, function(x){
  passed_isize = abs(x$isize[which(abs(x$isize) <= maximum_length & abs(x$isize) >= minimum_length)])
  region_coverage = length(passed_isize)
  chr=x$rname[1]
  # q_lowerbound = chromosomal_coverage_quantile[chr,1]
  # q_upperbound = chromosomal_coverage_quantile[chr,2]
  q_lowerbound = minimum_coverage_per_bin
  q_upperbound = max(region_coverage)
  trimmed = if (region_coverage>0 & q_lowerbound <=  region_coverage &  region_coverage <= q_upperbound ) FALSE else TRUE
  insert_short_fragments = passed_isize[which(passed_isize >= short_length_range[1] & passed_isize <= short_length_range[2])]
  insert_long_fragments = passed_isize[which(passed_isize >= long_length_range[1] & passed_isize <= long_length_range[2])]
  c("Total Fragments"=length(passed_isize),
    "Read Pairs in range"= sum(length(insert_short_fragments), length(insert_long_fragments)),
    "Mean"=mean(passed_isize),
    "short"=if(length(insert_short_fragments) > 0) length(insert_short_fragments) else 0,
    "long" = if(length(insert_long_fragments) > 0) length(insert_long_fragments) else 0,
    "SLRatio" = if(length(insert_long_fragments) > 0) {
      length(insert_short_fragments) / length(insert_long_fragments)
    } else {
      ifelse(length(insert_short_fragments) > 0,
             (length(insert_short_fragments)+1) / (length(insert_long_fragments)+1),
              NA)
    },  "Trimmed" = trimmed)
})))

coverage_df = data.frame(`Total Fragments` = reads_extracted_df_info$`Total Fragments`,Trimming = "Before Trimming")

coverage_df = rbind(coverage_df,data.frame(`Total Fragments` = reads_extracted_df_info$`Total Fragments`[which(reads_extracted_df_info$Trimmed==FALSE)],
                                           Trimming = "After Trimming"))


paste("Plot per-windows fragments ",output_folder,"/",output_file_name,"_fragment_per_windows_plots.png", sep="")
png(paste(output_folder,"/",output_file_name,"_fragment_per_windows_plots.png", sep=""), width=750, height=550, units="px")
ggplot(coverage_df, aes(Trimming,Total.Fragments)) +
  ggtitle("Number of Fragments per window Before/After Trimming ") +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1)  + coord_flip() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
dev.off()
paste("Plot per-windows fragments ",output_folder,"/",output_file_name,"_fragment_per_windows_plots.png (Done)", sep="")

## get regions with enough coverage
# all_region_name = rownames(reads_extracted_df_info[which(!is.na(reads_extracted_df_info$`Total Fragments`)),])
reads_extracted_df_info[which(reads_extracted_df_info$Trimmed==TRUE),]=NA
all_region_name_row = which(complete.cases(reads_extracted_df_info))
all_region_name = rownames(reads_extracted_df_info[all_region_name_row,])
SLRatio_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$SLRatio[all_region_name_row])))
totalread_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$`Total Fragments`[all_region_name_row])))
readpair_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$`Read Pairs in range`[all_region_name_row])))
shortread_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$short[all_region_name_row])))
longread_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$long[all_region_name_row])))


# totalread_df=as.data.frame(sapply(all_region_name,function(region_name){
#   reads_extracted_df_info[region_name,]$`Total Fragments`
# }))
# readpair_df=as.data.frame(sapply(all_region_name,function(region_name){
#   reads_extracted_df_info[region_name,]$`Read Pairs in range`
# }))
# shortread_df=as.data.frame(sapply(all_region_name,function(region_name){
#   reads_extracted_df_info[region_name,]$short
# }))
# longread_df=as.data.frame(sapply(all_region_name,function(region_name){
#   reads_extracted_df_info[region_name,]$long
# }))
# SLRatio = reads_extracted_df_info[all_region_name,]$SLRatio

totalread_df$corrected = main_functions.bias_correct(totalread_df[,1],sliding_windows[all_region_name,]$gc/100)
readpair_df$corrected = main_functions.bias_correct(readpair_df[,1],sliding_windows[all_region_name,]$gc/100)
shortread_df$corrected = main_functions.bias_correct(shortread_df[,1],sliding_windows[all_region_name,]$gc/100)
longread_df$corrected = main_functions.bias_correct(longread_df[,1],sliding_windows[all_region_name,]$gc/100)
# shortread_df$corrected = readpair_df$corrected * SLRatio
# longread_df$corrected = readpair_df$corrected * (1 - SLRatio)
SLRatio.corrected = shortread_df$corrected / longread_df$corrected


reads_extracted_df_info$`Total Fragments_corrected`= NA
reads_extracted_df_info$`Read Pairs in range_corrected`= NA
reads_extracted_df_info$`short_corrected` = NA
reads_extracted_df_info$`long_corrected` = NA
reads_extracted_df_info$`SLRatio_corrected`= NA

reads_extracted_df_info[all_region_name,]$`Total Fragments_corrected` = totalread_df$corrected
reads_extracted_df_info[all_region_name,]$`Read Pairs in range_corrected` = readpair_df$corrected
reads_extracted_df_info[all_region_name,]$`short_corrected` = shortread_df$corrected
reads_extracted_df_info[all_region_name,]$`long_corrected` = longread_df$corrected
reads_extracted_df_info[all_region_name,]$`SLRatio_corrected` = SLRatio.corrected
reads_extracted_df_info$gc=sliding_windows[rownames(reads_extracted_df_info),]$gc
paste("Writing reads_extracted_info to ",output_folder,"/",output_file_name,"_reads_extracted_info.csv",sep = "")
write.table(reads_extracted_df_info, file = paste(output_folder,"/",output_file_name,"_reads_extracted_info.csv", sep=""),
            row.names = TRUE,col.names = TRUE)
paste("Writing reads_extracted_info to ",output_folder,"/",output_file_name,"_reads_extracted_info.csv (Done)",sep = "")



paste("Plot qc-bias correction to ",output_folder,"/",output_file_name,"_bias.png",sep="")
png(paste(output_folder,"/",output_file_name,"_bias.png", sep=""), width=1800, height=1350, units="px")
gc_coverage_plot = ggplot(reads_extracted_df_info) + aes(gc,`Total Fragments`) +
  geom_point() + geom_smooth()
gc_coverage_corrected_plot = ggplot(reads_extracted_df_info) + aes(gc,`Total Fragments_corrected`) +
  geom_point() + geom_smooth()
gc_SLRatio_plot = ggplot(reads_extracted_df_info) + aes(gc,SLRatio) +
  geom_point() + geom_smooth() + ylim(0,0.8)
gc_SLRatio_corrected_plot = ggplot(reads_extracted_df_info) + aes(gc,SLRatio_corrected) +
  geom_point() + geom_smooth() + ylim(0,0.8)
grid.arrange(gc_coverage_plot, gc_coverage_corrected_plot,
             gc_SLRatio_plot,gc_SLRatio_corrected_plot, ncol=2)
dev.off()
# plot(sliding_windows$mappability,reads_extracted_df_info$`Total Fragments`,pch=16)
# plot(sliding_windows$mappability,reads_extracted_df_info$`Total Fragments_corrected`,pch=16)
# plot(sliding_windows$mappability,reads_extracted_df_info$SLRatio,pch=16)
# plot(sliding_windows$mappability,reads_extracted_df_info$SLRatio_corrected,ylim=c(0,0.6),pch=16)



test_regions_df_corrected = main_functions.get_short_ratio_corrected(reads_extracted_df_info, rownames(sliding_windows))


control_shortread_df = read.table(control_shortread_file,
                                  header = TRUE,stringsAsFactors = FALSE)
control_longread_df = read.table(control_longread_file,
                                 header = TRUE,stringsAsFactors = FALSE)

all_region_name = rownames(control_longread_df)
control_deltaF_df = as.data.frame(t(sapply(all_region_name,function(test_region_name){
  reference_shortread_df = control_shortread_df[!rownames(control_shortread_df) %in% test_region_name,]
  reference_longread_df = control_longread_df[!rownames(control_longread_df) %in% test_region_name,]
  SLRatio_corrected = control_shortread_df[test_region_name,] / control_longread_df[test_region_name,]
  SLRatio_reference = colSums(reference_shortread_df, na.rm = TRUE) / colSums(reference_longread_df, na.rm = TRUE)
  Delta_F= main_functions.calculate_deltaF(SLRatio_corrected, SLRatio_reference)
})))

# control_SLRatio_df = control_shortread_df / control_longread_df
# all_region_name = rownames(control_SLRatio_df)
# get_reference_region <- function(test_regions_df_corrected){
#   test_regions_df_reference =  test_regions_df_corrected[which(startsWith(rownames(test_regions_df_corrected),"X")),]
#   test_regions_df_reference =  test_regions_df_corrected[which(startsWith(rownames(test_regions_df_corrected),"X")),]
#   return (rownames(test_regions_df_reference))
# }
# reference_regions = get_reference_region(test_regions_df_corrected)
# control_shortread_reference_df = control_shortread_df[reference_regions,]
# control_longread_reference_df = control_longread_df[reference_regions,]
# reference_short_read = colSums(control_shortread_reference_df,na.rm = TRUE)
# reference_long_read = colSums(control_longread_reference_df,na.rm = TRUE)
# P150_reference =  reference_short_read / reference_long_read
#
# control_deltaF_df = as.data.frame(t(sapply(all_region_name,function(region_name){
#   SLRatio.corrected = control_SLRatio_df[region_name,]
#   Delta_F = as.numeric(SLRatio.corrected - P150_reference)
# }) ))

chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
chrLength_df = chrLength_df[which(chrLength_df$V1!="Y"),]
colnames(chrLength_df) = c("Chromosome", "Length")
chrNames = chrLength_df$Chromosome
chrLength = chrLength_df$Length
names(chrLength) = chrNames
chroffsets = cumsum(as.numeric(chrLength))
chroffsets <- c(0, chroffsets[0:(length(chroffsets)-1)])
names(chroffsets) <- names(chrLength)
chrMids <- cumsum(as.numeric(chrLength))
chrMids <- (chrMids + chroffsets)/2
names(chrMids) <- names(chrLength)

###
z_score_df_no_NaN = main_functions.get_zscore_per_binsize(test_regions_df_corrected,control_deltaF_df,binsize)
z_score_df_no_NaN = z_score_df_no_NaN[which(z_score_df_no_NaN$size_based_zscore != 'NaN' ),]
z_score_df_no_NaN$scaledPos = (z_score_df_no_NaN$Start + z_score_df_no_NaN$End)/2 +
  as.numeric(chroffsets[z_score_df_no_NaN$Chromosome])

# z_score_df_no_NaN_1000 = main_functions.get_zscore_per_binsize(test_regions_df_corrected,control_deltaF_df,1000)
# z_score_df_no_NaN_1000 = z_score_df_no_NaN_1000[which(z_score_df_no_NaN_1000$size_based_zscore != 'NaN' ),]
# z_score_df_no_NaN_1000$scaledPos = (z_score_df_no_NaN_1000$Start + z_score_df_no_NaN_1000$End)/2 +
#   as.numeric(chroffsets[z_score_df_no_NaN_1000$Chromosome])
#
# z_score_df_no_NaN_5000 = main_functions.get_zscore_per_binsize(test_regions_df_corrected,control_deltaF_df,5000)
# z_score_df_no_NaN_5000 = z_score_df_no_NaN_5000[which(z_score_df_no_NaN_5000$size_based_zscore != 'NaN' ),]
# z_score_df_no_NaN_5000$scaledPos = (z_score_df_no_NaN_5000$Start + z_score_df_no_NaN_5000$End)/2 +
#   as.numeric(chroffsets[z_score_df_no_NaN_5000$Chromosome])

z_score_df_no_NaN = test.smoothing_zscore.smooth_zscore(z_score_df_no_NaN)
# z_score_df_no_NaN_1000 = test.smoothing_zscore.smooth_zscore(z_score_df_no_NaN_1000)


z_score_df_no_NaN_file = paste0(output_folder,"/",output_file_name,"_z_score_",binsize,".csv")
write.table(x = z_score_df_no_NaN,
            file = z_score_df_no_NaN_file,
            quote = FALSE,row.names = FALSE,sep = "\t")
# write.table(x = z_score_df_no_NaN_1000,
#             file = paste0(output_folder,"/",output_file_name,"_z_score_1000.csv"),
#             quote = FALSE,row.names = FALSE,sep = "\t")
# write.table(x = z_score_df_no_NaN_5000,
#             file = paste0(output_folder,"/",output_file_name,"_z_score_5000.csv"),
#             quote = FALSE,row.names = FALSE,sep = "\t")

cat(paste0("z_score_df_file = \"",z_score_df_no_NaN_file,"\"\n"),file=runtime_var_file,append = TRUE)




chromosomes <- setNames(1:length(unique(z_score_df_no_NaN$Chromosome)),
                        unique(z_score_df_no_NaN$Chromosome))

paste("Plotting result to ",output_folder,"/",output_file_name,"_plots.png", sep="")
png(paste(output_folder,"/",output_file_name,"_plots.png", sep=""), width=1800, height=600, units="px")
par(mfrow=c(1,1))
main_functions.plot_zscore(z_score_df_no_NaN,binsize)
# main_functions.plot_zscore(z_score_df_no_NaN_1000,1000)
# main_functions.plot_zscore(z_score_df_no_NaN_5000,5000)
dev.off()
paste("Plotting result to ",output_folder,"/",output_file_name,"_plots.png (Done)", sep="")

paste("Plotting result to ",output_folder,"/",output_file_name,"_smoothedplots.png", sep="")
png(paste(output_folder,"/",output_file_name,"_smoothedplots.png", sep=""), width=1800, height=600, units="px")
par(mfrow=c(1,1))
main_functions.plot_smoothed_zscore(z_score_df_no_NaN,binsize)
# main_functions.plot_smoothed_zscore(z_score_df_no_NaN_1000,1000)
dev.off()
paste("Plotting result to ",output_folder,"/",output_file_name,"_smoothedplots.png (Done)", sep="")

paste0("Sized-based Fragment Analysis Completed")



region_names_non_NA = rownames(z_score_df_no_NaN)
if (calculate_gwzscore) {
  control_zscore_df = read.table(control_zscore_file,header = TRUE,stringsAsFactors = FALSE)
  genomewide_zscore = main_functions.calculate_gw_zscore(z_score_df_no_NaN$size_based_zscore,control_zscore_df[region_names_non_NA,])
}
if (calculate_mad) {
  mad_score = mad(z_score_df_no_NaN$size_based_zscore)
  mad_smoothed_score = mad(z_score_df_no_NaN$smoothed_zscore,na.rm = TRUE)
}

# genomewide_smoothed_zscore = calculate_gw_zscore(z_score_df_no_NaN$smoothed_zscore,control_zscore_df[region_names_non_NA,])

genomewide_scores = data.frame(Samfile=basename(readbam_rds_file),
                               MAD_Score=ifelse(calculate_mad,mad_score,NA),
                               MAD_Smoothed_Score=ifelse(calculate_mad,mad_smoothed_score,NA),
                               GW_ZScore=ifelse(calculate_gwzscore,genomewide_zscore,NA)
                               # ,GW_Smoothed_ZScore=genomewide_smoothed_zscore
                               )
paste("Writing genomewide_scores to ",output_folder,"/",output_file_name,"_genomewide_scores.csv",sep = "")
write.table(genomewide_scores,file = paste(output_folder,"/",output_file_name,"_genomewide_scores.csv", sep=""),
            row.names = FALSE,col.names = TRUE)
paste("Writing genomewide_scores to ",output_folder,"/",output_file_name,"_genomewide_scores.csv (Done)",sep = "")


