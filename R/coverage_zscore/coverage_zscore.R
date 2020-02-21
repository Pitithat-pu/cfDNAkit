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
  stop("Usage: coverage_zscore *runtime_variable_file(runtime_variables.R) *config_file")
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
# paste("Setting Min-Max length of Read Insert ",minimum_length," - ",maximum_length,sep="")
# paste("Setting Short Length interval ",short_length_range[1],"-",short_length_range[2],sep="")
# paste("Setting Long Length interval ",long_length_range[1],"-",long_length_range[2],sep="")
# paste("Setting z-score cutoff [",cutoff_zscore[2],",",cutoff_zscore[1],"]",sep="")

output_folder_temp = paste0(output_folder,"/",coverage_zscore_outdir_prefix,"_",
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

sliding_windows = main_function.get_sliding_windows(binsize)
output_file_name <- paste(tools::file_path_sans_ext(basename(readbam_rds_file)),
                          "_coverage_zscore",sep="")
paste0("Set output prefix to ",output_file_name)

paste0("Reading RDS file ",readbam_rds_file)
reads_extracted_df = readRDS(readbam_rds_file)

reads_extracted_df_info = as.data.frame(t(sapply(reads_extracted_df, function(x){
  ## include all fragment size
  passed_isize = abs(x$isize)
  region_coverage = length(passed_isize)
  # all_isize = x$isize
  # region_coverage = length(all_isize)
  chr=x$rname[1]
  # q_lowerbound = chromosomal_coverage_quantile[chr,1]
  # q_upperbound = chromosomal_coverage_quantile[chr,2]
  q_lowerbound = ifelse(chr!="X",minimum_coverage_per_bin,minimum_coverage_per_bin/2)
  q_upperbound = max(region_coverage)
  trimmed = if (region_coverage>0 & q_lowerbound <=  region_coverage &  region_coverage <= q_upperbound ) FALSE else TRUE
  c("Total Fragments"=region_coverage
    ,"Trimmed" = trimmed)
})))

# coverage_df = data.frame(`Total Fragments` = reads_extracted_df_info$`Total Fragments`,Trimming = "Before Trimming")
#
# coverage_df = rbind(coverage_df,data.frame(`Total Fragments` = reads_extracted_df_info$`Total Fragments`[which(reads_extracted_df_info$Trimmed==FALSE)],
#                                            Trimming = "After Trimming"))
reads_extracted_df_info[which(reads_extracted_df_info$Trimmed==TRUE),]=NA
all_region_name_row = which(complete.cases(reads_extracted_df_info))
all_region_name = rownames(reads_extracted_df_info[all_region_name_row,])

totalread_df = as.data.frame(do.call(rbind,as.list(reads_extracted_df_info$`Total Fragments`[all_region_name_row])))
rownames(totalread_df) = all_region_name
totalread_df$corrected = main_functions.bias_correct(totalread_df[,1],sliding_windows[all_region_name,]$gc/100)
totalread_df = totalread_df[all_region_name[-grep('^X:', all_region_name)],]
totalread_df$scaled_count = (totalread_df$corrected - mean(totalread_df$corrected,na.rm=TRUE)) / sd(totalread_df$corrected,na.rm=TRUE)

all_region_name = all_region_name[-grep('^X:', all_region_name)]
reads_extracted_df_info$`Total Fragments_corrected`= NA
reads_extracted_df_info = reads_extracted_df_info[all_region_name,]
reads_extracted_df_info[all_region_name,]$`Total Fragments_corrected` = totalread_df$corrected
reads_extracted_df_info$gc=sliding_windows[rownames(reads_extracted_df_info),]$gc
reads_extracted_df_info$Scaled_Fragment_Count = NA
reads_extracted_df_info[all_region_name,]$Scaled_Fragment_Count = totalread_df$scaled_count


paste("Plot qc-bias correction to ",output_folder,"/",output_file_name,"_bias.png",sep="")
png(paste(output_folder,"/",output_file_name,"_bias.png", sep=""), width=1800, height=1350, units="px")
gc_coverage_plot = ggplot(reads_extracted_df_info) + aes(gc,`Total Fragments`) +
  geom_point() + geom_smooth()
gc_coverage_corrected_plot = ggplot(reads_extracted_df_info) + aes(gc,`Total Fragments_corrected`) +
  geom_point() + geom_smooth()
grid.arrange(gc_coverage_plot, gc_coverage_corrected_plot, ncol=2)
dev.off()

if(main_functions.file_exists(output_regions_bed)) {
  output_regions = read.delim(output_regions_bed,header=FALSE,stringsAsFactors=FALSE)
  output_regions_gr = GRanges(seqnames = output_regions[,1],
                              ranges = IRanges(start = output_regions[,2],
                                               end=output_regions[,3]))
  output_regions_df = as.data.frame(output_regions_gr)
}

chr_position = strsplit(rownames(reads_extracted_df_info),"[:-]+")
chr = unlist(lapply(chr_position, '[[', 1))
start = unlist(lapply(chr_position, '[[', 2))
end = unlist(lapply(chr_position, '[[', 3))
sliding_windows_gr <- GRanges(seqnames = chr,
                              ranges = IRanges(start=as.numeric(start), end=as.numeric(end)))
overlapping_regions_gr = subsetByOverlaps(sliding_windows_gr,output_regions_gr)
overlapping_regions_df = mutate(as.data.frame(overlapping_regions_gr),
                                region_name=paste0(seqnames,":",start,"-",end))
overlapping_hits = findOverlaps(overlapping_regions_gr,output_regions_gr)
reads_extracted_df_info = reads_extracted_df_info[overlapping_regions_df$region_name,]
reads_extracted_df_info$output_region_row = overlapping_hits@to

reads_extracted_df_info = mutate(reads_extracted_df_info, output_region_row = paste0(output_regions_df[overlapping_hits@to,"seqnames"],
                                                                  ":",output_regions_df[overlapping_hits@to,"start"],
                                                                  "-",output_regions_df[overlapping_hits@to,"end"]))
merged_count_df = reads_extracted_df_info %>% group_by(output_region_row) %>% summarise(total_fragment=sum(`Total Fragments`,na.rm = TRUE),
                                                                             total_fragment_corrected = sum(`Total Fragments_corrected`,na.rm = TRUE)                                                                             )
merged_count_df$scaled_count = (merged_count_df$total_fragment_corrected -
                                  mean(merged_count_df$total_fragment_corrected,na.rm=TRUE)) / sd(merged_count_df$total_fragment_corrected,na.rm=TRUE)
merged_count_df=tibble::column_to_rownames(merged_count_df,var = "output_region_row")

paste("Writing reads_extracted_info to ",output_folder,"/",output_file_name,"_reads_extracted_info.csv",sep = "")
write.table(reads_extracted_df_info, file = paste(output_folder,"/",output_file_name,"_reads_extracted_info.csv", sep=""),
            row.names = TRUE,col.names = TRUE)
paste("Writing reads_extracted_info to ",output_folder,"/",output_file_name,"_reads_extracted_info.csv (Done)",sep = "")

control_total_fragment_df = read.table(control_total_fragment_file,
                                  header = TRUE,stringsAsFactors = FALSE)
totalread_df = totalread_df[all_region_name[-grep('^X:', all_region_name)],]
chr_position = strsplit(rownames(control_total_fragment_df),"[:-]+")
chr = unlist(lapply(chr_position, '[[', 1))
start = unlist(lapply(chr_position, '[[', 2))
end = unlist(lapply(chr_position, '[[', 3))
sliding_windows_gr <- GRanges(seqnames = chr,
                              ranges = IRanges(start=as.numeric(start), end=as.numeric(end)))
overlapping_regions_gr = subsetByOverlaps(sliding_windows_gr,output_regions_gr)
overlapping_regions_df = mutate(as.data.frame(overlapping_regions_gr),
                                region_name=paste0(seqnames,":",start,"-",end))
overlapping_hits = findOverlaps(overlapping_regions_gr,output_regions_gr)
control_total_fragment_df = control_total_fragment_df[overlapping_regions_df$region_name,]
control_total_fragment_df = mutate(control_total_fragment_df, output_region_row = paste0(output_regions_df[overlapping_hits@to,"seqnames"],
                                                                                     ":",output_regions_df[overlapping_hits@to,"start"],
                                                                                     "-",output_regions_df[overlapping_hits@to,"end"]))
merged_controls_count_df = control_total_fragment_df %>% group_by(output_region_row) %>% summarise_all(sum,na.rm = TRUE)
merged_controls_count_df = dplyr::filter(merged_controls_count_df, !grepl('^X', output_region_row))
merged_controls_count_df = tibble::column_to_rownames(merged_controls_count_df,var = "output_region_row")

control_scaled_fragment_count_df = sapply(merged_controls_count_df, function(fragment_counts){
  (fragment_counts-mean(fragment_counts,na.rm=TRUE))/
    sd(fragment_counts,na.rm=TRUE)})
rownames(control_scaled_fragment_count_df)=rownames(merged_controls_count_df)


chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
chrLength_df = chrLength_df[which(!chrLength_df$V1 %in% c("Y","X")),]
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


.coverage_zscore = function(merged_count_df,control_scaled_fragment_count_df){
  zscore_df = as.data.frame(t(sapply(rownames(merged_count_df), function(x){
    control_count = as.numeric(unlist(control_scaled_fragment_count_df[which(rownames(control_scaled_fragment_count_df)==x),]))
    chr_position = strsplit(x,"[:-]+")[[1]]
    chrN_ref_total_fragment_mean = mean(control_count,na.rm = TRUE)
    chrN_ref_total_fragment_sd = sd(control_count,na.rm = TRUE)
    chrN_sample_total_fragment = merged_count_df[which(rownames(merged_count_df)==x),
                                                         "scaled_count"]
    coverage_zscore = (chrN_sample_total_fragment - chrN_ref_total_fragment_mean) / chrN_ref_total_fragment_sd

    c("Chromosome" = chr_position[1],
      "Start" = chr_position[2],
      "End" = chr_position[3],
      "chrN_ref_total_fragment_mean"= chrN_ref_total_fragment_mean,
      "chrN_ref_total_fragment_sd" = chrN_ref_total_fragment_sd,
      "chrN_sample_total_fragment" = chrN_sample_total_fragment,
      "coverage_zscore" = coverage_zscore
    )
  })))
}



# .coverage_zscore = function(reads_extracted_df_info,control_scaled_fragment_count_df,binsize){
#   zscore_df = as.data.frame(t(sapply(rownames(reads_extracted_df_info), function(x){
#     region_deltaf = unlist(control_scaled_fragment_count_df[x,])
#     chr_position = strsplit(x,"[:-]+")[[1]]
#     chrN_ref_total_fragment_mean = mean(region_deltaf,na.rm = TRUE)
#     chrN_ref_total_fragment_sd = sd(region_deltaf,na.rm = TRUE)
#     chrN_sample_total_fragment = reads_extracted_df_info[which(rownames(reads_extracted_df_info)==x),
#                                                          "Scaled_Fragment_Count"]
#     coverage_zscore = (chrN_sample_total_fragment - chrN_ref_total_fragment_mean) / chrN_ref_total_fragment_sd
#
#     c("Chromosome" = chr_position[1],
#       "Start" = chr_position[2],
#       "End" = chr_position[3],
#       "chrN_ref_total_fragment_mean"= chrN_ref_total_fragment_mean,
#       "chrN_ref_total_fragment_sd" = chrN_ref_total_fragment_sd,
#       "chrN_sample_total_fragment" = chrN_sample_total_fragment,
#       "coverage_zscore" = coverage_zscore
#     )
#   })))
# }

z_score_df_no_NaN = .coverage_zscore(merged_count_df,control_scaled_fragment_count_df)

# z_score_df_no_NaN = .coverage_zscore(reads_extracted_df_info,
                                     # control_scaled_fragment_count_df,binsize)
z_score_df_no_NaN = z_score_df_no_NaN[which(z_score_df_no_NaN$coverage_zscore != 'NaN' ),]

z_score_df_no_NaN$Chromosome = as.character(levels(z_score_df_no_NaN$Chromosome))[z_score_df_no_NaN$Chromosome]
z_score_df_no_NaN$Start = as.numeric(levels(z_score_df_no_NaN$Start))[z_score_df_no_NaN$Start]
z_score_df_no_NaN$End = as.numeric(levels(z_score_df_no_NaN$End))[z_score_df_no_NaN$End]
z_score_df_no_NaN$chrN_ref_total_fragment_mean = as.numeric(levels(z_score_df_no_NaN$chrN_ref_total_fragment_mean))[z_score_df_no_NaN$chrN_ref_total_fragment_mean]
z_score_df_no_NaN$chrN_ref_total_fragment_sd = as.numeric(levels(z_score_df_no_NaN$chrN_ref_total_fragment_sd))[z_score_df_no_NaN$chrN_ref_total_fragment_sd]
z_score_df_no_NaN$chrN_sample_total_fragment = as.numeric(levels(z_score_df_no_NaN$chrN_sample_total_fragment))[z_score_df_no_NaN$chrN_sample_total_fragment]
z_score_df_no_NaN$coverage_zscore = as.numeric(levels(z_score_df_no_NaN$coverage_zscore))[z_score_df_no_NaN$coverage_zscore]

# z_score_df_no_NaN$scaledPos = (z_score_df_no_NaN$Start + z_score_df_no_NaN$End)/2 + as.numeric(chroffsets[z_score_df_no_NaN$Chromosome])
chromosomes <- setNames(1:length(unique(z_score_df_no_NaN$Chromosome)),
                        unique(z_score_df_no_NaN$Chromosome))

z_score_df_no_NaN_file = paste0(output_folder,"/",output_file_name,"_z_score_",binsize,".csv")
write.table(x = z_score_df_no_NaN,
            file = z_score_df_no_NaN_file,
            quote = FALSE,row.names = FALSE,sep = "\t")

#
#
# .smoothing_zscore.smooth_zscore = function(z_score_df){
#   z_score_df$chr_arm  = unlist(apply(z_score_df,1, function(z_score_df_position){
#     filter(merge_cytoband_df,chr==as.character(z_score_df_position["Chromosome"]),
#            start<=as.integer(z_score_df_position["Start"]),
#            end>=as.integer(z_score_df_position["End"])) %>% select(chr_arm)
#   }))
#
#   all_chr_arm = unique(z_score_df$chr_arm)
#
#   smoothed_z_score = sapply(all_chr_arm, function(chr.arm){
#     zscore = filter(z_score_df,chr_arm==chr.arm) %>% select(coverage_zscore)
#     # print(chr.arm)
#     if (length(unlist(zscore))>smooth_moving_window ) {
#       as.numeric(get_ma(unlist(zscore)))
#     } else as.numeric(unlist(zscore))
#
#   })
#   z_score_df$smoothed_zscore = unlist(smoothed_z_score)
#   return(z_score_df)
# }



#
#
# .plot_zscore <- function(z_score_df_no_NaN,binsize){
#   plot(
#     z_score_df_no_NaN$scaledPos,
#     z_score_df_no_NaN$coverage_zscore,
#     type=ifelse(binsize >= 1000, "l", "p"),
#     cex=0.3,
#     pch=19,
#     col=ifelse(chromosomes[z_score_df_no_NaN$Chromosome] %% 2 == 0, "red", "black"),
#     ylim=c(plot_lowerbound, plot_upperbound),
#     xaxt="n",
#     ylab=paste("z-score (binsize=",binsize," kb)",sep = ""),
#     xlab="",
#     cex.axis=2,
#     cex.lab=2
#   )
#   regionsOffTheChart <- z_score_df_no_NaN[z_score_df_no_NaN$coverage_zscore > plot_upperbound | z_score_df_no_NaN$coverage_zscore < plot_lowerbound,]
#   points(
#     regionsOffTheChart$scaledPos,
#     ifelse(regionsOffTheChart$coverage_zscore < plot_lowerbound, plot_lowerbound, plot_upperbound),
#     pch=ifelse(regionsOffTheChart$coverage_zscore < plot_lowerbound, 6, 2),
#     col=ifelse(chromosomes[regionsOffTheChart$Chromosome] %% 2 == 0, "red", "black")
#   )
#   text(chrMids, rep(plot_lowerbound, length(unique(z_score_df_no_NaN$Chromosome))), unique(z_score_df_no_NaN$Chromosome), cex=2)
#   for (x in c(0, cumsum(as.numeric(chrLength))))
#     lines(c(x, x), c(plot_lowerbound, plot_upperbound), col="grey", lty=2)
#   abline(h = 0 , col="blue")
#   abline(h = cutoff_zscore , col="red",lwd=3)
# }
# .plot_smoothed_zscore <- function(z_score_df_no_NaN,binsize){
#   plot(
#     z_score_df_no_NaN$scaledPos,
#     z_score_df_no_NaN$smoothed_zscore,
#     type=ifelse(binsize >= 1000, "l", "p"),
#     cex=0.2,
#     pch=19,
#     col=ifelse(chromosomes[z_score_df_no_NaN$Chromosome] %% 2 == 0, "red", "black"),
#     ylim=c(plot_lowerbound, plot_upperbound),
#     xaxt="n",
#     ylab=paste("z-score (binsize=",binsize," kb)",sep = ""),
#     xlab="",
#     cex.axis=2,
#     cex.lab=2
#   )
#   regionsOffTheChart <- z_score_df_no_NaN[z_score_df_no_NaN$coverage_zscore > plot_upperbound | z_score_df_no_NaN$coverage_zscore < plot_lowerbound,]
#   points(
#     regionsOffTheChart$scaledPos,
#     ifelse(regionsOffTheChart$coverage_zscore < plot_lowerbound, plot_lowerbound, plot_upperbound),
#     pch=ifelse(regionsOffTheChart$coverage_zscore < plot_lowerbound, 6, 2),
#     col=ifelse(chromosomes[regionsOffTheChart$Chromosome] %% 2 == 0, "red", "black")
#   )
#   text(chrMids, rep(plot_lowerbound, length(unique(z_score_df_no_NaN$Chromosome))), unique(z_score_df_no_NaN$Chromosome), cex=2)
#   for (x in c(0, cumsum(as.numeric(chrLength))))
#     lines(c(x, x), c(plot_lowerbound, plot_upperbound), col="grey", lty=2)
#   abline(h = 0 , col="blue")
#   abline(h = cutoff_zscore , col="red",lwd=3)
# }
# paste("Plotting result to ",output_folder,"/",output_file_name,"_plots.png", sep="")
# png(paste(output_folder,"/",output_file_name,"_plots.png", sep=""), width=1800, height=600, units="px")
# .plot_zscore(z_score_df_no_NaN,binsize)
# dev.off()
# paste("Plotting result to ",output_folder,"/",output_file_name,"_smoothedplots.png", sep="")
# png(paste(output_folder,"/",output_file_name,"_smoothedplots.png", sep=""), width=1800, height=600, units="px")
# z_score_df_no_NaN = .smoothing_zscore.smooth_zscore(z_score_df_no_NaN)
# .plot_smoothed_zscore(z_score_df_no_NaN,binsize)
# dev.off()
