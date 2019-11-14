library(PSCBS)
.file_exists <- function(files_vector){
  all_exist = TRUE
  for (file in files_vector) {
    if(!file.exists(file)){
      print(paste(file,"do not existed",sep=" "))
      all_exist = FALSE
    }else{
      print(paste(file,"existed",sep=" "))
    }
  }
  return(all_exist)
}


.segmentation_CBS <- function(z_score_df_no_NaN,delta_F_binsize){
  data <- z_score_df_no_NaN[,c("Chromosome","scaledPos","size_based_zscore")]
  colnames(data) <- c("chromosome","x","y")
  data <- dropSegmentationOutliers(data)
  gaps <- findLargeGaps(data, minLength = 1e+06)
  knownSegments <- gapsToSegments(gaps)
  fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
  fitP <- pruneByHClust(fit, h = 0.75, verbose = -10, distMethod = "euclidean",hclustMethod = "centroid")

  segment_df = getSegments(fitP, simplify = TRUE)
  segment_df$pass_cutoff = segment_df$mean >= cutoff_zscore[1] | segment_df$mean <= cutoff_zscore[2]
  print("Drawing Segmentation Lines")
  segments(x0 = segment_df[!segment_df$pass_cutoff,]$start, y0 = segment_df[!segment_df$pass_cutoff,]$mean,
           x1 = segment_df[!segment_df$pass_cutoff,]$end, y1 = segment_df[!segment_df$pass_cutoff,]$mean,
           col="yellow", lwd=5, lty=2)
  segments(x0 = segment_df[segment_df$pass_cutoff,]$start, y0 = segment_df[segment_df$pass_cutoff,]$mean,
           x1 = segment_df[segment_df$pass_cutoff,]$end, y1 = segment_df[segment_df$pass_cutoff,]$mean,
           col="green", lwd=7)
  print(paste("Performing PSCBS Segmentation for binsize",delta_F_binsize,"(Done)"))
  print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv",sep = ""))
  write.table(x = segment_df,file = paste(output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv", sep=""),quote = FALSE,row.names = FALSE)
  print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv"," (Done)",sep = ""))

}

.smooth_segmentation_CBS <- function(z_score_df_no_NaN,delta_F_binsize){
  data <- z_score_df_no_NaN[,c("Chromosome","scaledPos","smoothed_zscore")]
  colnames(data) <- c("chromosome","x","y")
  data <- dropSegmentationOutliers(data)
  gaps <- findLargeGaps(data, minLength = 500000)
  knownSegments <- gapsToSegments(gaps)
  fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
  fitP <- pruneByHClust(fit, h = 0.75, verbose = -10, distMethod = "euclidean",hclustMethod = "centroid")

  segment_df = getSegments(fitP, simplify = TRUE)
  segment_df$pass_cutoff = segment_df$mean >= cutoff_zscore[1] | segment_df$mean <= cutoff_zscore[2]
  print("Drawing Segmentation Lines")
  segments(x0 = segment_df[!segment_df$pass_cutoff,]$start, y0 = segment_df[!segment_df$pass_cutoff,]$mean,
           x1 = segment_df[!segment_df$pass_cutoff,]$end, y1 = segment_df[!segment_df$pass_cutoff,]$mean,
           col="yellow", lwd=5, lty=2)
  segments(x0 = segment_df[segment_df$pass_cutoff,]$start, y0 = segment_df[segment_df$pass_cutoff,]$mean,
           x1 = segment_df[segment_df$pass_cutoff,]$end, y1 = segment_df[segment_df$pass_cutoff,]$mean,
           col="green", lwd=7)
  print(paste("Performing PSCBS Segmentation for binsize",delta_F_binsize,"(Done)"))
  print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv",sep = ""))
  write.table(x = segment_df,file = paste(output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv", sep=""),quote = FALSE,row.names = FALSE)
  print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentation_",delta_F_binsize,".csv"," (Done)",sep = ""))

}



args <- commandArgs(T)
if (length(args) < 2)
  stop("Usage: PSCBS_segmentation.R *runtime_variable_file(runtime_variables.R) *config_file")
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

output_folder_temp = paste0(output_folder,"/",segmentation_outdir_prefix)

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

output_file_name <- paste(tools::file_path_sans_ext(basename(z_score_df_file)),
                          "_segmentation",sep="")
paste0("Set output prefix to ",output_file_name)

z_score_df_no_NaN = read.delim(z_score_df_file,sep="\t",stringsAsFactors = FALSE)


chromosomes <- setNames(1:length(unique(z_score_df_no_NaN$Chromosome)),
                        unique(z_score_df_no_NaN$Chromosome))


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

paste("Plotting result to ",output_folder,"/",output_file_name,"_plots.png", sep="")
png(paste(output_folder,"/",output_file_name,"_plots.png", sep=""), width=1800, height=1350, units="px")
par(mfrow=c(2,1))
main_functions.plot_zscore(z_score_df_no_NaN,binsize)
.segmentation_CBS(z_score_df_no_NaN,binsize)
main_functions.plot_smoothed_zscore(z_score_df_no_NaN,binsize)
# .smooth_segmentation_CBS(z_score_df_no_NaN,binsize)
# main_functions.plot_zscore(z_score_df_no_NaN_1000,1000)
# main_functions.plot_zscore(z_score_df_no_NaN_5000,5000)
dev.off()
paste("Plotting result to ",output_folder,"/",output_file_name,"_plots.png (Done)", sep="")


