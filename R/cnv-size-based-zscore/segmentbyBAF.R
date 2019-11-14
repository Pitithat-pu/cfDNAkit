source("/home/puranach/Documents/Rscripts/cfDNAkit/kit_functions.R")
source("/home/puranach/Documents/Rscripts/cfDNAkit/set_parameters.R")
args <- commandArgs(T)
if (length(args) != 4)
  stop("Usage: segmentbyBAF.R reads_extracted_info_file.tsv deltaf_control_file BAFLevel_file output_folder\n
       deltaF_controlfile locate at /icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/BH01_rmdup_paired_mapped_deltaf_table_n10.csv \n
       BAFLevel_file the result of BAF_segmentation.R \n")
reads_extracted_info_file <- args[1]
delta_f_control_file <- args[2]
BAFLevel_file <- args[3]

output_folder <- args[4]
output_file_name <- paste(tools::file_path_sans_ext(basename(BAFLevel_file)),sep="")

print(paste("Reading Control DeltaF matrix ",delta_f_control_file,sep=""))
delta_f_control_df = read.table(file = delta_f_control_file,sep = "\t")
print(paste("Reading Control DeltaF matrix ",delta_f_control_file," (Done)",sep=""))

BAFLevel_df = data.frame(get_middle_BAFlevel(BAFLevel_file))
BAF_region = as.character(BAFLevel_df$segmentID)
BAF_chr_position = strsplit(BAF_region,"[:-]+")
BAF_chromosome = sapply(BAF_chr_position, function(x){x[1]})
BAF_start = sapply(BAF_chr_position, function(x){x[2]})
BAF_end = sapply(BAF_chr_position, function(x){x[3]})
BAF_region_gr = GRanges(seqnames = BAF_chromosome,
                        ranges = IRanges(start = as.numeric(BAF_start), end=as.numeric(BAF_end)))   
print(paste("Reading read_info_file matrix ",reads_extracted_info_file,sep=""))
reads_extracted_df_info = read.table(reads_extracted_info_file,header = TRUE)
print(paste("Reading read_info_file matrix ",reads_extracted_info_file," (Done)",sep=""))

reads_extracted_df_position = strsplit(rownames(reads_extracted_df_info),"[:-]+")
reads_extracted_df_chromosome = sapply(reads_extracted_df_position, function(x){x[1]})
reads_extracted_df_start = sapply(reads_extracted_df_position, function(x){x[2]})
reads_extracted_df_end = sapply(reads_extracted_df_position, function(x){x[3]})
reads_extracted_df_gr = GRanges(seqnames = reads_extracted_df_chromosome,
                                ranges = IRanges(start = as.numeric(reads_extracted_df_start), end=as.numeric(reads_extracted_df_end)))
print("Getting read info with BAF segment (0.5 AF) by findOverlaps")
reference_overlap = findOverlaps(query = reads_extracted_df_gr, subject = BAF_region_gr, type = "within")
reference_region = reads_extracted_df_info[reference_overlap@from,]


all_region_name = rownames(reads_extracted_df_info)
test_regions_df = get_short_ratio_BAF(reads_extracted_df_info = reads_extracted_df_info,
                                  all_region_name = all_region_name,
                                  BAF_reference_region = reference_region)

print("Getting Test Regions info (Done)")

z_score_df = get_zscore_MAD_df(all_region_name,
                           test_regions_df = test_regions_df,
                           delta_f_control_df = delta_f_control_df)



print(paste("Writing size-based z-score analysis to ",output_folder,"/",output_file_name,"_final_zscore.csv",sep = ""))
write.table(x = z_score_df,file = paste(output_folder,"/",output_file_name,"_final_zscore.csv", sep=""),quote = FALSE,row.names = FALSE)
print(paste("Writing size-based z-score analysis to ",output_folder,"/",output_file_name,"_final_zscore.csv"," (Done)",sep = ""))

chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
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

chromosomes <- setNames(1:length(unique(z_score_df$Chromosome)), unique(z_score_df$Chromosome))
z_score_df_no_NaN = z_score_df[which(z_score_df$size_based_zscore != 'NaN' ),]
z_score_df_no_NaN$scaledPos = (z_score_df_no_NaN$Start + z_score_df_no_NaN$End)/2 + 
  as.numeric(chroffsets[z_score_df_no_NaN$Chromosome])

print(paste("Plotting result to ",output_folder,"/",output_file_name,"_segmentedbyBAF.png"))
png(paste(output_folder,"/",output_file_name,"_segmentedbyBAF.png", sep=""), width=2000, height=800, units="px")
plot(
  z_score_df_no_NaN$scaledPos,
  z_score_df_no_NaN$size_based_zscore,
  cex=0.2,
  pch=19,
  col=ifelse(chromosomes[z_score_df_no_NaN$Chromosome] %% 2 == 0, "red", "black"),
  ylim=c(plot_lowerbound, plot_upperbound),
  xaxt="n",
  ylab="Sized-based z-score",
  xlab="",
  cex.axis=1.5,
  cex.lab=1.5
)
title(main = "Sized-based Copy-number Analysis", sub = tools::file_path_sans_ext(basename(BAFLevel_file)))
regionsOffTheChart <- z_score_df_no_NaN[z_score_df_no_NaN$size_based_zscore > plot_upperbound | z_score_df_no_NaN$size_based_zscore < plot_lowerbound,]
points(
  regionsOffTheChart$scaledPos,
  #sign(regionsOffTheChart$smoothLog2) * ylim,
  ifelse(regionsOffTheChart$size_based_zscore < plot_lowerbound, plot_lowerbound, plot_upperbound),
  pch=ifelse(regionsOffTheChart$size_based_zscore < plot_lowerbound, 6, 2),
  col=ifelse(chromosomes[regionsOffTheChart$Chromosome] %% 2 == 0, "red", "black")
)
text(chrMids, rep(plot_lowerbound, length(unique(z_score_df_no_NaN$Chromosome))), unique(z_score_df_no_NaN$Chromosome), cex=1.5)
for (x in c(0, cumsum(as.numeric(chrLength))))
  lines(c(x, x), c(plot_lowerbound, plot_upperbound), col="grey", lty=2)

abline(h = 0 , col="blue")
abline(h = cutoff_zscore , col="red")

print("Performing PSCBS Segmentation")
library("PSCBS")
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
         col="yellow", lwd=7, lty=2)
segments(x0 = segment_df[segment_df$pass_cutoff,]$start, y0 = segment_df[segment_df$pass_cutoff,]$mean, 
         x1 = segment_df[segment_df$pass_cutoff,]$end, y1 = segment_df[segment_df$pass_cutoff,]$mean, 
         col="green", lwd=7)
print("Performing PSCBS Segmentation (Done)")
# segments(x0 = segment_df$start, y0 = segment_df$mean, 
#          x1 = segment_df$end, y1 = segment_df$mean, 
#          col="green", lwd=5)
dev.off()
print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentedbyBAF.csv",sep = ""))
write.table(x = segment_df,file = paste(output_folder,"/",output_file_name,"_segmentedbyBAF.csv", sep=""),quote = FALSE,row.names = FALSE)
print(paste("Writing segmentation to ",output_folder,"/",output_file_name,"_segmentedbyBAF.csv"," (Done)",sep = ""))

print("Sized-based Fragment Analysis Completed")
