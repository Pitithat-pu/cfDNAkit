## Check required Bioconductor packages
## install if required
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("QDNAseq","GenomicRanges","Rsamtools","PSCBS","dplyr","ggplot2","gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
## Load packages
library(QDNAseq)
library(GenomicRanges)
library(Rsamtools)
library(PSCBS)
# library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
##
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
main_functions.file_exists <- function(files_vector){
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
main_functions.create_blacklist_gr <- function(blacklist_regions){
  if(!main_functions.file_exists(blacklist_regions)) q('no')
  blacklist_targets_gr = GRanges()
  for (blacklist_region in blacklist_regions) {
    blacklist_targets = as.data.frame(read.table(gzfile(blacklist_region),
                                                 header = FALSE,sep = "\t"))[,1:3]
    colnames(blacklist_targets) = c("chromosome","start","end")
    blacklist_targets_gr_temp=GRanges(seqnames = blacklist_targets$chromosome,
                                 ranges = IRanges(start = as.numeric(blacklist_targets$start),
                                                  end=as.numeric(blacklist_targets$end)))

    blacklist_targets_gr = c(blacklist_targets_gr,blacklist_targets_gr_temp)
  }
  return(blacklist_targets_gr)
}
.filter_blacklist_regions <- function(which,blacklist_targets_gr){
  filtered_gr <- GenomicRanges::setdiff(which,blacklist_targets_gr)
  return(filtered_gr)
}
main_function.get_sliding_windows <- function(binsize){
  if (main_functions.file_exists(c(qdnaseq_sliding_windows_RDS))) {
    paste0("Reading QDNAseq sliding windows ",qdnaseq_sliding_windows_RDS)
    bins = readRDS(qdnaseq_sliding_windows_RDS)
  }else{
    paste0("Load new QDNAseq sliding windows")
    bins = getBinAnnotations(binSize=binsize)
    paste0("Save new QDNAseq sliding windows",qdnaseq_sliding_windows_RDS)
    saveRDS(bin,qdnaseq_sliding_windows_RDS)
  }
  sliding_windows = as.data.frame(bins@data)
  sliding_windows = sliding_windows[which(sliding_windows$chromosome!="Y" &
                                            sliding_windows$mappability>=1),]
  return(sliding_windows)
}

main_functions.readbam <- function(sliding_windows){
  chrom <- sliding_windows[1]
  start <- sliding_windows[2]
  end <- sliding_windows[3]
  # print(paste0(chrom," ",start," ",end))
  sliding_windows_gr <- GRanges(seqnames = chrom, ranges = IRanges(start=as.numeric(start), end=as.numeric(end)))

  if (exists('capture_targets_gr')) which <- subsetByOverlaps(capture_targets_gr,sliding_windows_gr)
  else which = sliding_windows_gr
  which_filtered <- .filter_blacklist_regions(which,blacklist_targets_gr)
  if (length(which_filtered)==0) {
    empty_list = list("qname" = as.character(), "rname" = factor(chrom),
                      "pos" = as.integer(),"isize"=as.integer())
    return(empty_list)
  }
  ## extract all forward-strand read

  flag=scanBamFlag(isPaired = TRUE,
                   isUnmappedQuery = FALSE,
                   isDuplicate = FALSE,isMinusStrand = FALSE,
                   hasUnmappedMate = FALSE, isSecondaryAlignment = FALSE,isMateMinusStrand = TRUE,
                   isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE)


  # flag=scanBamFlag(isPaired = TRUE, isProperPair=TRUE,
  #                  isUnmappedQuery = FALSE, isFirstMateRead = TRUE,
  #                  isDuplicate = FALSE,isSecondMateRead = FALSE,
  #                  hasUnmappedMate = FALSE, isNotPassingQualityControls = FALSE,
  #                  isSupplementaryAlignment = FALSE)

  param <- ScanBamParam(what = what, flag = flag, which = which_filtered)
  bam <- scanBam(file = sample_bamfile,
                 index = sample_bamindex,
                 param=param)
  bam <- unname(bam)
  elts <- setNames(bamWhat(param), bamWhat(param))
  lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))

}

main_functions.bias_correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias)
  coverage.model <- loess(predict(coverage.trend, i) ~ i)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <- coverage - coverage.pred + median(coverage)
}

get_middle_BAFlevel <- function(BAFLevel_file){
  print(paste("Reading Control BAFLevel matrix ",BAFLevel_file,sep=""))
  BAFLevel_df = read.table(file = BAFLevel_file,sep="\t", header = TRUE)
  print(paste("Reading Control BAFLevel matrix ",BAFLevel_file," (Done)",sep=""))
  BAFLevel_df[which(!is.na(BAFLevel_df$middleLevel)),]
}





main_functions.calculate_deltaF <- function(SLRatio, SLRatio_reference){
  SLRatio - SLRatio_reference
}
main_function.calculate_log2SLratio <- function(SLRatio, SLRatio_reference){
  log2(SLRatio/SLRatio_reference)
}

.get_delta_f_corrected <- function(test_region_name,reference_region){
  test_region = reads_extracted_df_info[test_region_name,]
  SLRatio_corrected = test_region$SLRatio_corrected
  SLRatio_reference = sum(reference_region$short_corrected, na.rm = TRUE) / sum(reference_region$long_corrected, na.rm = TRUE)
  c("ReadCount"=test_region$`Read Pairs in range_corrected`,
    "Average_Size"=mean(test_region$Mean),
    "SLRatio"= SLRatio_corrected,
    "SLRatio_reference"=SLRatio_reference,
    "Delta_F"= main_functions.calculate_deltaF(SLRatio_corrected, SLRatio_reference),# "Delta_F"= main_function.calculate_log2SLratio(SLRatio_corrected, SLRatio_reference)
    "log2SLratio"= main_function.calculate_log2SLratio(SLRatio_corrected, SLRatio_reference)
  )
}

main_functions.get_short_ratio_corrected <- function(reads_extracted_df_info, all_region_name){
  test_regions_df = as.data.frame(t(sapply(all_region_name,function(test_region_name){
    reference_region = reads_extracted_df_info
    reference_region = reference_region[!rownames(reference_region) %in% test_region_name,]
    ## Set X chromosome as reference regions
    # reference_region = reference_region[which(startsWith(rownames(reference_region),"X")),]
    .get_delta_f_corrected(test_region_name, reference_region)
  })))
  return(test_regions_df)
}
get_short_ratio_BAF <- function(reads_extracted_df_info, all_region_name, BAF_reference_region){
  test_regions_df = as.data.frame(t(sapply(all_region_name,function(test_region_name){
    get_delta_f(test_region_name, reference_region = BAF_reference_region)
  })))

  return(test_regions_df)
}

.get_zscore_df <- function(all_region_name,test_regions_df,delta_f_control_df){
  zscore_df = as.data.frame(t(sapply(all_region_name, function(x){
    region_deltaf = unlist(delta_f_control_df[x,])
    chr_position = strsplit(x,"[:-]+")[[1]]
    chrN_ref_delta_F_mean = mean(region_deltaf,na.rm = TRUE)
    chrN_ref_delta_F_sd = sd(region_deltaf,na.rm = TRUE)
    chrN_sample_delta_F = test_regions_df[which(rownames(test_regions_df)==x),"Delta_F"]
    chrN_sample_log2SLratio = test_regions_df[which(rownames(test_regions_df)==x),"log2SLratio"]
    size_based_zscore = (chrN_sample_delta_F - chrN_ref_delta_F_mean) / chrN_ref_delta_F_sd
    # size_based_zscore = (chrN_sample_log2SLratio - chrN_ref_delta_F_mean) / chrN_ref_delta_F_sd

    c("Chromosome" = chr_position[1],
      "Start" = chr_position[2],
      "End" = chr_position[3],
      "chrN_ref_delta_F_mean"= chrN_ref_delta_F_mean,
      "chrN_ref_delta_F_sd" = chrN_ref_delta_F_sd,
      "chrN_sample_delta_F" = chrN_sample_delta_F,
      "chrN_sample_log2SLratio" = chrN_sample_log2SLratio,
      "size_based_zscore" = size_based_zscore
    )
  })))
  zscore_df$Chromosome = as.character(levels(zscore_df$Chromosome))[zscore_df$Chromosome]
  zscore_df$Start = as.numeric(levels(zscore_df$Start))[zscore_df$Start]
  zscore_df$End = as.numeric(levels(zscore_df$End))[zscore_df$End]
  zscore_df$chrN_ref_delta_F_mean = as.numeric(levels(zscore_df$chrN_ref_delta_F_mean))[zscore_df$chrN_ref_delta_F_mean]
  zscore_df$chrN_ref_delta_F_sd = as.numeric(levels(zscore_df$chrN_ref_delta_F_sd))[zscore_df$chrN_ref_delta_F_sd]
  zscore_df$chrN_sample_delta_F = as.numeric(levels(zscore_df$chrN_sample_delta_F))[zscore_df$chrN_sample_delta_F]
  zscore_df$size_based_zscore = as.numeric(levels(zscore_df$size_based_zscore))[zscore_df$size_based_zscore]
  zscore_df$chrN_sample_log2SLratio = as.numeric(levels(zscore_df$chrN_sample_log2SLratio))[zscore_df$chrN_sample_log2SLratio]
  return(zscore_df)
}

.get_zscore_df <- function(all_region_name,test_regions_df,delta_f_control_df){
  zscore_df = as.data.frame(t(sapply(all_region_name, function(x){
    region_deltaf = unlist(delta_f_control_df[x,])
    chr_position = strsplit(x,"[:-]+")[[1]]
    chrN_ref_delta_F_mean = mean(region_deltaf,na.rm = TRUE)
    chrN_ref_delta_F_sd = sd(region_deltaf,na.rm = TRUE)
    chrN_sample_delta_F = test_regions_df[which(rownames(test_regions_df)==x),"Delta_F"]
    chrN_sample_log2SLratio = test_regions_df[which(rownames(test_regions_df)==x),"log2SLratio"]
    size_based_zscore = (chrN_sample_delta_F - chrN_ref_delta_F_mean) / chrN_ref_delta_F_sd
    # size_based_zscore = (chrN_sample_log2SLratio - chrN_ref_delta_F_mean) / chrN_ref_delta_F_sd

    c("Chromosome" = chr_position[1],
      "Start" = chr_position[2],
      "End" = chr_position[3],
      "chrN_ref_delta_F_mean"= chrN_ref_delta_F_mean,
      "chrN_ref_delta_F_sd" = chrN_ref_delta_F_sd,
      "chrN_sample_delta_F" = chrN_sample_delta_F,
      "chrN_sample_log2SLratio" = chrN_sample_log2SLratio,
      "size_based_zscore" = size_based_zscore
    )
  })))
  zscore_df$Chromosome = as.character(levels(zscore_df$Chromosome))[zscore_df$Chromosome]
  zscore_df$Start = as.numeric(levels(zscore_df$Start))[zscore_df$Start]
  zscore_df$End = as.numeric(levels(zscore_df$End))[zscore_df$End]
  zscore_df$chrN_ref_delta_F_mean = as.numeric(levels(zscore_df$chrN_ref_delta_F_mean))[zscore_df$chrN_ref_delta_F_mean]
  zscore_df$chrN_ref_delta_F_sd = as.numeric(levels(zscore_df$chrN_ref_delta_F_sd))[zscore_df$chrN_ref_delta_F_sd]
  zscore_df$chrN_sample_delta_F = as.numeric(levels(zscore_df$chrN_sample_delta_F))[zscore_df$chrN_sample_delta_F]
  zscore_df$size_based_zscore = as.numeric(levels(zscore_df$size_based_zscore))[zscore_df$size_based_zscore]
  # zscore_df$chrN_sample_log2SLratio = as.numeric(levels(zscore_df$chrN_sample_log2SLratio))[zscore_df$chrN_sample_log2SLratio]
  return(zscore_df)
}


main_functions.get_zscore_per_binsize <- function(test_regions_df_corrected,control_deltaF_df,delta_F_binsize){
  # if (delta_F_binsize != binsize){
  #   test_regions_df_corrected = .merge_regions(test_regions_df_corrected,delta_F_binsize)
  #   test_regions_df_corrected = test_regions_df_corrected[,-1]
  #   control_deltaF_df = .merge_regions(control_deltaF_df,delta_F_binsize)
  #   control_deltaF_df = control_deltaF_df[,-1]
  # }
  all_region_name = rownames(test_regions_df_corrected)

  # print(all_region_name)
  z_score_df_corrected = .get_zscore_df(all_region_name,
                                       test_regions_df = test_regions_df_corrected,
                                       delta_f_control_df = control_deltaF_df)
  # if (!delta_F_binsize == binsize){
  #   z_score_df_corrected = .merge_regions(z_score_df_corrected[,c(-1,-2,-3)],delta_F_binsize)
  #   z_score_df_corrected = z_score_df_corrected[,-1]
  #   chr_position = strsplit(rownames(z_score_df_corrected),"[:-]+")
  #   z_score_df_corrected$Chromosome = as.character(lapply(chr_position, `[[`, 1))
  #   z_score_df_corrected$Start = as.numeric(lapply(chr_position, `[[`, 2))
  #   z_score_df_corrected$End = as.numeric(lapply(chr_position, `[[`, 3))
  # }
  # z_score_df_no_NaN = z_score_df_corrected[which(z_score_df_corrected$size_based_zscore != 'NaN' ),]
  # z_score_df_no_NaN$scaledPos = (z_score_df_no_NaN$Start + z_score_df_no_NaN$End)/2 +
  #   as.numeric(chroffsets[z_score_df_no_NaN$Chromosome])

  return(z_score_df_corrected)
}

.n.colmedians <- function(df, n = 10){
  aggregate(x = df,
            by = list(gl(ceiling(nrow(df)/n), n)[1:nrow(df)]),
            function(x) median(x, na.rm=TRUE))
}

.seqlast <- function (from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}
.sum_square <- function(x){
  sum(x^2,na.rm = TRUE)
}
.merge_regions <- function (test_regions_corrected,delta_F_binsize){
  test_regions_expand_df = data.frame()
  all_chr=unique(sapply(strsplit(rownames(test_regions_corrected),":"), `[`, 1))
  for (chr in all_chr) {
    test_regions_by_chr = test_regions_corrected[which(startsWith(rownames(test_regions_corrected), paste(chr,":",sep=""))),]
    test_regions_temp = .n.colmedians(test_regions_by_chr,delta_F_binsize/binsize)
    start_region_index = seq(from = 1, to = nrow(test_regions_by_chr),
                             by = min(nrow(test_regions_by_chr),delta_F_binsize/binsize))
    # end_region_index = ifelse(nrow(test_regions_by_chr)<(delta_F_binsize/binsize), nrow(test_regions_by_chr),
                              # .seqlast(from = delta_F_binsize/binsize, to = nrow(test_regions_by_chr), by = delta_F_binsize/binsize))
    if (nrow(test_regions_by_chr)<(delta_F_binsize/binsize))
      end_region_index = nrow(test_regions_by_chr)
    else
      end_region_index = .seqlast(from = delta_F_binsize/binsize, to = nrow(test_regions_by_chr),
                               by = delta_F_binsize/binsize)

    # print(end_region_index)
    test_region_start = sapply(strsplit(rownames(test_regions_by_chr)[start_region_index],"-"), `[`, 1)
    test_region_end = sapply(strsplit(rownames(test_regions_by_chr)[end_region_index],"-"), `[`, 2)
    new_test_regions = paste(test_region_start,test_region_end,sep="-")
    rownames(test_regions_temp) = new_test_regions
    test_regions_expand_df = rbind(test_regions_expand_df,test_regions_temp)
  }
  return(test_regions_expand_df)
}

main_functions.calculate_gw_zscore <- function(sample_z_score,control_zscore_df){
  control_sum_square = sapply(control_zscore_df, .sum_square)
  sample_sum_square = .sum_square(sample_z_score)
  gw_zscore = (sample_sum_square - mean(control_sum_square)) / sd(control_sum_square)
}


main_functions.plot_smoothed_zscore <- function(z_score_df_no_NaN,delta_F_binsize){
  plot(
    z_score_df_no_NaN$scaledPos,
    z_score_df_no_NaN$smoothed_zscore,
    type=ifelse(delta_F_binsize >= 1000, "l", "p"),
    cex=delta_F_binsize/500,
    pch=19,
    col=ifelse(chromosomes[z_score_df_no_NaN$Chromosome] %% 2 == 0, "red", "black"),
    ylim=c(plot_lowerbound, plot_upperbound),
    xaxt="n",
    ylab=paste("z-score (binsize=",delta_F_binsize," kb)",sep = ""),
    xlab="",
    cex.axis=2,
    cex.lab=2
  )
  regionsOffTheChart <- z_score_df_no_NaN[z_score_df_no_NaN$smoothed_zscore > plot_upperbound | z_score_df_no_NaN$smoothed_zscore < plot_lowerbound,]
  points(
    regionsOffTheChart$scaledPos,
    ifelse(regionsOffTheChart$smoothed_zscore < plot_lowerbound, plot_lowerbound, plot_upperbound),
    pch=ifelse(regionsOffTheChart$smoothed_zscore < plot_lowerbound, 6, 2),
    col=ifelse(chromosomes[regionsOffTheChart$Chromosome] %% 2 == 0, "red", "black")
  )
  text(chrMids, rep(plot_lowerbound, length(unique(z_score_df_no_NaN$Chromosome))), unique(z_score_df_no_NaN$Chromosome), cex=1.5)
  for (x in c(0, cumsum(as.numeric(chrLength))))
    lines(c(x, x), c(plot_lowerbound, plot_upperbound), col="grey", lty=2)
  abline(h = 0 , col="blue")
  abline(h = cutoff_zscore , col="red")

  # if (delta_F_binsize <= 100) {
  #   print(paste("Performing PSCBS Segmentation for binsize",delta_F_binsize))
  #   .smooth_segmentation_CBS(z_score_df_no_NaN,delta_F_binsize)
  # }
  #
}

main_functions.plot_zscore <- function(z_score_df_no_NaN,delta_F_binsize){
  plot(
    z_score_df_no_NaN$scaledPos,
    z_score_df_no_NaN$size_based_zscore,
    type=ifelse(delta_F_binsize >= 1000, "l", "p"),
    cex=delta_F_binsize/500,
    pch=19,
    col=ifelse(chromosomes[z_score_df_no_NaN$Chromosome] %% 2 == 0, "red", "black"),
    ylim=c(plot_lowerbound, plot_upperbound),
    xaxt="n",
    ylab=paste("z-score (binsize=",delta_F_binsize," kb)",sep = ""),
    xlab="",
    cex.axis=2,
    cex.lab=2
  )
  regionsOffTheChart <- z_score_df_no_NaN[z_score_df_no_NaN$size_based_zscore > plot_upperbound | z_score_df_no_NaN$size_based_zscore < plot_lowerbound,]
  points(
    regionsOffTheChart$scaledPos,
    ifelse(regionsOffTheChart$size_based_zscore < plot_lowerbound, plot_lowerbound, plot_upperbound),
    pch=ifelse(regionsOffTheChart$size_based_zscore < plot_lowerbound, 6, 2),
    col=ifelse(chromosomes[regionsOffTheChart$Chromosome] %% 2 == 0, "red", "black")
  )
  text(chrMids, rep(plot_lowerbound, length(unique(z_score_df_no_NaN$Chromosome))), unique(z_score_df_no_NaN$Chromosome), cex=2)
  for (x in c(0, cumsum(as.numeric(chrLength))))
    lines(c(x, x), c(plot_lowerbound, plot_upperbound), col="grey", lty=2)
  abline(h = 0 , col="blue")
  abline(h = cutoff_zscore , col="red",lwd=3)
}
