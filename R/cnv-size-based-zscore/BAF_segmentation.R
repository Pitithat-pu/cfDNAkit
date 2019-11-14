#BAF_File = "~/Documents/H021-Z60C.plasma_BAF.tsv.gz"
# BAF_File = "/abi/data2/puranach/S006-C41MJC.plasma03_BAF.tsv.gz"
#BAF_df = read.table(gzfile(BAF_File), header=T)
source("/home/puranach/Documents/Rscripts/cfDNAkit/set_parameters.R")
args <- commandArgs(T)
if (length(args) != 3)
  stop("Usage: BAF_segmentation.R BAF.tsv.gz segmentation.csv output_folder\n")
BAF_file <- args[1]
segments_file <- args[2]
output_folder <- args[3]


output_file_name <- paste(tools::file_path_sans_ext(basename(segments_file)),
                          "_baf_segmentation",sep="")
if(!dir.exists(output_folder)) {
  print(paste("Creating output folder ",output_folder,sep=""))
  dir.create(output_folder,recursive = TRUE)
}  

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
# source("https://bioconductor.org/biocLite.R")
# library(QDNAseq)
# bins <- getBinAnnotations(binSize=binsize)
# segments <- as.data.frame(bins@data)
# segments <- segments[which(segments$chromosome!="Y"),]
# segments$segmentID <- rownames(segments)
segments <- read.table(segments_file, header=T, sep=" ")
segments = segments[which(!is.na(segments$chromosome)),]
segments$chromosome = trimws(segments$chromosome)
segments$start = segments$start - as.numeric(chroffsets[as.character(segments$chromosome)])
segments$end = segments$end - as.numeric(chroffsets[as.character(segments$chromosome)])
segments$segmentID <- paste(segments$chromosome,":",as.integer(segments$start),"-",as.integer(segments$end),sep="")


baf = read.table(gzfile(BAF_file), header=T)
chromosomes <- setNames(1:length(unique(baf$chromosome)), unique(baf$chromosome))
#baf$scaledPos = baf$position + as.numeric(chroffsets[as.character(baf$chromosome)])
baf$baf.control = baf$control_alt_freq/baf$controlCov
#baf.GR = makeGRangesFromDataFrame(baf, start.field = "position", end.field = "position")

baf.forLevelEstimation = baf[baf$baf.control>0.3 & baf$baf.control<0.7,]
# segments$segmentID = rownames(segments)
# segments$width = round((segments$end-segments$start)/1e6,2)
baf$bafID = as.integer(rownames(baf))
baf.forLevelEstimation$bafID = as.integer(rownames(baf.forLevelEstimation))
# d <- density(baf$baf) 
# plot(d)
print("Determine BAF Levels")
BAFLevels = apply(segments[,c("chromosome","start","end","segmentID")], 1, function(segment) {
  baf.part = baf.forLevelEstimation[
    as.character(baf.forLevelEstimation$chromosome)==as.character(segment["chromosome"]) & 
      baf.forLevelEstimation$position>=as.integer(segment["start"]) & 
      baf.forLevelEstimation$position<=as.integer(segment["end"]),]
  # only use BAFdata/SNPs which have a weight (because they are within an exon target region)

  
  # Define tcnNbrOfHets to be the number of heterozygous SNPs
  # That is, which have meet 0.3<=BAF<=0.7
  # heterozygousIndices = which(baf.part$baf>=0.3 & baf.part$baf<=0.7)
  # redundant criterion check...
  heterozygousIndices = which(baf.part$baf.control>=0.3 & baf.part$baf.control<=0.7)
  tcnNbrOfHets = length(heterozygousIndices)
  
  # only accept segments which are represented by at least 5 heterozygous SNPs
  if (nrow(baf.part) < min_heterozygous_SNPs) {
    return (list(segmentID=as.character(segment["segmentID"]), leftLevel=NA, rightLevel=NA, middleLevel=NA, minBAFId=NA, maxBAFId=NA, nrOfPeaks=NA, tcnNbrOfHets=tcnNbrOfHets, meanCovB=NA, meanCovB2=NA))
    # return (-99)
  }

  dens = density(baf.part$baf, bw = 0.05)
  LevelAroundPoint5 = mean(dens$y[dens$x>=0.45 & dens$x<=0.55])
  dens.left = dens$y[dens$x>0.05 & dens$x<0.5]
  dens.right = dens$y[dens$x>0.5  & dens$x<0.95]
  leftPeakBAF = max(dens.left)
  rightPeakBAF = max(dens.right)
  leftLevel = dens$x[which.max(dens.left)+length(which(dens$x<0.05))]
  rightLevel = dens$x[which.max(dens.right) + length(which(dens$x<0.5))]
  if (length(leftLevel)==0 || length(rightLevel)==0) {
    return (list(segmentID=as.character(segment["segmentID"]), leftLevel=NA, rightLevel=NA, middleLevel=NA, minBAFId=NA, maxBAFId=NA, nrOfPeaks=NA, tcnNbrOfHets=NA, meanCovB=NA, meanCovB2=NA))
    # return (-99)
  }
  
  majorAlleleCount = sapply(heterozygousIndices, function(index) {
    count.alt = baf.part[index,"tumor_alt_freq"]
    count.ref = baf.part[index,"tumorCov"] - count.alt
    return( ifelse(count.alt>count.ref,count.alt,count.ref) )
  })
  if (rightLevel-leftLevel < 0.05) {
    # upper and lower BAF cloud are very close to each other --> continue with mean BAF level and check
    # whether it is close to 0.5
    middleLevel = mean(c(leftLevel, rightLevel))
    meanCovB = mean( baf.part[heterozygousIndices,"tumor_alt_freq"] )
    # use meanCovB2 in order to have "correct" meanCovB availabel if segment is set to 2 peaks later in downstream analyses
    meanCovB2 = mean( majorAlleleCount/baf.part[heterozygousIndices,"tumorCov"] ) * segments[segment["segmentID"],"depth"]
    # old formula:          
    # meanCovB = weighted.mean(baf.part[heterozygousIndices,"tumor_alt_freq"], baf.part[heterozygousIndices,"weight"])        
    if (middleLevel-0.5 < 0.05) {
      # if so, set BAF level to 0.5 manually
      return (list(segmentID=as.character(segment["segmentID"]), leftLevel=NA, rightLevel=NA, middleLevel=0.5, minBAFId=min(as.integer(baf.part$bafID)), maxBAFId=max(as.integer(baf.part$bafID)), nrOfPeaks=1, tcnNbrOfHets=tcnNbrOfHets, meanCovB=meanCovB, meanCovB2=meanCovB2))
    } else {
      # otherwise, take the determined mean BAF level
      return (list(segmentID=as.character(segment["segmentID"]), leftLevel=NA, rightLevel=NA, middleLevel=middleLevel, minBAFId=min(as.integer(baf.part$bafID)), maxBAFId=max(as.integer(baf.part$bafID)), nrOfPeaks=1, tcnNbrOfHets=tcnNbrOfHets, meanCovB=meanCovB, meanCovB2=meanCovB2))
    }
  } else {
    # upper and lower BAF cloud are NOT very close to each other
    # is there a local minimum around BAF of 0.5?
    
    meanCovB = mean( majorAlleleCount/baf.part[heterozygousIndices,"tumorCov"] ) * segments[segment["segmentID"],"depth"]
    meanCovB2 = meanCovB
    # old formula:
    # meanCovB = weighted.mean(majorAlleleCount, baf.part[heterozygousIndices,"weight"])
    
    if (LevelAroundPoint5<leftPeakBAF*0.9 & LevelAroundPoint5<rightPeakBAF*0.9) {
      # plot(density(majorAlleleCount), xlim=c(0,100))
      # abline(v=mean(majorAlleleCount))
      # abline(v=segment$depth, col="green")
      return (list(segmentID=as.character(segment["segmentID"]), leftLevel=leftLevel, rightLevel=rightLevel, middleLevel=NA, minBAFId=min(as.integer(baf.part$bafID)), maxBAFId=max(as.integer(baf.part$bafID)), nrOfPeaks=2, tcnNbrOfHets=tcnNbrOfHets, meanCovB=meanCovB, meanCovB2=meanCovB2))
    } else {
      # BAF level of 0.5 has no local minimum --> no two BAF peaks surrounding 0.5 --> only 1 peak with left/right shift
      # --> heterozygous SNPs prefer one allel over the other allel (why?) 
      # something wrong here...?
      # indicate weird nature of BAF peak by return -1 as number of BafPeaks
      if (leftPeakBAF > rightPeakBAF) {
        # BAF peak is on the left side
        return (list(segmentID=as.character(segment["segmentID"]), leftLevel=leftLevel, rightLevel=NA, middleLevel=NA, minBAFId=min(as.integer(baf.part$bafID)), maxBAFId=max(as.integer(baf.part$bafID)), nrOfPeaks=-1, tcnNbrOfHets=tcnNbrOfHets, meanCovB=meanCovB, meanCovB2=meanCovB2))
      } else {
        return (list(segmentID=as.character(segment["segmentID"]), leftLevel=NA, rightLevel=rightLevel, middleLevel=NA, minBAFId=min(as.integer(baf.part$bafID)), maxBAFId=max(as.integer(baf.part$bafID)), nrOfPeaks=-1, tcnNbrOfHets=tcnNbrOfHets, meanCovB=meanCovB, meanCovB2=meanCovB2))
      }
    }
  }
})
BAFLevels = as.data.frame(do.call("rbind", BAFLevels))
segments$nrOfPeaks = unlist(BAFLevels$nrOfPeaks)
segments$tcnNbrOfHets = unlist(BAFLevels$tcnNbrOfHets)
baf$scaledPos <- baf$position + chroffsets[as.character(baf$chromosome)]
print(paste("Plotting BAF segments to ",output_folder,"/",output_file_name,"plot.png"))
png(paste(output_folder,"/",output_file_name,".png", sep=""), width=2000, height=800, units="px")
plot(
  baf$scaledPos,
  baf$baf,
  cex=0.2,
  pch=19,
  col=ifelse(chromosomes[as.character(baf$chromosome)] %% 2 == 0, "red", "black"),
  ylim=c(-0.2, 1.2),
  xaxt="n",
  ylab="BAF",
  xlab="",
  cex.axis=1.5,
  cex.lab=1.5
)
for (x in c(0, cumsum(as.numeric(chrLength))))
  lines(c(x, x), c(-100, 100), col="grey", lty=2)
text(chrMids, rep(-0.2, length(unique(baf$chromosome))), unique(baf$chromosome), cex=1.5)
apply(BAFLevels, 1, function(level) {
  # level = as.vector(BAFLevels[1,])
  scaledPos.start = baf[as.integer(baf$bafID)==as.integer(level["minBAFId"]),"scaledPos"]
  scaledPos.end = baf[as.integer(baf$bafID)==as.integer(level["maxBAFId"]),"scaledPos"]      
  if (!is.na(level["leftLevel"])) { 
    levelValue = as.double(level["leftLevel"])
    lines(c(scaledPos.start,scaledPos.end), c(levelValue,levelValue), col="yellow", lwd=6) 
  }
  if (!is.na(level["rightLevel"])){
    levelValue = as.double(level["rightLevel"])
    lines(c(scaledPos.start,scaledPos.end), c(levelValue,levelValue), col="yellow", lwd=6) 
  } 
  if (!is.na(level["middleLevel"])) {
    levelValue = as.double(level["middleLevel"])
    lines(c(scaledPos.start,scaledPos.end), c(levelValue,levelValue), col="purple", lwd=6) 
  } 
})
dev.off()

BAFLevels_outfile = data.frame(sapply(BAFLevels, function(x){
  as.character(x)
}))
print(paste("Writing BAF Level segments to ",output_folder,"/",output_file_name,".csv",sep = ""))
write.table(x = BAFLevels_outfile, file = paste(output_folder,"/",output_file_name,".csv", sep=""),sep = "\t",
            col.names = TRUE,row.names = FALSE,quote = FALSE)
print(paste("Writing BAF Level segments to ",output_folder,"/",output_file_name,".csv"," (Done)",sep = ""))
