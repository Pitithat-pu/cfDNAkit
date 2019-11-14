chrLength_file = "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_chrTotalLength.tsv"
control_density_file = "/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/BH01_rmdup_paired_mapped_Fragment-length_report_50_insert_size_density.csv"
duke_blacklist_region = "/abi/data2/puranach/ct_DNA/cfDNAkit/resources/wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz"
dac_blacklist_region = "/abi/data2/puranach/ct_DNA/cfDNAkit/resources/wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz"
# control_min_coverage = 100
# sampling_size = 20
flag=scanBamFlag(isPaired = TRUE, isProperPair=TRUE,
                 isUnmappedQuery = FALSE, isFirstMateRead = TRUE,
                 isDuplicate = FALSE,isSecondMateRead = FALSE,
                 hasUnmappedMate = FALSE, isMinusStrand = FALSE,
                 isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE)
maximum_length = 250
minimum_length = 1
min_coverage = 20
binsize = 100 #kb
cutoff_zscore = c(1.96,-1.96)
#cutoff_zscore = c(1.65,-1.65)
# cutoff_zscore = c(1.28, -1.28)
short_length_range = c(100,150)
long_length_range = c(151,maximum_length)
plot_upperbound = 15
plot_lowerbound = -15
min_heterozygous_SNPs = 5
blacklist_targets_gr = create_blacklist_gr(c(duke_blacklist_region,dac_blacklist_region))
delta_F_binsize = 5000
