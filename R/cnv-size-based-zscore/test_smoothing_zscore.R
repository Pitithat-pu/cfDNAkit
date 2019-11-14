smooth_moving_window = 20 ## moving 20 windows at a time
cytoband_file = "/abi/data2/puranach/ct_DNA/cfDNAkit/resources/cytoBand_GRCh37.txt.gz"
get_ma <- function(x, window_size=smooth_moving_window){stats::filter(x, rep(1 / window_size, window_size), sides = 2)}

cytoband_df = read.table(cytoband_file,header = FALSE,stringsAsFactors = FALSE)
colnames(cytoband_df) = c("chr","start","end","cytoband","annotation")
chr_arm = sapply(strsplit(cytoband_df$cytoband,split=""),function(x){
  unlist(x)[1]
})
cytoband_df$chr_arm = paste(cytoband_df$chr,chr_arm,sep="")
cytoband_df_temp = cytoband_df
merge_cytoband_df = data.frame("chr"=as.character(),"start"=as.integer(),
                               "end"=as.integer(),"chr_arm"=as.character())
for (chr_arm in unique(cytoband_df_temp$chr_arm)) {
  cytoband_df_temp = cytoband_df[which(cytoband_df$chr_arm == chr_arm),
                                 c("chr","start","end","chr_arm" )]
  merge_cytoband_df = rbind(merge_cytoband_df,data.frame(cytoband_df_temp$chr[1],cytoband_df_temp$start[1],
             cytoband_df_temp$end[nrow(cytoband_df_temp)],chr_arm))
}

colnames(merge_cytoband_df)=c("chr","start","end","chr_arm" )
merge_cytoband_df$chr_arm = as.character(levels(merge_cytoband_df$chr_arm)[merge_cytoband_df$chr_arm])

test.smoothing_zscore.smooth_zscore = function(z_score_df){
  z_score_df$chr_arm  = unlist(apply(z_score_df,1, function(z_score_df_position){
    filter(merge_cytoband_df,chr==as.character(z_score_df_position["Chromosome"]),
           start<=as.integer(z_score_df_position["Start"]),
           end>=as.integer(z_score_df_position["End"])) %>% select(chr_arm)
  }))

  all_chr_arm = unique(z_score_df$chr_arm)

  smoothed_z_score = sapply(all_chr_arm, function(chr.arm){
    zscore = filter(z_score_df,chr_arm==chr.arm) %>% select(size_based_zscore)
    # print(chr.arm)
    if (length(unlist(zscore))>smooth_moving_window ) {
      as.numeric(get_ma(unlist(zscore)))
    } else as.numeric(unlist(zscore))

  })
  z_score_df$smoothed_zscore = unlist(smoothed_z_score)
  return(z_score_df)
}



