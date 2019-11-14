# control_shortread_file = "/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/deltaf_whole_genome/BH01_rmdup_paired_mapped_downsampled_deltaf_table_n20_100_mincov1000_shortread.corrected.csv"
# control_longread_file = "/icgc/dkfzlsdf/analysis/G200/puranach/ctDNA/1_alignment/BH01/deltaf_whole_genome//BH01_rmdup_paired_mapped_downsampled_deltaf_table_n20_100_mincov1000_longread.corrected.csv"
control_shortread_df = read.table(control_shortread_file,
                                header = TRUE,stringsAsFactors = FALSE)
control_longread_df = read.table(control_longread_file,
                                  header = TRUE,stringsAsFactors = FALSE)
control_SLRatio_df = control_shortread_df / control_longread_df
all_region_name = rownames(control_shortread_df)
# control_SLRatio_df_no_NA = control_SLRatio_df[complete.cases(control_SLRatio_df),]
# test_regions_df_corrected
get_reference_region <- function(test_regions_df_corrected){
  test_regions_df_reference =  test_regions_df_corrected[which(startsWith(rownames(test_regions_df_corrected),"X")),]
  # test_regions_df_reference_non_NA = test_regions_df_reference[complete.cases(test_regions_df_reference),]
  return (rownames(test_regions_df_reference))
}
reference_regions = get_reference_region(test_regions_df_corrected)
control_shortread_reference_df = control_shortread_df[reference_regions,]
control_longread_reference_df = control_longread_df[reference_regions,]
reference_short_read = colSums(control_shortread_reference_df,na.rm = TRUE)
reference_long_read = colSums(control_longread_reference_df,na.rm = TRUE)
P150_reference =  reference_short_read / reference_long_read

delta_f_df.corrected = as.data.frame(t(sapply(all_region_name,function(region_name){
    SLRatio.corrected = control_SLRatio_df[region_name,]

    Delta_F = as.numeric(round(SLRatio.corrected - P150_reference,4))
}) ))


