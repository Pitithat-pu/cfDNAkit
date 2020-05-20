list.of.packages <- c("ensembldb","EnsDb.Hsapiens.v75")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
library(ensembldb)
library(EnsDb.Hsapiens.v75)


# segment_file = "/icgc/dkfzlsdf/analysis/hipo2/hipo_K34R/fragment_length_analysis/K34R-H7ARKZ/cfDNAkit/tumor15-b1_K34R-H7ARKZ_merged.mdup/coverage_zscore_PON_100k/tumor15-b1_K34R-H7ARKZ_merged.mdup_100_readbam_coverage_zscore_reads_extracted_info.csv"
# segment_file = ""
# segment_df=read.table(segment_file,header = TRUE,stringsAsFactors = FALSE)
edb <- EnsDb.Hsapiens.v75
.annotate_by_region = function(segment_chr,segment_start,segment_end){
  # edb <- EnsDb.Hsapiens.v75
  grf <- GRangesFilter(GRanges(segment_chr,
                               ranges = IRanges(segment_start , segment_end)), type = "any")
  gn <- genes(edb,
              filter = ~ grf & (gene_biotype %in% c("protein_coding",
                                                    "miRNA",
                                                    "non_coding")),
              columns = c("gene_name", "entrezid", "gene_biotype"),
              return.type = "DataFrame")
  return(gn)
}

segment_annotation.annotate_segment = function(segment_df){
  query_results = apply(segment_df, 1,function(segment){
    .annotate_by_region(segment["Chromosome"],segment["Start"] , segment["End"])
  })
  gene_name= c()
  gene_id=c()
  for (query_result in query_results) {
    gene_name = c(gene_name, paste(query_result$gene_name,collapse = ","))
    gene_id = c(gene_id, paste(query_result$gene_id,collapse = ","))
    # print(query_result$gene_name)
  }
  segment_df$gene_name = gene_name
  segment_df$gene_id = gene_id
  return(segment_df)
}

