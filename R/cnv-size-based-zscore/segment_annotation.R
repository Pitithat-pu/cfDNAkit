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
  stop("Usage: segment_annotation.R *runtime_variable_file(runtime_variables.R) *config_file")
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

paste0("Read runtime variable file ",runtime_var_file)
source(runtime_var_file)
source(paste0(package_dir,"/kit_functions/gene_annotation.R"))

z_score_segment_df = read.delim(z_score_df_file,sep="\t",stringsAsFactors = FALSE)
annotated_z_score_segment_df=gene_annotation.annotate_segment(z_score_segment_df)

paste0("Writing gene annotation result to ",z_score_df_file)
write.table(annotated_z_score_segment_df,
            z_score_df_file,quote = FALSE,
            row.names = FALSE,sep = "\t")
paste0("Done segment annotation")
