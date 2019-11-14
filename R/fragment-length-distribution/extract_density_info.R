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
.get_insert_size <- function (x){
  insert_size = abs(x$isize)
}

args = commandArgs(T)
if (length(args) < 2)
  stop("Usage: Rscript extract_density_info.R *runtime_variable_file(runtime_variables.R) *config_file ")
runtime_var_file = args[1]
config_file = args[2]

if(!.file_exists(c(runtime_var_file,config_file))) q('no')
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


sliding_windows = main_function.get_sliding_windows(binsize)
paste0("Read runtime variable file ",runtime_var_file)
source(runtime_var_file)
paste0("Read alignment info file ",readbam_rds_file)
reads_extracted_df = readRDS(readbam_rds_file)
output_file_name = paste(tools::file_path_sans_ext(basename(readbam_rds_file)),
                         "_insertsize_density_info",sep="")




output_folder_temp = paste0(output_folder,"/",extract_density_outdir_prefix)

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


nfragment = sapply(reads_extracted_df, function(sliding_window){
  isize = sliding_window$isize
  length(isize[which(isize >= short_length_range[1] & isize <= long_length_range[2])])
})

short = sapply(reads_extracted_df, function(sliding_window){
  isize = sliding_window$isize
  length(isize[which(isize >= short_length_range[1] & isize <= short_length_range[2])])
})
long = sapply(reads_extracted_df, function(sliding_window){
  isize = sliding_window$isize
  length(isize[which(isize >= long_length_range[1] & isize <= long_length_range[2])])
})
total.corrected=main_functions.bias_correct(nfragment, sliding_windows$gc/100)
short.corrected=main_functions.bias_correct(short, sliding_windows$gc/100)
long.corrected=main_functions.bias_correct(long, sliding_windows$gc/100)


insertsize_df = as.numeric(unlist(sapply(reads_extracted_df,.get_insert_size)))
insertsize_df_inrange = insertsize_df[which(insertsize_df <= maximum_length & insertsize_df >= short_length_range[1])]

insert_short_fragments = insertsize_df_inrange[which(insertsize_df_inrange >= short_length_range[1] & insertsize_df_inrange <= short_length_range[2])]
insert_long_fragments = insertsize_df_inrange[which(insertsize_df_inrange >= long_length_range[1] & insertsize_df_inrange <= long_length_range[2])]

.getmode <- function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

insert_info = data.frame("Total Fragments" = length(insertsize_df),
                         "Read Pairs in range" = sum(nfragment),
                         "Read Pairs in range_corrected" = sum(total.corrected,na.rm=T),
                         "Mode" = .getmode(insertsize_df_inrange),
                         "Median" = median(insertsize_df_inrange, na.rm = TRUE),
                         "Mean" = round(mean(insertsize_df_inrange, na.rm = TRUE), 2),
                         "SD" = round(sd(insertsize_df, na.rm = TRUE), 2),
                         "short" = sum(short,na.rm=T),
                         "long" = sum(long,na.rm=T),
                         "short_corrected" = sum(short.corrected,na.rm=T),
                         "long_corrected" = sum(long.corrected,na.rm=T),
                         "S/L Ratio" = round(length(insert_short_fragments) / length(insert_long_fragments),2),
                         "S/L Ratio_corrected" = sum(short.corrected,na.rm=T) / sum(long.corrected,na.rm=T),
                         "Bin Size"=binsize)

d = density(na.omit(insertsize_df_inrange)) # returns the density data
d = d[c("x","y")]

print(paste("Writing density to ",output_folder,"/",output_file_name,"_density_df.tsv",sep = ""))
write.table(x = d, file = paste(output_folder,"/",output_file_name,"_density_df.tsv", sep=""),sep = "\t",
            col.names = TRUE,row.names = FALSE,quote = FALSE)
print(paste("Writing density to ",output_folder,"/",output_file_name,"_density_df.tsv"," (Done)",sep = ""))

print(paste("Writing isize info to ",output_folder,"/",output_file_name,".tsv",sep = ""))
write.table(x = insert_info, file = paste(output_folder,"/",output_file_name,".tsv", sep=""),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
print(paste("Writing isize info to ",output_folder,"/",output_file_name,".tsv (Done)",sep = ""))


cat(paste0("density_file = \"",output_folder,"/",output_file_name,"_density_df.tsv","\"\n"),file=runtime_var_file,append = TRUE)
