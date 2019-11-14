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
  stop("Usage: Rscript plot_density_info.R *runtime_variable_file(runtime_variables.R) *config_file")
runtime_var_file = args[1]
config_file = args[2]

if(!.file_exists(c(runtime_var_file,config_file))) q('no')
paste("source(",config_file,")",sep="")
source(config_file)

paste0("Read runtime variable file ",runtime_var_file)
source(runtime_var_file)

output_folder = paste0(output_folder,"/",extract_density_outdir_prefix)
if(!dir.exists(output_folder)) {
  q(paste0("no output directory for density file was found! Something wrong.\n
    please check ",runtime_var_file, " for density_file or validity of ", output_folder))
}
if(!.file_exists(c(density_file))) q('no')

.plot.multi.dens <- function(s,labels){
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, data.frame(s[i])$x)
    junk.y = c(junk.y, data.frame(s[i])$y)
  }

  # xr <- range(junk.x)
  yr <- range(junk.y)
  legend_text= c("Healthy Control")
  plot(x = data.frame(s[1])$x,
       y= data.frame(s[1])$y, xlim = c(0,maximum_length),lwd=2, ylim = yr, main = "", type="l",lty=2, col=c("gray"),
       xlab="Insert Size (bases)",
       ylab="Density" )
  for(i in 2:length(s)) {
    lines(x = data.frame(s[i])$x,
          y= data.frame(s[i])$y, xlim = c(0,maximum_length), ylim = yr,
          col = i,type="l",pch=8,lwd=2,cex.lab=2,cex.axis=2)
    legend_text = c(legend_text,paste(sep="",labels[i]))
    legend(legend = legend_text, "topright", col=c(1:length(s)),cex=0.8,lty=1)
  }

  title(main = "Insert Size Distribution")
  # legend(legend = legend_text, "topright", col=c("black","red"), lty=2:1, cex=0.8)
  # legend(legend = legend_text, "topright", col=c("black","red"),lty=2:1,cex=0.6)
}

.plot.multi.dens_single <- function(s,labels){
  # xr <- range(junk.x)
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, data.frame(s[i])$x)
    junk.y = c(junk.y, data.frame(s[i])$y)
  }
  yr <- range(junk.y)
  plot(x = data.frame(s[1])$x,
       y= data.frame(s[1])$y, xlim = c(0,maximum_length), ylim = yr,
       main = "", type="l",lty=1,lwd=2,
       xlab="Insert Size (bases)",
       ylab="Density", col = c("red"),cex.lab=2,cex.axis=2)

  title(main = "Insert Size Distribution",cex.main=2)
}
d <- read.table(density_file, header = TRUE)
# isize_info_df = read.table(isize_info_file,header=TRUE)
output_plot_file <- paste(tools::file_path_sans_ext(basename(density_file)),
                          "_plot.png",sep="")

png(paste(output_folder,"/",output_plot_file,sep=""), width=900, height=520, units="px")
par(mar=c(4, 4, 2, 0.5)) # no spacing between CNV anf BAF plots
par(mfrow=c(1,1)) # CNV and BAF plots above each other

#par(mar=c(4, 4, 2, 0.5)) # no spacing between CNV anf BAF plots
#par(mfrow=c(2, 1)) # CNV and BAF plots above each other
if (plot_insertsize_density_with_PON == "TRUE") {
  d_control <- read.table(control_density_file, header = TRUE)
  colnames(d_control)=c("x","y")
  .plot.multi.dens( list(d_control,d), c("Healthy control","sample"))
} else {
  .plot.multi.dens_single(list(d), "sample")
}
paste("Produce density plot file ",output_folder,output_plot_file,sep="/")
dev.off()



cat(paste0("density_plot_file = \"",output_folder,"/",output_plot_file,"\"\n"),file=runtime_var_file,append = TRUE)
