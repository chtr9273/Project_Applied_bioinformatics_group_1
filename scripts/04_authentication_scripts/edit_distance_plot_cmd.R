# define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#add arguments to respective variables
input_path <- args[1]
output_path <-args[2]

#loar Rsamtools library
library("Rsamtools")

#create parameter object for finding the number of mismatches for each read
param <- ScanBamParam(tag = "NM")

#identify parameter and read bam file
bam <- scanBam(input_path, param = param)

png(filename = output_path, width = 800, height = 600, res = 150, type = "cairo")

#create a bar plot for the number of mismatches
barplot(table(bam[[1]]$tag$NM), 
        ylab ="Number of reads",
        xlab ="Number of mismatches",
        main = "Edit distance")
dev.off()