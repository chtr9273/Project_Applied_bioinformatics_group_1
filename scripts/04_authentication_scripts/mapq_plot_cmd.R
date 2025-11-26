# define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#add arguments to respective variables
input_path <- args[1]
output_path <-args[2]

mapq_vals <- as.numeric(readLines(input_path))

png(filename = output_path, width = 800, height = 600, res = 150, type = "cairo")
hist(mapq_vals,
     col = 'darkred',
     breaks = 100, 
     main = "Histogram of MAPQ Values",
     xlab = "MAPQ value",
     ylab = "Frequency")
dev.off()