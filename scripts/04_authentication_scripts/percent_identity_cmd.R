#define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#arguments to respective variables
input_path <- args[1]
output_path <- args[2]

#load Rsamtools
library(Rsamtools)

#path to BAM
bam_file <- input_path 

#minimum MAPQ filter
min_mapq <- 0                  

#specify which parameters to extract
param <- ScanBamParam(
  what = c("qname", "flag", "mapq", "cigar"),
  tag = c("MD")
)

#read bam file with respect to specified parameters
bam <- scanBam(bam_file, param = param)[[1]]

#filter by MAPQ, only reads with positive or zero mapq values are kept
keep <- bam$mapq >= min_mapq
cigar <- bam$cigar[keep]
md    <- bam$tag$MD[keep]

#helper functions

#define function header
count_mismatches_md <- function(md_string) {
  if (is.na(md_string)) return(NA)
  # extract mismatch runs (letters only)
  mismatch_runs <- regmatches(md_string, gregexpr("[A-Za-z]+", md_string))[[1]]
  if (length(mismatch_runs) == 0) return(0)
  sum(nchar(mismatch_runs))
}

#define function header
parse_cigar <- function(cigar_string) {
  matches <- gregexpr("[0-9]+[A-Z]", cigar_string)
  ops <- regmatches(cigar_string, matches)[[1]]
  
  nums <- as.numeric(gsub("[A-Z]", "", ops))
  letters <- gsub("[0-9]", "", ops)
  
  matches_count <- sum(nums[letters %in% c("M", "=", "X")])
  inserts <- sum(nums[letters == "I"])
  deletes <- sum(nums[letters == "D"])
  
  return(list(matches = matches_count, inserts = inserts, deletes = deletes))
}

#compute percent identity
percent_identity <- numeric(length(cigar))

for (i in seq_along(cigar)) {
  cg <- parse_cigar(cigar[i])
  mm <- count_mismatches_md(md[i])
  
  total_aligned <- cg$matches + cg$deletes
  if (total_aligned == 0) {
    percent_identity[i] <- NA
    next
  }
  
  correct_bases <- cg$matches - mm
  percent_identity[i] <- (correct_bases / total_aligned) * 100
}

#create barplot

#round to nearest integer
perc_id_round <- round(percent_identity)

#remove NAs
perc_id_round <- perc_id_round[!is.na(perc_id_round)]

#create a frequency table
freq_table <- table(perc_id_round)

#plot
png(filename = output_path, width = 800, height = 600, res = 150, type = "cairo")

barplot(
  freq_table,
  col = "lightblue",
  border = "white",
  xlab = "Percent Identity to Reference (integer)",
  ylab = "Number of Reads",
  main = "Percent Identity Distribution of Aligned Reads",
  las = 2)

#average identity
mtext(paste0(round(mean(percent_identity, na.rm = TRUE), 2), 
             "% average identity"), cex = 0.8)

dev.off()
