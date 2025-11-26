#define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#add arguments to respective variables
input_path <- args[1]
output_path <-args[2]

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
  #return NA if md_string is missing
  if (is.na(md_string)){
    return(NA)}
  #find all matches to a pattern, "[A-Za-z]" means all lower- and uppercase letters 
  loc <- gregexpr("[A-Za-z]", md_string)[[1]]
  #return 0 if no mismatches are found
  if (loc[1] == -1){
    return(0)}
  #return the length of loc, which corresponds to number of mismatches
  return(length(loc))
}

#define function header
parse_cigar <- function(cigar_string) {
  #separate various cigar operations
  #find number-letter pairs in cigar string
  matches <- gregexpr("[0-9]+[A-Z]", cigar_string)
  #extract the actual pairs in matches
  ops <- regmatches(cigar_string, matches)[[1]]
  
  #extract numbers from operation as numerical vector
  nums <- as.numeric(gsub("[A-Z]", "", ops))
  #extract letters from operations in a vector
  letters <- gsub("[0-9]", "", ops)
  #calculate number of aligned bases (matches or mismatches: M, =, X)
  matches <- sum(nums[letters %in% c("M", "=", "X")])
  #calculate number of inserts("I")
  inserts <- sum(nums[letters == "I"])
  #calculate number of deletions("D")
  deletes <- sum(nums[letters == "D"])
  #return a numeric list with matches, inserts and deletions from cigar string
  return(list(matches = matches, inserts = inserts, deletes = deletes))
}

#compute percent identity

#create vector for storing percent identity values for each read
#the vector is numeric with zeros, and number of elements corresponding to the
#initial number of reads
percent_identity <- numeric(length(cigar))

#loop over all reads
for (i in seq_along(cigar)) {
  
  #extract all aligned bases, inserts and deletions with parse_cigar function
  cg <- parse_cigar(cigar[i])
  
  #extract mismatches with count_mismatches_md function
  mm <- count_mismatches_md(md[i])
  
  #calculate the total alignment length
  total_aligned <- cg$matches + cg$inserts + cg$deletes
  
  #calculate number of correctly aligned bases
  correct_bases <- cg$matches - mm
  
  #calculate percentage identity for each read
  percent_identity[i] <- (correct_bases / total_aligned) * 100
}

#create barplot

#round down to closest integer
perc_id_round <- floor(percent_identity)

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

dev.off()
