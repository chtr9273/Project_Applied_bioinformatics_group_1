#define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#add arguments to respective variables
input_path <- args[1]
output_path <-args[2]

#Import output of samtools depth command
df <- read.delim(input_path, header = FALSE, sep = "\t")
names(df) <- c("Ref", "Pos", "N_reads")

#Split reference genome in tiles, compute breadth of coverage for each tile
N_tiles <- 500

#define window width
start <- min(df$Pos)
end   <- max(df$Pos)

#Setup tile borders
tiles <- seq(start, end, length.out = N_tiles + 1) 

#Create numeric vector for BoC, aligned to positions
boc_list <- vector("numeric", nrow(df))

for(i in 1:N_tiles) {
  
  #indices of positions in this tile
  idx <- which(df$Pos >= tiles[i] & df$Pos <= tiles[i+1])
  
  if(length(idx) > 0) {
    
    #breadth of coverage for the specific tile
    boc_value <- sum(df$N_reads[idx] > 0) / length(idx)
    
    #store the value for each position in this tile
    boc_list[idx] <- boc_value
    
  }
}

# add boc to the original data
df$boc <- boc_list

png(filename = output_path, width = 1600, height = 1200, res = 300, type = "cairo")

#plot boc agianst postions
plot(df$Pos, df$boc, type = "s", xlab = "Genome position", ylab = "Coverage",
     main = "Breadth of Coverage Across Genome")

# add abline to 0% coverage
abline(h = 0, col = "red", lty = 2)

#add % coverage in plot
mtext(paste0(round((sum(df$N_reads > 0) / length(df$N_reads)) * 100, 2), 
             "% of genome covered"), cex = 0.8)

dev.off()
