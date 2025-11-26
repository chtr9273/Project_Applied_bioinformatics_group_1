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

#Create an empty list with one element per tile
boc_list <- vector("list", N_tiles)

for(i in 1:N_tiles) {
  
  #positions in this tile
  df_loc <- df[df$Pos >= tiles[i] & df$Pos < tiles[i+1], ]
  
  if (nrow(df_loc) > 0) {
    
    #breadth of coverage for the specific tile
    boc_value <- sum(df_loc$N_reads > 0) / nrow(df_loc)
    
    #store the value for each position in each tile
    boc_list[[i]] <- rep(boc_value, nrow(df_loc))
    
  } else {
    #Store an empty vector if this tile contains no positions at all
    boc_list[[i]] <- numeric(0)
  }
}

#flaten boc_list
boc <- unlist(boc_list)


#make sure boc has the same number of elements as the rows in the dataframe
len_diff <- nrow(df) - length(boc)
if (len_diff > 0) {
  boc <- c(boc, rep(0, len_diff))}


# add boc to the original data
df$boc <- boc


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