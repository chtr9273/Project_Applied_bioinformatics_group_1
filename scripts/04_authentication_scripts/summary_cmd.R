#define arguments from command line
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_depth_file> <output_png>")
}

#add arguments to respective variables
input_path <- args[1]
output_path <- args[2]

library(magick)
library(grid)

# Your 4 PNG files
png_files <- c("boc_plot.png", "mapq_plot.png", "edit_distance_plot.png", "percent_identity_plot.png")
png_files <- file.path(input_path, png_files)

# Read images with magick
images <- lapply(png_files, image_read)
rasters <- lapply(images, as.raster)

# Grid settings
nrow_grid <- 2
ncol_grid <- 2
total_cells <- nrow_grid * ncol_grid

# Fill remaining cells with NULL
if(length(rasters) < total_cells) {
  rasters <- c(rasters, rep(list(NULL), total_cells - length(rasters)))
}

# Save grid as a PNG
png(output_path, width = 2000, height = 1200, res = 600, type = "cairo")  # adjust size as needed
grid.newpage()
#grid.text(plot_title,
 #         y = unit(0.98, "npc"),
  #        gp = gpar(fontsize = 20, fontface = "bold"))
pushViewport(viewport(layout = grid.layout(nrow_grid, ncol_grid)))

# Function to place image
place_image <- function(raster_img, row, col) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  if(!is.null(raster_img)) grid.raster(raster_img, interpolate = TRUE)
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 2))
  popViewport()
}

# Loop over grid and place images
counter <- 1
for(i in 1:nrow_grid){
  for(j in 1:ncol_grid){
    place_image(rasters[[counter]], i, j)
    counter <- counter + 1
  }
}

dev.off()  # Close PNG device
