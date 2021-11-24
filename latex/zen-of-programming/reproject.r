library(raster)

# this raster is the present-natural distribution
# of the African lion.
r <- raster("lion.gri")
crs(r)

reproj <- projectRaster(r, crs = "+proj=moll +datum=WGS84 +units=m +no_defs", method = "ngb")

r[r == 0] <- NA
r <- trim(r)
reproj[reproj == 0] <- NA
reproj <- trim(reproj)

png("lion.png", width = 600, height = 600)
plot(r, col = "grey20", box = FALSE, axes = FALSE, legend = FALSE)
plot(reproj, add = TRUE, col = adjustcolor("tomato", alpha.f = 0.75), legend = FALSE)
dev.off()