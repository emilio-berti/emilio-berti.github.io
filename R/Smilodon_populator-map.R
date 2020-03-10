library(tidyverse)
library(raster)
library(sf)


d <- read_rds("Smilodon_data.rds") #data is derived from the PHYLACINE database.

ggplot() +
  geom_sf(data = d[[2]], fill = "gainsboro") +
  geom_sf(data = d[[1]], fill = "green3", alpha = 0.5) +
  coord_sf(xlim = c(-15721463, -3467427),
           ylim = c(247164, 7327763)) +
  ggtitle(substitute(paste("Potential distribution of ", italic("Smilodon fatalis")))) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("../Figures/Smilodon_fatalis.png", width = 10, height = 6)
