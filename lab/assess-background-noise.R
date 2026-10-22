library(ggplot2)
library(tidyr)

BACKGROUND <- 0.03442129 # mm
MASS <- 5 # kg

shakers <- data.frame(
  model = c("4810", "4809", "4808"),
  force = c(10, 45, 112)
)

# Peak displacement in mm
peak_displacement <- function(
  frequency,
  force,
  mass
) {
  acceleration = force / mass
  displacement <- acceleration / (2 * pi * frequency) ^ 2
  return(displacement * 1e3)
}

freq <- seq(0, 100, by = 5)
freq <- freq[freq > 0]

d <- lapply(shakers$force, \(force) {
  peak_displacement(freq, force, mass = MASS)
})
names(d) <- shakers$model
d <- as.data.frame(d)
d$frequency <- freq

# figure ------------
d |> 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "model",
    values_to = "peakDisplacement"
) |> 
  ggplot() +
    aes(x = frequency, y = peakDisplacement, col = model, fill = model) +
    geom_hline(yintercept = BACKGROUND, linetype = "dashed") +
    geom_hline(yintercept = BACKGROUND * 10, linetype = "dotted") +
    geom_line() +
    geom_point(shape = 21) +
    scale_x_continuous("Frequency (Hz)", breaks = seq(0, 200, by = 10)) +
    scale_y_log10("Peak displacement (mm)") +
    scale_color_manual("Model", values = c("dodgerblue3", "indianred", "goldenrod")) +
    scale_fill_manual("Model", values = c("dodgerblue3", "indianred", "goldenrod")) +
    ggtitle(sprintf("Peak displacement for a weight of %d kg", MASS)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    )
ggsave("hbk.pdf", width = 5, height = 3)
