force_rating <- 10 # for Mini-shaker Type 4810

peak_displacement <- function(force_rating, mass, frequency) {
  peak_acceleration <- force_rating / mass
  ans <- peak_acceleration / (2 * pi * frequency) ^ 2
  return(ans)
}

frequency <- seq(10, 100, by = 10)  # Hz
mass <- seq(5, 20, by = 5)          # Kg
d <- expand.grid(frequency, mass)
names(d) <- c("frequency", "mass")
d$peakDisplacement <- peak_displacement(
  force_rating,
  d$mass,
  d$frequency
) * 1e3

# Figure -----------
op <- par(mar = c(5, 5, 2, 8) - 1, xpd = FALSE)
plot(
  NULL,
  xlab = "Mass (kg)",
  xlim = range(d$mass),
  ylab = ("Peak displacement (mm)"),
  ylim = c(1e-3, 1),
  log = "y",
  yaxt = "n",
  frame = FALSE
)
br <- 10 ^ seq(-3, 1)
axis(2, br, br)
pal <- hcl.colors(length(unique(d$frequency)), "Spectral", rev = TRUE)
for (freq in unique(d$frequency)) {
  color <- pal[which(unique(d$frequency) == freq)]
  with(
    subset(d, frequency == freq),
    lines(mass, peakDisplacement, col = color)
  )
  with(
    subset(d, frequency == freq),
    points(mass, peakDisplacement, pch = 21, bg = color)
  )
}
legend(
  "topright",
  title = "Frequency (Hz)",
  inset = c(-0.25, 0), xpd = TRUE, bty = "n",
  legend = unique(d$frequency),
  fill = pal
)
par(op)
