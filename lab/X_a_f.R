a <- seq(0, 1, by = 0.1)
f <- seq(1, 3, length.out = length(a))
X <- matrix(NA, length(a), length(a))
for (i in seq_along(a)) {
  for (j in seq_along(f)) {
    X[i, j] <- a[i] / (2 * pi * f[j])^2
  }
}
image(
  a, f, X,
  xlab = expression(Acceleration (m/s^2)),
  ylab = expression(Frequency (Hz)),
  main = "Peak displacement (m)",
  col = hcl.colors(100, "Spectral", rev = TRUE)
)

png("~/emilio-berti.github.io/lab/peak-displacement.png", 1800, 1200, res = 300)
par(mar = c(4, 4, 2, 2))
for (i in seq_along(a)) {
  col <- hcl.colors(length(a), "Spectral", rev = TRUE)[i]
  if (i == 1) {
    plot(
      a, X[, i], pch = 21, bg = col,
      xlab = expression(PeakAcceleration (m/s^2)),
      ylab = expression(PeakDisplacement (m))
    )
  } else {
    points(a, X[, i], pch = 21, bg = col)
  }
  lines(a, X[, i], col = col, lwd = 2)
}
legend(
  x = 0, y = 0.025,
  legend = f[seq(1, length(f), by = 2)],
  fill = hcl.colors(length(a), "Spectral", rev = TRUE)[seq(1, length(f), by = 2)],
  box.lwd = 0,
  title = expression(Frequency (Hz))
)
dev.off()
