library(dplyr)

top_freqs <- function(x, fs, n_peaks = 2, pad_factor = 4, min_sep = 1, t = NULL) {
  if (is.null(t)) t <- (seq_along(x) - 1) / fs

  x0 <- x - mean(x)
  n <- length(x0)
  w <- 0.5 - 0.5 * cos(2*pi*(0:(n-1)) / (n-1))       # Hann window
  xw <- x0 * w

  n_fft <- 2^ceiling(log2(n)) * pad_factor
  X <- fft(c(xw, rep(0, n_fft - n)))
  k <- 0:(n_fft %/% 2)
  freqs <- k * fs / n_fft
  mag <- Mod(X)[k + 1]

  # peak detection in the spectrum
  is_peak <- c(FALSE, diff(sign(diff(mag))) == -2, FALSE)
  pk <- which(is_peak)
  pk <- pk[order(mag[pk], decreasing = TRUE)]
  chosen <- c()
  for (i in pk) {
    if (all(abs(freqs[i] - freqs[chosen]) >= min_sep)) chosen <- c(chosen, i)
    if (length(chosen) == n_peaks) break
  }
  f_hat <- freqs[chosen]

  # least-squares fit: x ~ c0 + sum_i [a_i sin(2*pi*f_i*t) + b_i cos(2*pi*f_i*t)]
  D <- do.call(cbind, lapply(f_hat, function(f)
         cbind(sin(2*pi*f*t), cos(2*pi*f*t))))
  fit <- lm(x ~ D)
  cf <- coef(fit)[-1]
  a <- cf[seq(1, length(cf), by = 2)]
  b <- cf[seq(2, length(cf), by = 2)]
  amp   <- unname(sqrt(a^2 + b^2))
  
  # phase conventions
  phase_sin <- unname(atan2(b, a))      # A*sin(2*pi*f*t + phase_sin)
  phase_cos <- unname(atan2(-a, b))     # A*cos(2*pi*f*t + phase_cos)

  # each component separately, evaluated at the sample times
  comps <- lapply(seq_along(f_hat), function(i)
    a[i] * sin(2*pi*f_hat[i]*t) + b[i] * cos(2*pi*f_hat[i]*t))

  summary_df <- do.call(rbind, lapply(seq_along(f_hat), function(i) {
    ci <- comps[[i]]
    data.frame(
          freq          = f_hat[i],
          max_x         = max(ci),
          min_x         = min(ci),
          max_abs_x     = max(abs(ci)),
          time_of_max   = t[which.max(ci)],
          amplitude     = amp[i],
          phase_sin_rad = phase_sin[i],
          phase_sin_deg = phase_sin[i] * 180 / pi,
          phase_cos_rad = phase_cos[i],
          phase_cos_deg = phase_cos[i] * 180 / pi,
          rms           = sqrt(mean(ci^2))
        )
  }))

  rec <- fitted(fit)
  list(
    components = summary_df,
    component_signals = setNames(comps, sprintf("%.2f Hz", f_hat)),
    offset = unname(coef(fit)[1]),
    reconstruction_max = max(rec),
    observed_max = max(x),
    r2 = summary(fit)$r.squared
  )
}

units <- paste0("U", seq_len(6))
files <- paste0(units, ".csv")
dec <- as.list(rep(NA, length(units)))
names(dec) <- units

for (f in files) {

  message(f)

  # Load skipping comments
  d <- read.csv(f, skip = 7)
  # Remove first 50 records
  d <- tail(d, n = nrow(d) - 30)
  # Leave out last 50 records
  d <- head(d, n = nrow(d) - 30)

  # Vertical acceleration
  z <- d$az

  # Dominant frequency of `z`
  dt <- d$time[-1] - d$time[-length(d$time)]
  sampling_frequency <- 1 / median(dt)
  res <- top_freqs(z, sampling_frequency, n_peaks = 2, min_sep = 1)

  # Reconstruct fitted acceleration from the two components
  fitted_acc <- res$offset + Reduce(`+`, res$component_signals)
  residual <- z - fitted_acc
  old_par <- par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
  plot(
    d$time, z,
    col = "dodgerblue3",
    type = "l",
    xlab = "",
    ylab = "Acceleration",
    main = paste(f, "- measured acceleration")
  )
  lines(d$time, fitted_acc, col = "indianred")
  legend(
    "topright",
    legend = c("Measured", "Two-frequency fit"),
    fill = c("dodgerblue3", "indianred"),
    bty = "n"
  )
  plot(
    d$time, residual,
    col = "dodgerblue3",
    type = "l",
    xlab = "Time (s)",
    ylab = "Residual",
    main = "Residual: measured - fitted"
  )
  abline(h = 0, col = "indianred")
  par(old_par)

  x <- with(
    res$components,
    - amplitude[1] / (2 * pi * freq[1]) ^ 2 * sin(2 * pi * freq[1] * d$time + phase_sin_rad[1]) -
    amplitude[2] / (2 * pi * freq[2]) ^ 2 * sin(2 * pi * freq[2] * d$time + phase_sin_rad[2])
  )
  peak_to_peak_displacement <- max(x) - min(x)

  dec[[sub("[.]csv", "", f)]] <- data.frame(
    frequency = res$components$freq,
    peak_acceleration = res$components$max_abs_x,
    peak_to_peak_displacement = peak_to_peak_displacement * 1e3
  )

  null <- readline("Press Enter to continue")

}

d <- do.call("rbind", dec)
d$unit <- gsub("[.][0-9]", "", rownames(d))
rownames(d) <- NULL

# Select only largest dominant frequency
peak_displacement <- d |>
  group_by(unit) |>
  arrange(desc(frequency)) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(peak_displacement = peak_acceleration / (2 * pi * frequency) ^ 2 * 1e3) |>
  pull(peak_displacement) |>
  mean()
