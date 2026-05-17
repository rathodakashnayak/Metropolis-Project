# ============================================================
# Metropolis Algorithm Simulation
# Group 11 — Metropolis Algorithm
# Indian Statistical Institute
# May 2026
# ============================================================

set.seed(42)

# ============================================================
# SECTION 1: GENERAL METROPOLIS FUNCTION
# ============================================================

metropolis <- function(S, n_iter, x_init, c_step) {
  x <- numeric(n_iter)
  x[1] <- x_init
  naccept <- 0
  
  for (k in 2:n_iter) {
    x_curr <- x[k-1]
    dx     <- runif(1, -c_step, c_step)
    x_prop <- x_curr + dx
    
    S_curr <- S(x_curr)
    S_prop <- S(x_prop)
    
    if (is.finite(S_prop)) {
      log_ratio <- S_curr - S_prop
      if (log(runif(1)) < log_ratio) {
        x[k]    <- x_prop
        naccept <- naccept + 1
      } else {
        x[k] <- x_curr
      }
    } else {
      x[k] <- x_curr
    }
  }
  
  list(chain           = x,
       acceptance_rate = naccept / n_iter)
}

# ============================================================
# SECTION 2: DEFINE DISTRIBUTIONS
# ============================================================

# Gaussian: S(x) = x^2/2
# True mean = 0, true <x^2> = 1
S_gaussian <- function(x) {
  0.5 * x^2
}

# Exponential: S(x) = x for x >= 0, Inf for x < 0
# True mean = 1, true <x^2> = 2
S_exponential <- function(x) {
  if (x < 0) return(Inf)
  return(x)
}

# ============================================================
# SECTION 3: GAUSSIAN HISTOGRAM — 5 PANELS
# Extended to K = 10^2 through 10^8 (Dipendu's suggestion)
# ============================================================

plot_histogram_convergence_5panel <- function(S_func,
                                              true_density,
                                              x_range,
                                              filename,
                                              c_step,
                                              x_init) {
  K_values <- c(1e2, 1e3, 1e5, 1e7, 1e8)
  
  png(filename, width = 2000, height = 450, res = 100)
  par(mfrow = c(1, 5), mar = c(4, 4, 3, 1))
  
  for (K in K_values) {
    cat("Running K =", format(K, scientific = TRUE), "\n")
    result  <- metropolis(S_func, K, x_init, c_step)
    burnin  <- ceiling(K * 0.1)
    x_post  <- result$chain[burnin:K]
    
    hist(x_post,
         breaks      = 60,
         probability = TRUE,
         col         = "steelblue",
         border      = "white",
         main        = bquote(K == .(format(K,
                                            scientific = TRUE,
                                            digits = 1))),
         xlab        = "x",
         ylab        = "Density",
         xlim        = x_range)
    
    x_seq <- seq(x_range[1], x_range[2], length.out = 500)
    lines(x_seq, true_density(x_seq), col = "red", lwd = 2)
    cat("Done K =", format(K, scientific = TRUE), "\n")
  }
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_histogram_convergence_5panel(
  S_func       = S_gaussian,
  true_density = dnorm,
  x_range      = c(-4, 4),
  filename     = "plot1a_gaussian_histogram_5panel.png",
  c_step       = 2.0,
  x_init       = 0
)

# ============================================================
# SECTION 4: EXPONENTIAL HISTOGRAM — 4 PANELS
# Idea: Sharath Kumar Yerramada
# ============================================================

plot_histogram_convergence_4panel <- function(S_func,
                                              true_density,
                                              x_range,
                                              filename,
                                              c_step,
                                              x_init) {
  K_values <- c(1e2, 1e3, 1e5, 1e7)
  
  png(filename, width = 1600, height = 450, res = 100)
  par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
  
  for (K in K_values) {
    cat("Running K =", format(K, scientific = TRUE), "\n")
    result  <- metropolis(S_func, K, x_init, c_step)
    burnin  <- ceiling(K * 0.1)
    x_post  <- result$chain[burnin:K]
    
    hist(x_post,
         breaks      = 50,
         probability = TRUE,
         col         = "steelblue",
         border      = "white",
         main        = bquote(K == .(format(K,
                                            scientific = TRUE,
                                            digits = 1))),
         xlab        = "x",
         ylab        = "Density",
         xlim        = x_range)
    
    x_seq <- seq(x_range[1], x_range[2], length.out = 500)
    lines(x_seq, true_density(x_seq), col = "red", lwd = 2)
    cat("Done K =", format(K, scientific = TRUE), "\n")
  }
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_histogram_convergence_4panel(
  S_func       = S_exponential,
  true_density = dexp,
  x_range      = c(0, 7),
  filename     = "plot1b_exponential_histogram.png",
  c_step       = 1.0,
  x_init       = 1
)

# ============================================================
# SECTION 5: RUNNING MEAN CONVERGENCE
# ============================================================

plot_running_mean <- function(S_func,
                              true_mean,
                              true_second_moment,
                              dist_name,
                              filename,
                              c_step,
                              x_init,
                              n_iter = 1e6) {
  result <- metropolis(S_func, n_iter, x_init, c_step)
  x      <- result$chain
  
  K_seq      <- unique(round(exp(
    seq(log(10), log(n_iter), length.out = 600)
  )))
  run_mean   <- sapply(K_seq, function(k) mean(x[1:k]))
  run_second <- sapply(K_seq, function(k) mean(x[1:k]^2))
  
  y_min <- min(c(run_mean, run_second,
                 true_mean, true_second_moment)) - 0.3
  y_max <- max(c(run_mean, run_second,
                 true_mean, true_second_moment)) + 0.3
  
  png(filename, width = 900, height = 500, res = 100)
  par(mar = c(5, 4, 3, 2))
  
  plot(K_seq, run_second,
       type = "l", col = "black", lwd = 1.5,
       log  = "x",
       ylim = c(y_min, y_max),
       xlab = "K",
       ylab = "Running average",
       main = paste(dist_name,
                    "— Convergence of Expected Values"))
  
  lines(K_seq, run_mean,
        col = "gray50", lwd = 1.5, lty = 2)
  
  abline(h = true_second_moment,
         col = "black",  lwd = 1, lty = 3)
  abline(h = true_mean,
         col = "gray50", lwd = 1, lty = 3)
  
  legend("bottomright",
         legend = c(expression(langle*x^2*rangle),
                    expression(langle*x*rangle)),
         col    = c("black", "gray50"),
         lty    = c(1, 2),
         lwd    = 2,
         bg     = "white")
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_running_mean(
  S_func             = S_gaussian,
  true_mean          = 0,
  true_second_moment = 1,
  dist_name          = "Gaussian",
  filename           = "plot2a_gaussian_running_mean.png",
  c_step             = 2.0,
  x_init             = 0
)

plot_running_mean(
  S_func             = S_exponential,
  true_mean          = 1,
  true_second_moment = 2,
  dist_name          = "Exponential",
  filename           = "plot2b_exponential_running_mean.png",
  c_step             = 1.0,
  x_init             = 1
)

# ============================================================
# SECTION 6: BURN-IN CHAIN HISTORY
# ============================================================

plot_burnin <- function(S_func,
                        dist_name,
                        filename,
                        c_step,
                        x_init_far,
                        true_typical,
                        n_iter = 5000) {
  result <- metropolis(S_func, n_iter, x_init_far, c_step)
  x      <- result$chain
  
  png(filename, width = 900, height = 450, res = 100)
  par(mar = c(5, 4, 3, 2))
  
  plot(1:n_iter, x,
       type = "l",
       col  = "steelblue",
       xlab = "k",
       ylab = "x",
       main = paste0(dist_name,
                     " — Chain History from x(0) = ",
                     x_init_far))
  
  abline(h   = true_typical,
         lty = 2,
         col = "red",
         lwd = 1.5)
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_burnin(
  S_func       = S_gaussian,
  dist_name    = "Gaussian",
  filename     = "plot3a_gaussian_burnin.png",
  c_step       = 2.0,
  x_init_far   = 20,
  true_typical = 0
)

plot_burnin(
  S_func       = S_exponential,
  dist_name    = "Exponential",
  filename     = "plot3b_exponential_burnin.png",
  c_step       = 1.0,
  x_init_far   = 15,
  true_typical = 1
)

# ============================================================
# SECTION 7: AUTOCORRELATION FUNCTION
# ============================================================

plot_acf <- function(S_func,
                     dist_name,
                     filename,
                     c_step,
                     x_init,
                     n_iter      = 1e5,
                     burnin_frac = 0.1) {
  result <- metropolis(S_func, n_iter, x_init, c_step)
  burnin <- ceiling(n_iter * burnin_frac)
  x_post <- result$chain[burnin:n_iter]
  
  png(filename, width = 900, height = 450, res = 100)
  par(mar = c(5, 4, 3, 2))
  
  acf(x_post,
      lag.max = 100,
      col     = "steelblue",
      lwd     = 2,
      main    = paste0(dist_name,
                       " — Autocorrelation Function",
                       " (c = ", c_step, ")"))
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_acf(
  S_func    = S_gaussian,
  dist_name = "Gaussian",
  filename  = "plot4a_gaussian_acf.png",
  c_step    = 2.0,
  x_init    = 0
)

plot_acf(
  S_func    = S_exponential,
  dist_name = "Exponential",
  filename  = "plot4b_exponential_acf.png",
  c_step    = 1.0,
  x_init    = 1
)

# ============================================================
# SECTION 8: JACKKNIFE ERROR vs BLOCK SIZE
# ============================================================

jackknife_error <- function(x, w) {
  n_groups <- floor(length(x) / w)
  if (n_groups < 2) return(NA)
  
  group_means <- sapply(1:n_groups, function(l) {
    mean(x[((l-1)*w + 1):(l*w)])
  })
  
  f_bar <- mean(group_means)
  sqrt(sum((group_means - f_bar)^2) /
         (n_groups * (n_groups - 1)))
}

plot_jackknife <- function(S_func,
                           dist_name,
                           filename,
                           c_step,
                           x_init,
                           n_iter      = 1e5,
                           burnin_frac = 0.1) {
  result <- metropolis(S_func, n_iter, x_init, c_step)
  burnin <- ceiling(n_iter * burnin_frac)
  x_post <- result$chain[burnin:n_iter]
  fx     <- x_post^2
  
  w_values <- 1:150
  delta_w  <- sapply(w_values,
                     function(w) jackknife_error(fx, w))
  
  png(filename, width = 900, height = 450, res = 100)
  par(mar = c(5, 4, 3, 2))
  
  plot(w_values, delta_w,
       type = "l",
       col  = "steelblue",
       lwd  = 2,
       xlab = "Block size w",
       ylab = expression(Delta[w]),
       main = paste(dist_name,
                    "— Jackknife Error vs Block Size"))
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_jackknife(
  S_func    = S_gaussian,
  dist_name = "Gaussian",
  filename  = "plot5a_gaussian_jackknife.png",
  c_step    = 2.0,
  x_init    = 0
)

plot_jackknife(
  S_func    = S_exponential,
  dist_name = "Exponential",
  filename  = "plot5b_exponential_jackknife.png",
  c_step    = 1.0,
  x_init    = 1
)

# ============================================================
# SECTION 9: STEP SIZE COMPARISON
# ============================================================

plot_stepsize_comparison <- function(S_func,
                                     dist_name,
                                     filename,
                                     c_values,
                                     true_second_moment,
                                     x_init,
                                     n_iter = 2e5) {
  colors <- c("red", "blue", "darkgreen")
  
  K_seq <- unique(round(exp(
    seq(log(100), log(n_iter), length.out = 400)
  )))
  
  png(filename, width = 900, height = 500, res = 100)
  par(mar = c(5, 4, 3, 2))
  
  plot(NULL,
       xlim = c(100, n_iter),
       ylim = c(true_second_moment - 0.5,
                true_second_moment + 0.5),
       log  = "x",
       xlab = "K",
       ylab = expression(langle*x^2*rangle),
       main = paste(dist_name,
                    "— Effect of Step Size c"))
  
  abline(h   = true_second_moment,
         lty = 2,
         col = "gray50",
         lwd = 1.5)
  
  for (i in seq_along(c_values)) {
    result     <- metropolis(S_func, n_iter,
                             x_init, c_values[i])
    x          <- result$chain
    run_second <- sapply(K_seq,
                         function(k) mean(x[1:k]^2))
    lines(K_seq, run_second,
          col = colors[i], lwd = 1.5)
    
    ar <- round(result$acceptance_rate * 100, 1)
    cat(dist_name, "| c =", c_values[i],
        "| Acceptance rate:", ar, "%\n")
  }
  
  legend("bottomright",
         legend = paste("c =", c_values),
         col    = colors,
         lwd    = 2,
         bg     = "white")
  
  dev.off()
  cat("Saved:", filename, "\n")
}

plot_stepsize_comparison(
  S_func             = S_gaussian,
  dist_name          = "Gaussian",
  filename           = "plot6a_gaussian_stepsize.png",
  c_values           = c(0.5, 2.0, 8.0),
  true_second_moment = 1,
  x_init             = 0
)

plot_stepsize_comparison(
  S_func             = S_exponential,
  dist_name          = "Exponential",
  filename           = "plot6b_exponential_stepsize.png",
  c_values           = c(0.2, 1.0, 4.0),
  true_second_moment = 2,
  x_init             = 1
)

# ============================================================
# SECTION 10: ACCEPTANCE RATE TABLES
# ============================================================

cat("\n========================================\n")
cat("Acceptance Rate Table — Gaussian\n")
cat("========================================\n")
cat(sprintf("%-12s %-18s %-15s\n",
            "Step size c", "Acceptance rate", "c x A.R."))
cat(rep("-", 48), "\n", sep = "")

for (c_step in c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0)) {
  result <- metropolis(S_gaussian, 10000, 0, c_step)
  ar     <- result$acceptance_rate
  cat(sprintf("%-12.1f %-18.4f %-15.4f\n",
              c_step, ar, c_step * ar))
}

cat("\n========================================\n")
cat("Acceptance Rate Table — Exponential\n")
cat("========================================\n")
cat(sprintf("%-12s %-18s %-15s\n",
            "Step size c", "Acceptance rate", "c x A.R."))
cat(rep("-", 48), "\n", sep = "")

for (c_step in c(0.2, 0.5, 1.0, 2.0, 3.0, 4.0)) {
  result <- metropolis(S_exponential, 10000, 1, c_step)
  ar     <- result$acceptance_rate
  cat(sprintf("%-12.1f %-18.4f %-15.4f\n",
              c_step, ar, c_step * ar))
}

cat("\n========================================\n")
cat("All plots saved successfully.\n")
cat("========================================\n")