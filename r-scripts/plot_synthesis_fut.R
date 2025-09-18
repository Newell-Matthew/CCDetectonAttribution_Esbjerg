plot_synthesis_fut <- function(synth, xlim, lwd = 10, xlab = "", main = "", add_space = TRUE, 
                                  log = NA, hide_labels = FALSE) {
  # Define colors (only models and synth)
  gcols = c(
    models = adjustcolor("orange", 0.75),
    ens = adjustcolor("red", 0.5),
    synth = "red"
  )
  
  # Handle log axis
  if (is.na(log)) {
    if (!is.null(synth$synth)) {
      logaxs <- if (synth$synth_type == "PR") "x" else ""
    } else {
      logaxs <- ""
    }
  } else {
    logaxs <- if (log) "x" else ""
  }
  
  # If synth is a list, extract the dataframe
  if (is(synth, "list")) 
    synth <- synth$df
  
  # Set xlim if missing
  if (missing(xlim)) {
    if (logaxs == "x") {
      xlim <- exp(range(pretty(log(as.numeric(unlist(synth[, c("lower", "upper", "l_wb", "u_wb")]))))))
    } else {
      xlim <- range(pretty(as.numeric(unlist(synth[, c("lower", "upper", "l_wb", "u_wb")]))))
    }
  }
  
  # Make sure group labels are correct
  if (is.numeric(synth$group)) 
    synth$group <- names(gcols)[synth$group]
  
  # Only count models now
  nmod <- sum(synth$group == "models")
  if (add_space) {
    spacing <- 2  # space between groups
    yy <- c()
    offset <- 0
    
    for (g in names(gcols)) {
      n <- sum(synth$group == g)
      yy <- c(yy, seq_len(n) + offset)
      offset <- max(yy) + spacing
    }
    
    yy <- max(yy) - yy + 1
  } else {
    yy <- nrow(synth):1
  }
  
  # Set vertical reference line
  vline <- if (logaxs == "x") 1 else 0
  
  # Plot empty canvas
  plot(0, type = "n", xlim = xlim, ylim = range(yy) + c(-0.5, 0.5),
       log = logaxs, yaxt = "n", ylab = "", xlab = xlab, main = main)
  grid(ny = NA, col = adjustcolor("black", 0.1), lty = 1)
  abline(v = vline, lty = 2)
  
  # Match colors
  gcols <- gcols[synth$group]
  
  # Plot white background segments
  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd, col = "black", lend = 2)
  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd - 2, col = "white", lend = 2)
  
  # Plot colored uncertainty bars
  segments(y0 = yy, x0 = synth$lower, x1 = synth$upper, lwd = lwd, col = gcols, lend = 2)
  
  # Plot points
  points(synth$est, yy, pch = 21, bg = gcols, lwd = 2, cex = lwd / 10)
  
  # Add labels if not hidden
  if (!hide_labels) 
    axis(2, at = yy, labels = synth$model, las = 1)
}