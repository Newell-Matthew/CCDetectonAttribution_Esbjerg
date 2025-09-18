plot_synthesis_mod <- function (synth, xlim, lwd = 10, xlab = "", main = "", add_space = T, 
          log = NA, hide_labels = F) 
{
  gcols = c(obs = adjustcolor("blue", 0.5), obs_synth = "blue", 
            models = adjustcolor("orange", 0.75), ens = adjustcolor("red", 0.5), ens_average = "red", 
            synth_model = adjustcolor("magenta", 0.25),
            synth_ens= "magenta")
  print(names(gcols))
  if (is.na(log)) {
    if (!is.null(synth$synth)) {
      if (synth$synth_type == "PR") {
        logaxs <- "x"
      }
      else {
        logaxs <- ""
      }
    }
    else {
      logaxs <- ""
    }
  }
  else {
    if (log) {
      logaxs <- "x"
    }
    else {
      logaxs <- ""
    }
  }
  if (is(synth, "list")) 
    synth <- synth$df
  if (missing(xlim)) {
    if (logaxs == "x") {
      xlim <- exp(range(pretty(log(as.numeric(unlist(synth[, 
                                                           c("lower", "upper", "l_wb", "u_wb")]))))))
    }
    else {
      xlim <- range(pretty(as.numeric(unlist(synth[, c("lower", 
                                                       "upper", "l_wb", "u_wb")]))))
    }
  }
  if (is.numeric(synth$group)) 
    synth$group <- names(gcols)[synth$group]
  nobs <- sum(synth$group == "obs")
  nmod <- sum(synth$group == "models")
  nens <- sum(synth$group == "ens")
  if (add_space) {
    spacing <- 2   # space between groups
    yy <- c()
    offset <- 0
    gnames <- names(gcols)
    
    for (i in seq_along(gnames)) {
      g <- gnames[i]
      n <- sum(synth$group == g)
      
      # If it's one of the last two groups, use the same offset
      if (i >= length(gnames) - 1) {
        yy <- c(yy, seq_len(n) + offset)
        # Don't update offset so they stack
      } else {
        yy <- c(yy, seq_len(n) + offset)
        offset <- max(yy) + spacing
      }
    }
    
    yy <- max(yy) - yy + 1
  } else {
    yy <- nrow(synth):1
  }
  if (logaxs == "x") {
    vline <- 1
  }
  else {
    vline <- 0
  }
  plot(0, type = "n", xlim = xlim, ylim = range(yy) + c(-0.5, 
                                                        0.5), log = logaxs, yaxt = "n", ylab = "", xlab = xlab, 
       main = main)
  grid(ny = NA, col = adjustcolor("black", 0.1), lty = 1)
  abline(v = vline, lty = 2)
  gcols <- gcols[synth$group]
  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd, 
           col = "black", lend = 2)
  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd - 
             2, col = "white", lend = 2)
  segments(y0 = yy, x0 = synth$lower, x1 = synth$upper, lwd = lwd, 
           col = gcols, lend = 2)
  points(synth$est, yy, pch = 21, bg = gcols, lwd = 2, cex = lwd/10)
  if (!hide_labels) 
    axis(2, at = yy, labels = synth$model, las = 1)
}
