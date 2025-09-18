plot_returnlevels_fix <- function (mdl, cov_f, cov_cf, ev, seed = 42, nsamp = 1000, model_desc = T, 
                                   xlim = c(1, 10000), ylim = NA, pch = 20, xlab = "Return period (years)", 
                                   ylab = NA, main = "", legend_pos = "topright", legend_labels = c("Present climate", 
                                                                                                    "Counterfactual climate")) 
{
  x <- mdl$x
  if (missing(ev)) {
    ev <- mdl$ev
  }
  rp_x <- unique(c(seq(1.1, 2, 0.1), seq(2, 100, 1), seq(100, 
                                                         1000, 10), seq(100, 1000, 100), seq(1000, 10000, 1000)))
  rp_th <- 1/seq(1, 0, length.out = length(x) + 2)[2:(length(x) + 
                                                        1)]
  if (nrow(cov_f) > 1) {
    print("cov_f has more than one row: only first row will be used as factual covariates")
    cov_f <- cov_f[1, , drop = F]
  }
  if (nrow(cov_cf) > 1) {
    print("cov_cf has more than one row: only first row will be used as counterfactual covariates")
    cov_cf <- cov_cf[1, , drop = F]
  }
  rl_curve_pres <- map_from_u(mdl, 1/rp_x, fixed_cov = cov_f)
  rl_curve_cf <- map_from_u(mdl, 1/rp_x, fixed_cov = cov_cf)
  rl_obs_pres <- map_from_u(mdl, map_to_u(mdl), fixed_cov = cov_f)
  rl_obs_cf <- map_from_u(mdl, map_to_u(mdl), fixed_cov = cov_cf)
  rp_event_pres <- 1/map_to_u(mdl, ev, fixed_cov = cov_f)
  rp_event_cf <- 1/map_to_u(mdl, ev, fixed_cov = cov_cf)
  if (is.na(ylim[1])) {
    ylim <- range(pretty(c(x, rl_curve_pres, rl_curve_cf)))
  }
  if (is.na(ylab)) {
    ylab <- mdl$varnm
  }
  plot(0, type = "n", xlim = xlim, ylim = ylim, log = "x", 
       xlab = "", ylab = "", main = main)
  mtext(xlab, side = 1, line = 2.5, cex = par("cex.lab"))
  mtext(ylab, side = 2, line = 2.5, cex = par("cex.lab"))
  if (model_desc) {
    legend_title <- paste0(mdl$varnm, " ~ ", paste0(mdl$covnm, 
                                                    collapse = " + "), " (", mdl$dist, ", ", mdl$type, 
                           ")")
  }
  else {
    legend_title <- ""
  }
  legend(legend_pos, legend = c(legend_labels, "Observed event"), 
         col = c("firebrick", "blue", "magenta"), lty = 1, pch = c(pch, 
                                                                   pch, NA), bty = "n", cex = par()$cex.lab, title = legend_title)
  lines(rp_x, rl_curve_pres, lwd = 2, col = "firebrick", lty = 1)
  lines(rp_x, rl_curve_cf, lwd = 2, col = "blue", lty = 1)
  points(rp_th, sort(rl_obs_pres, decreasing = mdl$lower), 
         col = "firebrick", pch = pch)
  points(rp_th, sort(rl_obs_cf, decreasing = mdl$lower), col = "blue", 
         pch = pch)
  abline(h = ev, col = "magenta", lty = 2)
  suppressWarnings(rug(rp_event_pres, lwd = 3, col = "firebrick"))
  suppressWarnings(rug(rp_event_cf, lwd = 3, col = "blue"))
  if (!is.na(nsamp)) {
    x_ci <- c(5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 
              10000)
    set.seed(seed)
    ctab = names(mdl$data)
    mdl_df <- mdl$data[, ctab]
    boot_res <- sapply(1:nsamp, function(i) {
      boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), 
                               replace = T), ]
      tryCatch({
        boot_mdl <- refit(mdl, boot_df)
        c(map_from_u(boot_mdl, 1/x_ci, fixed_cov = cov_f), 
          map_from_u(boot_mdl, 1/x_ci, fixed_cov = cov_cf))
      }, error = function(cond) {
        return(rep(NA, length(x_ci) * 2))
      })
    })
    est_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), 
                    na.rm = T)
    polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1, 1:length(x_ci)], 
                                          rev(est_ci[2, 1:length(x_ci)])), density = NULL, 
            border = NA, col = adjustcolor("firebrick", 0.1))
    polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1, -(1:length(x_ci))], 
                                          rev(est_ci[2, -(1:length(x_ci))])), density = NULL, 
            border = NA, col = adjustcolor("blue", 0.1))
  }
}

process_ci_df <- function(ci_df, invert = "PR", negate = c("dI_abs", "dI_rel")) {
  # Ensure it's a data frame
  ci_df <- as.data.frame(ci_df)
  
  # 1. Invert selected rows (e.g., PR) and reorder columns (est, upper, lower)
  for (row in invert) {
    if (row %in% rownames(ci_df)) {
      ci_df[row, ] <- 1 / ci_df[row, c("est", "97.5%", "2.5%")]
    }
  }
  
  # 2. Negate selected rows (e.g., dI_abs, dI_rel) and reorder columns
  for (row in negate) {
    if (row %in% rownames(ci_df)) {
      ci_df[row, ] <- -ci_df[row, c("est", "97.5%", "2.5%")]
    }
  }
  
  
  return(ci_df)
}

