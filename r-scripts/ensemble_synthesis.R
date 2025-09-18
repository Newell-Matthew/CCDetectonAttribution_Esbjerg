
average_ensemble <- function(data) {
  # Simple mean of estimates
  s1 <- mean(data$est)
  
  # Estimate average standard errors from the CI bounds
  se_lower <- mean((data$est - data$lower) / 1.96)
  se_upper <- mean((data$upper - data$est) / 1.96)
  
  # Compute symmetric 95% CI using averaged SEs
  lower <- s1 - 1.96 * se_lower
  upper <- s1 + 1.96 * se_upper
  
  return(setNames(c(s1, lower, upper), c("est", "lower", "upper")))
}


ensemble_synthesis <- function (obs_in = NA, models_in, synth_type = "abs") 
{
  if (is.na(unlist(obs_in))[1]) {
    no_obs <- T
    obs_in <- data.frame(est = 0, lower = 0, upper = 0)
    rownames(obs_in) <- "dummy"
  }
  else {
    no_obs <- F
  }
  colnames(obs_in) <- colnames(models_in) <- c("est", "lower", 
                                               "upper")
  if (!("model" %in% colnames(obs_in))) 
    obs_in$model <- rownames(obs_in)
  if (!("model" %in% colnames(models_in))) 
    models_in$model <- rownames(models_in)
  if (!synth_type %in% c("abs", "rel", "PR")) {
    cat(paste0("Synthesis type '", synth_type, "' not implemented - must be abs, rel or PR"), 
        "\n")
  }
  if (synth_type == "PR") {
    obs_in[, c("est", "lower", "upper")] <- log(obs_in[, 
                                                       c("est", "lower", "upper")])
    models_in[, c("est", "lower", "upper")] <- log(models_in[, 
                                                             c("est", "lower", "upper")])
  }
  else if (synth_type == "rel") {
    obs_in[, c("est", "lower", "upper")] <- log(1 + obs_in[, 
                                                           c("est", "lower", "upper")]/100)
    models_in[, c("est", "lower", "upper")] <- log(1 + models_in[, 
                                                                 c("est", "lower", "upper")]/100)
  }
  nobs = nrow(obs_in)
  obs <- apply(obs_in[, c("est", "lower", "upper"), drop = F], 
               2, mean)
  if (nobs == 1) {
    sig_obs = 0
  }
  else {
    s2 = sum((obs_in$est - obs[1])^2)
    sig_obs = sqrt(s2/(nobs - 1))
  }
  obs_in$l_wb <- obs_in$est - sqrt((obs_in$est - obs_in$lower)^2 + 
                                     (1.96 * sig_obs)^2)
  obs_in$u_wb <- obs_in$est + sqrt((obs_in$est - obs_in$upper)^2 + 
                                     (1.96 * sig_obs)^2)
  obs[2] <- obs[1] - sqrt((obs[1] - obs[2])^2 + (1.96 * sig_obs)^2)
  obs[3] <- obs[1] + sqrt((obs[1] - obs[3])^2 + (1.96 * sig_obs)^2)
  chi2 <- getsynchi2(models_in, sig_mod = 0)
  mdof <- nrow(models_in) - 1
  if (chi2/mdof > 1) {
    sig_mod <- optim(0, function(x) {
      (getsynchi2(models_in, sig_mod = x) - (nrow(models_in) - 
                                               1))^2
    }, method = "Brent", lower = 0, upper = 5)$par
  }
  else {
    sig_mod <- 0
  }
  models <- average_ensemble(models_in)
  models_in$l_wb <- models_in$est - sqrt((models_in$est - 
                                            models_in$lower)^2 + (1.96 * sig_mod)^2)
  models_in$u_wb <- models_in$est + sqrt((models_in$est - 
                                            models_in$upper)^2 + (1.96 * sig_mod)^2)
  w_obs <- unname((obs["upper"] - obs["lower"])^{
    -2
  })
  w_mod <- unname((models["upper"] - models["lower"])^{
    -2
  })
  wmean <- (w_obs * obs["est"] + w_mod * models["est"])/(w_obs + 
                                                           w_mod)
  sig_lower = sqrt((w_obs * ((obs["est"] - obs["lower"])/1.96)^2 + 
                      w_mod * ((models["est"] - models["lower"])/1.96)^2)/(w_obs + 
                                                                             w_mod))
  sig_upper = sqrt((w_obs * ((obs["est"] - obs["upper"])/1.96)^2 + 
                      w_mod * ((models["est"] - models["upper"])/1.96)^2)/(w_obs + 
                                                                             w_mod))
  synth <- setNames(c(wmean, wmean - 1.96 * sig_lower, wmean + 
                        1.96 * sig_upper), c("est", "lower", "upper"))
  umean <- (obs["est"] + models["est"])/2
  synth["l_wb"] <- umean - sqrt(((obs["est"] - obs["lower"])^2 + 
                                   (models["est"] - models["lower"])^2)/2)
  synth["u_wb"] <- umean + sqrt(((obs["est"] - obs["upper"])^2 + 
                                   (models["est"] - models["upper"])^2)/2)
  obs_in <- cbind(obs_in, group = "obs")
  obs <- data.frame(t(c(model = "Observations", group = "obs_synth", 
                        obs)))
  models_in <- cbind(models_in, group = "models")
  models <- data.frame(t(c(model = "Models", group = "model_synth", 
                           models)))
  synth <- data.frame(t(c(model = "Synthesis", group = "synth", 
                          synth)))
  res <- rbind.fill(obs_in, obs, models_in, models, synth)[, 
                                                           c("group", "model", "est", "lower", "upper", "l_wb", 
                                                             "u_wb")]
  for (cnm in c("est", "lower", "upper", "l_wb", "u_wb")) {
    res[, cnm] <- as.numeric(res[, cnm])
  }
  if (no_obs) {
    res <- res[grepl("model", res$group), ]
    sig_obs <- NA
  }
  if (synth_type == "PR") {
    res[, c("est", "lower", "upper", "l_wb", "u_wb")] <- exp(res[, 
                                                                 c("est", "lower", "upper", "l_wb", "u_wb")])
    sig_obs <- exp(sig_obs)
    sig_mod <- exp(sig_mod)
    umean <- exp(umean)
  }
  else if (synth_type == "rel") {
    res[, c("est", "lower", "upper", "l_wb", "u_wb")] <- 100 * 
      (exp(res[, c("est", "lower", "upper", "l_wb", "u_wb")]) - 
         1)
    sig_obs <- 100 * (exp(sig_obs) - 1)
    sig_mod <- 100 * (exp(sig_mod) - 1)
    umean <- 100 * (exp(umean) - 1)
  }
  return(list(synth_type = synth_type, sig_obs = sig_obs, 
              `chi2/dof` = chi2/mdof, sig_mod = sig_mod, df = res, 
              uw_mean = umean))
}

