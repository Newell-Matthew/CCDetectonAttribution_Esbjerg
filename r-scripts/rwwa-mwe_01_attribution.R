# install & load required packages
devtools::install_github("WorldWeatherAttribution/rwwa")
library(rwwa)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ANALYSIS OF A SINGLE TIME SERIES                                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Load & plot data                                                                  ####

# The model fitting function takes as its main input a data.frame with named columns
# containing the variable and covariates of interest, and a column 'year' 
# containing integer-valued years.

# load covariate data
gmst <- read.csv("\\..\\GMST_anomalies.csv", col.names = c("year", "GMST"), sep = ";", header = TRUE)


# load time series of interest
ts_station <- read.csv("\\..\\annual_max_pr_station_annual.csv", col.names = c("year", "pr"), sep = ";", header = TRUE)
ts_era5 <- read.csv("\\..\\annual_max_pr_era5land_annual.csv", col.names = c("year", "pr"), sep = ";", header = TRUE)
ts_cerra <- read.csv("\\..\\annual_max_pr_cerra_annual.csv", col.names = c("year", "pr"), sep = ";", header = TRUE)


#load historical model time series
model_climex_hist <- read.csv("\\..\\annual_max_pr_CLIMEXI_hist.csv", col.names = c("year", "pr", "GMST"), sep = ";", header = TRUE)

#load future model time series
model_climex_fut <- read.csv("\\..\\annual_max_pr_CLIMEXI_fut.csv", col.names = c("year", "pr", "GMST"), sep = ";", header = TRUE)


# combine into single dataframe
df_station <- merge(gmst, ts_station, by = "year", all = TRUE)
df_era5 <- merge(gmst[77:151, ], ts_era5, by = "year", all = TRUE)
df_cerra <- merge(gmst[112:147, ], ts_cerra, by = "year", all = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot the yearly time series
plot(df_station$year, df_station$pr, type = "s", lwd = 2, xlab = "year", ylab = "pr", bg = "white")
abline(h = 0, lty = 2)

# add a point for the event of interest
points(df_station[df_station$year == 2024,c("year", "pr")],pch = 22, bg = "magenta", lwd = 2, cex = 1.2)

# It's useful to add a simple Loess smoother to see if there is any trend over time
lines(df_station$year, fitted(loess(pr ~ year, df_station)), col = "forestgreen", lty = "32", lwd = 3)

legend("bottomleft", c("pr", "Loess smoothed"), lty = c("solid","32"), lwd = 2, col = c("black", "forestgreen"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Fitting and choosing a model                                                      ####

# We use the `fit_ns` function to fit a nonstationary distribution to the data

# We need to specify the parametric form of the distribution; 
# whether the distribution shifts linearly with GMST or scales exponentially with fixed dispersion; 
# specify the data.frame containing the data,
# the name of the columns containing the predictand and covariates,
# and whether we're interested in the lower tail of extremes.

# All methods are fully documented so that you can check the optional arguments
?fit_ns 


mdl_stn <- fit_ns("gev", "fixeddisp", df_station, "pr", c("GMST"), lower = F, ev_year =2024, ev = 144.6) #c("GMST") 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Fitted trend vs covariates                                                        ####

# We can also plot the fitted trend against each covariate to check the fit visually

# Use the `ci_cov` argument to add confidence intervals at x-values of interest
# (these are bootstrapped so will take a few seconds to run)

plot_covtrend(mdl_stn, "GMST", rp = c(20,100), add_loess = F, legend = "topleft", ylab = "pr (mm/day)", xlab = "GMST (°C)", ci_cov = data.frame(GMST = c(0.1,1.3)), ev_x =2024)
abline(h = 0, lty = 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Return level plots                                                                ####

# Return level plots are used to check whether the chosen distribution 
# is a good choice for the data.

# First, we define a 'factual' and 'counterfactual' climate to compare.
# and visualise the difference between the observations if they'd
# occurred at those values of the covariates. These must be data.frames
# with column names corresponding to the covariates in the specified model.

# our factual climate is determined by the values of the covariates in 2024
cov_2024 <- gmst[gmst$year == 2024,c("GMST"),drop = F]



# we can compare several counterfactual climates:
#     - 'hist': Preindustrial climate - 1.2C cooler than today
#     - '2deg': 2C climate above preindustrial
#     - '3deg': 3C climate above preindustrial
cov_1900 <- gmst[gmst$year == 1900,c("GMST"),drop = F]
cov_2deg <- cov_2024 + 0.6717234
cov_3deg <- cov_2024 + 1.6717234


plot_returnlevels_fix(mdl_stn, cov_f = cov_2024, cov_1900, nsamp = 1000, ylab = "pr (mm/day)",
                  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Estimate return period (pre-industrial and current climate)                        ####

rp_stn_1900 <-  return_period(mdl_stn, x = 144.6, fixed_cov = cov_1900)
rp_stn_2024 <-  return_period(mdl_stn, x = 144.6, fixed_cov = cov_2024)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Fit reanalysis dataset                                                             ####

mdl_era5 <- fit_ns("gev", "fixeddisp", df_era5, "pr", c("GMST"), lower = F) #c("GMST")
mdl_cerra <- fit_ns("gev", "fixeddisp", df_cerra, "pr", c("GMST"), lower = F) #c("GMST")

rl_era5_2024 <- eff_return_level(mdl_era5, rp = rp_stn_2024, fixed_cov = cov_2024)
rl_cerra_2024 <- eff_return_level(mdl_cerra, rp = rp_stn_2024, fixed_cov = cov_2024)


mdl_era5 <- fit_ns("gev", "fixeddisp", df_era5, "pr", c("GMST"), lower = F, ev_year = 2024, ev = rl_era5_2024) #c("GMST")
mdl_cerra <- fit_ns("gev", "fixeddisp", df_cerra, "pr", c("GMST"), lower = F, ev_year = 2024, ev = rl_cerra_2024) #c("GMST")

plot_returnlevels_fix(mdl_cerra, cov_f = cov_2024, cov_cf = cov_1900, nsamp = 1000, ylab = "pr (mm/day)", ylim= c(0,140)
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Bootstrap results (reanalysis)                                         ####

res_era5 <- boot_ci(mdl_era5, cov_f = cov_2024, cov_cf = cov_1900, ev = rl_era5_2024, nsamp = 1000)
res_cerra <- boot_ci(mdl_cerra, cov_f = cov_2024, cov_cf = cov_1900, ev = rl_cerra_2024, nsamp = 1000)

write.csv(res_era5, "\\..\\res_ERA5-Land.csv")
write.csv(res_cerra, "\\..\\res_CERRA.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > Bootstrap results (station)                                                       ####

# Bootstrapping is used to estimate a 95% confidence interval for the quantities
# of interest: the model parameters, return period, and probability ratio (PR) 
# and absolute and relative changes in intensity (dI\_abs, dI\_rel)
# associated with each of the counterfactual climates defined.

res_stn <- boot_ci(mdl_stn, cov_f = cov_2024, cov_cf = cov_1900, ev = 144.6, nsamp = 1000)
res_stn

# save the results for synthesis with the climate models later
write.csv(res_stn, "\\..\\res_Station.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > compile the results from all observations into a single file                      ####

# after bootstrapping the results for all available observational datasets, 
# compile the results into a single .csv to use in the synthesis

# Compile results
res <- t(sapply(
  list.files(path = "\\..\\/", pattern = "\\.csv$", full.names = TRUE),
  read.csv,
  row.names = "X"
))

# Use filename without path and without ".csv"
rownames(res) <- sapply(
  list.files(path = "\\..\\",pattern = "\\.csv$", full.names = TRUE),
  function(fnm) sub("res_", "", sub("\\.csv$", "", basename(fnm)))
)


# Save to CSV
write.csv(res, file = "\\..\\res_obs.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit the same nonstationary model for climate models                                 ####

# We can only attribute any trends we've identified if they also occur
# in climate models. We therefore calculate the time series of 
# the variable of interest and all relevant covariates for each climate model,
# and repeat the analysis.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# > CLIMEX Historical and Future                                                     ####
mdl_climex_hist <- fit_ns("gev", "fixeddisp", model_climex_hist, "pr", c("GMST"), lower = F) #c("GMST")

rl_climex_2024 <- eff_return_level(mdl_climex_hist, rp = rp_stn_2024, fixed_cov = cov_2024)
mdl_climex_hist <- fit_ns("gev", "fixeddisp", model_climex_hist, "pr", c("GMST"), lower = F, ev_year = 2024, ev = rl_climex_2024) #c("GMST")


plot_returnlevels_fix(mdl_climex_hist, cov_f = cov_2024, cov_cf = cov_1900, nsamp = 1000, ylab = "pr (mm/day)",
)

res_climex <- boot_ci(mdl_climex_hist, cov_f = cov_2024, cov_cf = cov_1900, ev = rl_climex_2024, nsamp = 1000)

write.csv(res_climex, "\\..\\res_CLIMEX.csv")

mdl_climex_fut <- refit(mdl_climex_hist, model_climex_fut)

#Can't use cmodel_results because we have duplicate years (50 realizations)
res_climex_2deg <- boot_ci(mdl_climex_fut, cov_f = cov_2024, cov_cf = cov_2deg, ev = rl_climex_2024, nsamp = 1000)
res_climex_3deg <- boot_ci(mdl_climex_fut, cov_f = cov_2024, cov_cf = cov_3deg, ev = rl_climex_2024, nsamp = 1000)

res_climex_2deg <- process_ci_df(res_climex_2deg)
res_climex_3deg <- process_ci_df(res_climex_3deg)

write.csv(res_climex_2deg, "\\..\\res_CLIMEX_2deg.csv")
write.csv(res_climex_3deg, "\\..\\res_CLIMEX_3deg.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We can now loop over the climate models and repeat the statistical analysis

# list all the available climate model time series
for(fnm in list.files(path = "\\..\\", pattern = "\\.csv$", full.names = T)) {
  print(fnm)
  
  # extract filenames etc
  model_name <- sub(".*pr_(.*)\\.csv", "\\1", fnm)
  print(model_name)

  # load data
  mdl_cordex <- read.csv(fnm, col.names = c("year", "GMST", "pr"), sep = ";", header = TRUE)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SPECIFY MODEL TO BE ANALYSED
  
  # value of covariate in present/factual climate
  # our factual climate is determined by the values of the covariates in 2024
  cov_2024 <- gmst[gmst$year == 2024,c("GMST"),drop = F]
  
  # we can compare several counterfactual climates:
  #     - 'hist': Preindustrial climate - 1.2C cooler than today
  #     - '2000': Year 2000 climate - 0.5C cooler than today
  #     - 'neut': 2023 climate with neutral ENSO phase
  cov_1900 <- gmst[gmst$year == 1900,c("GMST"),drop = F]
  cov_2deg <- cov_2024 + 0.6717234
  cov_3deg <- cov_2024 + 1.6717234
  cov_4deg <- cov_2024 + 2.6717234
  
  # fit the model - this is used as a template for the attribution
  mdl <- fit_ns("gev", "fixeddisp", mdl_cordex, "pr", c("GMST"), lower = F)
  
  #print(mdl$par)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #minimum year of cordex model (may vary slightly)
  min_year = min(mdl_cordex$year)
  
  # compute the climate model results
  res_df <- cmodel_results(mdl, rp = rp_stn_2024, cov_f = cov_2024, cov_hist = cov_1900, cov_fut = cov_2deg, nsamp = 1000, y_start = min_year, y_now = 2024, y_fut = 2050)
  
  #extract parameters and alpha values
  alpha_gmst <- mdl$par
  
  #calculate magnitude 
  rl_mdl_2024 <- eff_return_level(mdl, rp = rp_stn_2024, fixed_cov = cov_2024)
  
  #bootstrap to obtain ci of trend
  mdl_res <- boot_ci(mdl, cov_f = cov_2024, cov_cf = cov_1900, ev = rl_mdl_2024, nsamp = 1000)
  
  #extract trend and CI
  mdl_trend <- mdl_res[3,]
  
  print(mdl_trend)
  
  # Build the new save path
  save_path_2deg <- paste0("\\..\\", model_name, ".csv")
  save_path_pars <- paste0("\\..\\", model_name, ".csv")
  save_path_trends <- paste0("\\..\\", model_name, ".csv")
  
  write.csv(res_df, file = save_path_2deg)
  write.csv(alpha_gmst, file = save_path_pars)
  write.csv(mdl_trend, file = save_path_trends)
  
}

# Once the model fitting has been repeated for all climate models, 
# we can compile the results into a single .csv for easier handling

# Compile results
res <- t(sapply(
  list.files(path = "\\..\\", pattern = "\\.csv$", full.names = TRUE),
  read.csv,
  row.names = "X"
))

# Use filename without path and without ".csv"
rownames(res) <- sapply(
  list.files(path = "\\..\\", pattern = "\\.csv$", full.names = TRUE),
  function(fnm) sub("\\.csv$", "", basename(fnm))
)


# Save to CSV

write.csv(res, file = "\\..\\res_cordex_2deg.csv")
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Synthesis

# Having evaluated the changes in likelihood and frequency of 
# the specified event in both the observations and the climate models, 
# we now combine the two sets of results to produce a single 
# overarching attribution statement.

 # This script uses output from the functions above, but the function only 
# requires that you provide, for the models and (optionally) observations,
# a data.frame with row.names indicating the model names, and columns 
# containing the best estimate and lower and upper confidence bounds 
# of the quantity to be synthesised (either a PR or a change in intensity).

# The synthesis algorithm produces a precision-weighted average of the models and 
# (if using) the observations, then provides a precision weighted average (purple bar)
# and unweighted average (white bar) of the two sets of data. 
# For a more detailed explanation of the calculation, see Otto et al., 
# 'Formally combining different lines of evidence in extreme event attribution';
# Advances in Statistical Climatology, Meteorology and Oceanography (2024).

# load the obs 
obs_res <- read.csv("\\..\\res_obs.csv", row.names = "X")
model_res <- read.csv("\\..\\res_cordex_2deg.csv", row.names = "X")
model_ensemble_res <- read.csv("\\..\\res_CLIMEX_CORDEX_ensembles.csv", row.names = "X")

#load the cordex
model_res_3deg <- read.csv("\\..\\res_cordex_3deg.csv", row.names = "X")

#load the ensembles (climex & cordex)
model_ensemble_2deg <- read.csv("\\..\\res_CLIMEX_CORDEX_ensembles_2deg.csv", row.names = "X")
model_ensemble_3deg <- read.csv("\\..\\res_CLIMEX_CORDEX_ensembles_3deg.csv", row.names = "X")

####################################################
# run the synthesis

#Current climate
synth_pr_attr <- synthesis(obs_in = obs_res[, grepl("PR", colnames(obs_res))], model_res[, grepl("attr_PR", colnames(model_res))], synth_type = "PR")
synth_dI_attr <- synthesis(obs_in = obs_res[, grepl("dI_rel", colnames(obs_res))], model_res[, grepl("attr_dI.rel", colnames(model_res))], synth_type = "rel")

#Current climate
synth_pr_attr_ensemble <- ensemble_synthesis(obs_in = obs_res[, grepl("PR", colnames(obs_res))], model_ensemble_res[, grepl("attr_PR", colnames(model_ensemble_res))], synth_type = "PR")
synth_dI_attr_ensemble <- ensemble_synthesis(obs_in = obs_res[, grepl("dI_rel", colnames(obs_res))], model_ensemble_res[, grepl("attr_dI.rel", colnames(model_ensemble_res))], synth_type = "rel")

#Future warming 2 deg
synth_pr_proj_2deg <- synthesis(obs_in = NA, model_res[, grepl("proj_PR", colnames(model_res))], synth_type = "PR")
synth_dI_proj_2deg  <- synthesis(obs_in = NA, model_res[, grepl("proj_dI.rel", colnames(model_res))], synth_type = "rel")

#Future warming 3 deg
synth_pr_proj_3deg <- synthesis(obs_in = NA, model_res_3deg[, grepl("proj_PR", colnames(model_res))], synth_type = "PR")
synth_dI_proj_3deg <- synthesis(obs_in = NA, model_res_3deg[, grepl("proj_dI.rel", colnames(model_res))], synth_type = "rel")

#Future warming 2 deg (ensemble)
synth_pr_proj_2deg_ensemble <- ensemble_synthesis(obs_in = NA, model_ensemble_2deg[, grepl("proj_PR", colnames(model_ensemble_2deg))], synth_type = "PR")
synth_dI_proj_2deg_ensemble  <- ensemble_synthesis(obs_in = NA, model_ensemble_2deg[, grepl("proj_dI.rel", colnames(model_ensemble_2deg))], synth_type = "rel")

#Future warming 3 deg (ensemble)
synth_pr_proj_3deg_ensemble <- ensemble_synthesis(obs_in = NA, model_ensemble_3deg[, grepl("proj_PR", colnames(model_ensemble_3deg))], synth_type = "PR")
synth_dI_proj_3deg_ensemble <- ensemble_synthesis(obs_in = NA,model_ensemble_3deg[, grepl("proj_dI.rel", colnames(model_ensemble_3deg))], synth_type = "rel")

####################################################
#plotting

#Prepare plotting (pre-industrial to current) (PR)
synth_pr_attr$df[39,] <- synth_pr_attr$df[37,]
synth_pr_attr$df[36,] <- synth_pr_attr_ensemble$df[6,]
synth_pr_attr$df[37,] <- synth_pr_attr_ensemble$df[5,]
synth_pr_attr$df[38,] <- synth_pr_attr_ensemble$df[7,]
synth_pr_attr$df[40,] <- synth_pr_attr_ensemble$df[8,]

synth_pr_attr$df[36,]$group <- "ens"
synth_pr_attr$df[37,]$group <- "ens"
synth_pr_attr$df[38,]$group <- "ens_average"
synth_pr_attr$df[40,]$group <- "synth_ens"
synth_pr_attr$df[39,]$group <- "synth_model"

synth_pr_attr$df[37,]$model <- "ClimEx Ensemble"
synth_pr_attr$df[38,]$model <- "Ensembles"

synth_pr_attr$df[39,]$l_wb <- NA
synth_pr_attr$df[39,]$u_wb <- NA
synth_pr_attr$df[39,]$est <- NA

#Plot
par(mar = c(5, 27, 4, 2)) 
plot_synthesis_mod(synth_pr_attr)

####################################################
#Prepare plotting (pre-industrial to current) (dI)

synth_dI_attr$df[39,] <- synth_dI_attr$df[37,]
synth_dI_attr$df[36,] <- synth_dI_attr_ensemble$df[6,]
synth_dI_attr$df[37,] <- synth_dI_attr_ensemble$df[5,]
synth_dI_attr$df[38,] <- synth_dI_attr_ensemble$df[7,]
synth_dI_attr$df[40,] <- synth_dI_attr_ensemble$df[8,]

synth_dI_attr$df[36,]$group <- "ens"
synth_dI_attr$df[37,]$group <- "ens"
synth_dI_attr$df[38,]$group <- "ens_average"
synth_dI_attr$df[40,]$group <- "synth_ens"
synth_dI_attr$df[39,]$group <- "synth_model"

synth_dI_attr$df[37,]$model <- "ClimEx Ensemble"
synth_dI_attr$df[38,]$model <- "Ensembles"

synth_dI_attr$df[39,]$l_wb <- NA
synth_dI_attr$df[39,]$u_wb <- NA
synth_dI_attr$df[39,]$est <- NA

#Plot
par(mar = c(5, 27, 4, 2))
plot_synthesis_mod(synth_dI_attr)

####################################################
#Prepare plotting (pre-industrial to current) (PR) #2deg or 3deg

synth_pr_proj_2deg$df[33,] <- synth_pr_proj_2deg_ensemble$df[1,]
synth_pr_proj_2deg$df[34,] <- synth_pr_proj_2deg_ensemble$df[3,]

synth_pr_proj_2deg$df[32,]$group <- "ens"
synth_pr_proj_2deg$df[33,]$group <- "ens"
synth_pr_proj_2deg$df[34,]$group <- "synth"

synth_pr_proj_2deg$df[32,]$model <- "CORDEX Ensemble"
synth_pr_proj_2deg$df[33,]$model <- "ClimEx Ensemble"
synth_pr_proj_2deg$df[34,]$model <- "Ensembles"

#Plot
par(mar = c(5, 27, 4, 2)) 
plot_synthesis_fut(synth_pr_proj_2deg)

####################################################
#Prepare plotting (current to to future) (dI) #2deg or 3deg

synth_dI_proj_2deg$df[33,] <- synth_dI_proj_2deg_ensemble$df[1,]
synth_dI_proj_2deg$df[34,] <- synth_dI_proj_2deg_ensemble$df[3,]

synth_dI_proj_2deg$df[32,]$group <- "ens"
synth_dI_proj_2deg$df[33,]$group <- "ens"
synth_dI_proj_2deg$df[34,]$group <- "synth"

synth_dI_proj_2deg$df[32,]$model <- "CORDEX Ensemble"
synth_dI_proj_2deg$df[33,]$model <- "ClimEx Ensemble"
synth_dI_proj_2deg$df[34,]$model <- "Ensembles"

#Plot
par(mar = c(5, 27, 4, 2))
plot_synthesis_fut(synth_dI_proj_2deg, hide_labels = TRUE)




