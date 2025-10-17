# install.packages(c("metafor","mixmeta","dplyr","ggplot2","splines"))
library(metafor)
library(mixmeta)
library(dplyr)
library(ggplot2)
library(splines)

# 1) Assumed input data frame
dat <- read.csv("meta_input.csv")

# Required columns (example names):
# study_id      - unique study identifier
# esm           - effect measure reported (HR, RR or OR) numeric
# lower_ci      - lower bound of 95% CI for esm
# upper_ci      - upper bound of 95% CI for esm
# esm_type      - "HR", "RR", or "OR" (character)  (optional but useful for sensitivity)
# n_events      - number of dementia events (optional)
# n_total       - total sample (optional)
# follow_up_yrs - follow-up duration in years (numeric)
# risk_bias     - "low","moderate","high" or numeric risk-of-bias score (optional)
# moderatorX    - study-level moderator(s): e.g. "social_engagement_mean", "physical_activity", etc.

# 2)Prepare effect estimates
# Treat RR and OR as HRs (per protocol) as outcome number may be low.
dat2 <- dat %>%
  mutate(
    esm = as.numeric(esm),
    lower_ci = as.numeric(lower_ci),
    upper_ci = as.numeric(upper_ci),
    yi = log(esm),
    sei = (log(upper_ci) - log(lower_ci)) / (2 * 1.96), 
    vi = sei^2
  )

# 3) Primary random-effects meta-analysis using REML and HKSJ
res_reml_kn <- rma(yi = yi, vi = vi, data = dat2,
                   method = "REML",
                   test = "knha",
                   slab = study_id)
print(res_reml_kn)

pooled_log <- coef(res_reml_kn)
pooled_se  <- res_reml_kn$se
pooled_ci <- confint(res_reml_kn, level = 0.95, digits=4)  
pooled_HR <- exp(pooled_log)
pooled_CI_HR <- exp(pooled_ci$random)  
cat("Pooled HR (REML, HKSJ): ", round(pooled_HR,3), " (",
    round(pooled_CI_HR[1],3), ", ", round(pooled_CI_HR[2],3), ")\n", sep="")

# 4) Heterogeneity: tau^2, I^2, and 95% prediction interval
tau2_hat <- res_reml_kn$tau2
tau_hat  <- sqrt(tau2_hat)

# Compute I^2 (Higgins & Thompson) — using common formula;I^2 = tau^2 / (tau^2 + mean(vi))
mean_vi <- mean(dat2$vi, na.rm=TRUE)
I2 <- 100 * tau2_hat / (tau2_hat + mean_vi)
cat("tau^2 =", signif(tau2_hat,3), "; tau =", signif(tau_hat,3), "; I^2 =", signif(I2,3), "%\n")

# Prediction interval (95%): metafor predict provides pi.lb/pi.ub typically
pred <- predict(res_reml_kn, transf = exp)  
print(pred)  

# If predict didn't produce pi.lb/pi.ub (older versions), compute manually; On log scale: pi = theta_hat +- t_{k-2} * sqrt(se^2 + tau^2)
k <- res_reml_kn$k
t_crit <- qt(0.975, df = k - 2)
pred_log <- coef(res_reml_kn)
pred_se <- res_reml_kn$se
pi_log_lb <- pred_log - t_crit * sqrt(pred_se^2 + tau2_hat)
pi_log_ub <- pred_log + t_crit * sqrt(pred_se^2 + tau2_hat)
cat("Prediction interval (HR): ", round(exp(pi_log_lb),3), " to ", round(exp(pi_log_ub),3), "\n")

# 5) Meta-regression with study-level moderators (primary research question: model lifestyle indicators)
# Example model:e.g., dat2$social_engagement_mean
if ("social_engagement_mean" %in% names(dat2)) {
  res_mod1 <- rma(yi = yi, vi = vi, mods = ~ social_engagement_mean, data = dat2,
                  method = "REML", test = "knha")
  print(res_mod1)
  cat("Meta-regression: coefficient for social_engagement_mean (log HR):",
      signif(coef(res_mod1)[2],3), "\n")
} else {
  message("Moderator 'social_engagement_mean' not found in data — skip meta-regression example.")
}

# For multiple moderators:
# res_mod_multi <- rma(yi=yi, vi=vi, mods = ~ social_engagement_mean + mean_age + percent_male, data=dat2,
#                      method="REML", test="knha")
# print(res_mod_multi)

# 6) Dose-response meta-analysis (continuous lifestyle exposure)
#    Example using mixmeta for linear and restricted cubic spline
# ---------------------------
# This requires each study to report exposure-specific log HR and variances across multiple exposure levels
# If you have one effect per study with a continuous study-level mean only, you can model that with rma (meta-regression)
#
# If you have multivariate estimates per study (e.g., several exposure categories per study), transform into
# a long format with columns: study, x (dose), yi (logHR), vi (variance), then use mixmeta.

if (all(c("dose_x","yi","vi","study_id") %in% names(dat2))) {
  # Example: linear dose-response model
  # Long-format: dat_dr <- dat2_long with multiple rows per study
  dat_dr <- dat2 %>% rename(study = study_id, x = dose_x)
  mix_lin <- mixmeta(y = yi, S = vi, X = model.matrix(~ x, data = dat_dr), id = study, data = dat_dr, method = "reml")
  summary(mix_lin)
  
  # Spline (3 knots) example: use natural cubic splines
  knots <- quantile(dat_dr$x, probs = c(0.10, 0.5, 0.90), na.rm=TRUE)
  Xs <- model.matrix(~ ns(x, knots = knots), data = dat_dr)
  mix_spline <- mixmeta(y = yi, S = vi, X = Xs, id = study, data = dat_dr, method = "reml")
  summary(mix_spline)
} else {
  message("Dose-response analysis skipped: no 'dose_x' or multiple-exposure rows found. If you have per-study dose-category estimates, restructure to long format and re-run.")
}

# 7) Leave-one-out sensitivity analysis
loo <- leave1out(res_reml_kn)
print(loo)   # shows pooled estimates with each study removed

# Quick plot of leave-one-out pooled HRs
loo_df <- data.frame(study = rownames(loo), 
                     yi_leave = loo$estimate, 
                     ci.lb = loo$ci.lb, 
                     ci.ub = loo$ci.ub) %>%
  mutate(HR = exp(yi_leave), HR_lb = exp(ci.lb), HR_ub = exp(ci.ub))

ggplot(loo_df, aes(x = reorder(study, HR), y = HR)) +
  geom_point() + geom_errorbar(aes(ymin = HR_lb, ymax = HR_ub), width = 0.2) +
  coord_flip() + labs(y = "Pooled HR (leave-one-out)", x = "Study removed")

# 8) Additional sensitivity analyses
#    (a) exclude OR/RR studies
#    (b) exclude follow-up < 5 yrs
#    (c) exclude high risk of bias

# (a) Exclude studies reporting OR/RR
if ("esm_type" %in% names(dat2)) {
  dat_hr_only <- dat2 %>% filter(toupper(esm_type) == "HR")
  res_hronly <- rma(yi=yi, vi=vi, data=dat_hr_only, method="REML", test="knha")
  cat("Results excluding non-HR studies:\n"); print(res_hronly)
} else {
  message("Column 'esm_type' not found; cannot perform 'exclude OR/RR' sensitivity analysis automatically.")
}

# (b) Exclude follow-up < 5 years
if ("follow_up_yrs" %in% names(dat2)) {
  dat_longfu <- dat2 %>% filter(is.na(follow_up_yrs) | follow_up_yrs >= 5)
  res_longfu <- rma(yi=yi, vi=vi, data=dat_longfu, method="REML", test="knha")
  cat("Results excluding follow-up < 5 years:\n"); print(res_longfu)
} else {
  message("Column 'follow_up_yrs' not found; cannot perform follow-up sensitivity analysis automatically.")
}

# (c) Exclude high risk of bias
if ("risk_bias" %in% names(dat2)) {
  dat_lowrob <- dat2 %>% filter(risk_bias != "high")
  res_lowrob <- rma(yi=yi, vi=vi, data=dat_lowrob, method="REML", test="knha")
  cat("Results excluding high risk-of-bias studies:\n"); print(res_lowrob)
} else {
  message("Column 'risk_bias' not found; cannot perform risk-of-bias sensitivity analysis automatically.")
}

# 9) Publication bias: funnel plot, Egger's test, trim-and-fill
# Funnel plot (log scale)
funnel(res_reml_kn, xlab = "Log(HR)", main = "Funnel plot (log scale)")

# Egger's test (regression test for small-study effects)
# Use metafor::regtest; for random-effects model, use model="rma"
egger <- regtest(res_reml_kn, model = "rma")
print(egger)   # p-value for asymmetry

# Trim-and-fill to assess impact of missing studies
tf <- trimfill(res_reml_kn)
print(tf)
# Plot trim-and-fill results on funnel
funnel(tf, xlab = "Log(HR)", main = "Trim-and-fill funnel plot (log scale)")
cat("Trim-and-fill pooled HR (with imputed studies):", exp(coef(tf)), "\n")

# 10) Save outputs / create result tables
# Create a results summary table
results_summary <- data.frame(
  pooled_logHR = coef(res_reml_kn),
  pooled_logHR_se = res_reml_kn$se,
  pooled_HR = exp(coef(res_reml_kn)),
  pooled_HR_lb = exp(res_reml_kn$ci.lb),
  pooled_HR_ub = exp(res_reml_kn$ci.ub),
  tau2 = tau2_hat,
  tau = tau_hat,
  I2 = I2,
  k = res_reml_kn$k
)
print(results_summary)

# Save results if desired
# write.csv(results_summary, "meta_pooled_summary.csv", row.names = FALSE)
