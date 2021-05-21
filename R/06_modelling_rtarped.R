#####################################
##### MODELLING lreg_tarpedarea #####
#####################################

### IMPORTANT NOTE: to avoid creating notes for unquoted variables, I must add the following code at
# the beginning of every source file (e.g. R/myscript.R) that uses unquoted variables, so in front of
# all my scripts doing any kind of analyses (otherwise, I should always assign each variable, e.g.
# mydata$myvariable, which is quite time consuming and wearisome).
# if(getRversion() >= "2.29.1")  utils::globalVariables(c(
#   "manager_id", "xp_id", "lreg_tarpedarea", "eff_eradication",
#   "latitude", "elevation", "slope", "coarse_env", "obstacles", "flood",
#   "geomem", "maxveg", "uprootexcav", "stand_surface", "fully_tarped", "distance", "tarping_duration",
#   "stripsoverlap_ok", "tarpfix_pierced", "tarpfix_multimethod", "sedicover_height",
#   "trench_depth", "plantation", "add_control",
#   "pb_fixation", "pb_durability")) # WRONG LIST! Should be updated if I keep this chunk of code!
# -------------------------------------------------------------------------------------------------------- #





# --------------------------------------------------- #
##### Data preparation for modelling 'lreg_tarpedarea' #####
# --------------------------------------------------- #

# List of used packages (for publication or package building): here, readr, MuMIn, glmmTMB, ggplot2,
# (broom.mixed), stats, DHARMa, performance

.pardefault <- par() # To save the default graphical parameters (in case I want to restore them).

library(jk.dusz.tarping)
readr::read_csv(here::here("mydata", "rtarped.csv"), col_names = TRUE, col_types =
                  readr::cols(
                    manager_id = readr::col_factor(),
                    xp_id = readr::col_factor(),
                    lreg_tarpedarea = readr::col_factor(),
                    geomem = readr::col_factor(c("0", "1")),
                    woven_geotex = readr::col_factor(c("0", "1")),
                    maxveg = readr::col_factor(c("0", "1")),
                    uprootexcav = readr::col_factor(c("0", "1")),
                    levelling = readr::col_factor(c("0", "1")),
                    fully_tarped = readr::col_factor(c("0", "1")),
                    stripsoverlap_ok = readr::col_factor(c("0", "1")),
                    tarpfix_multimethod = readr::col_factor(c("0", "1")),
                    tarpfix_pierced = readr::col_factor(c("0", "1")),
                    plantation = readr::col_factor(c("0", "1")),
                    repairs = readr::col_factor(c("0", "1")),
                    add_control = readr::col_factor(c("0", "1")),
                    pb_fixation = readr::col_factor(c("0", "1")),
                    pb_durability = readr::col_factor(c("0", "1")),
                    reg_edges = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(latitude = jitter(x = latitude, factor = 0.1)) %>%
  dplyr::mutate(longitude = jitter(x = longitude, factor = 0.1)) -> rtarped # Added a very small amount of
# noise to coordinates to avoid groups with exactly similar coordinates (related to low Lat/Long resolution)
# which prevent the proper use of the DHARMa package autocorrelation test!
summary(rtarped)





# ---------------------------------------- #
##### Building of the candidate models #####
# ---------------------------------------- #

Cand.mod <- list()
R.ajust <- data.frame(Model=integer(0), R2=numeric(0)) # Creates an empty data.frame with 2 variables



### Testing the relevance of the random effect structure:
m0.glm <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~1, data = rtarped,
                           family = stats::binomial(link = "logit"), REML = FALSE)
m0.glmer <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~(1|manager_id), data = rtarped,
                             family = stats::binomial(link = "logit"), REML = FALSE)
m0.glmer1 <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~(1|manager_id) + (1|xp_id), data = rtarped,
                             family = stats::binomial(link = "logit"), REML = FALSE)
aic.glm <- AIC(logLik(m0.glm))
aic.glmer <- AIC(logLik(m0.glmer))
aic.glmer1 <- AIC(logLik(m0.glmer1))

# Likelihood Ratio Test:
null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer1)
pchisq(as.numeric(null.id), df=1, lower.tail=F) # The Likelihood Ratio Test is NOT significant suggesting
# that the use of the random effect structure is not necessary! HOWEVER, model diagnostics for subsequent
# models have shown that failing to include "manager_id" as a random effect leads to model misspecification.
# Consequently, and as initially planned, we included a random structure within our candiate models.
rm(m0.glm, m0.glmer, m0.glmer1, aic.glm, aic.glmer, aic.glmer1, null.id)



###########################
###########################
###########################
###########################
###########################

### Model diagnostics:
# One plot of residuals (https://www.r-bloggers.com/2011/07/model-validation-interpreting-residual-plots/):
plot(fitted(Cand.mod[[1]]), residuals(Cand.mod[[1]]),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Cand.mod[[1]]), residuals(Cand.mod[[1]])))

# If mixed model:
# Check for residual pattern within groups (levels of random factor) and difference between groups
xyplot(residuals(glmm1) ~ fitted(glmm1) | Count$plot, main = "glmm1 – full model by plot",
       panel=function(x, y){
         panel.xyplot(x, y)
         panel.loess(x, y, span = 0.75)
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       })


### EN cas de models mixtes, relire la vignette DHARMa, car ça change!

###########################
###########################
###########################
###########################
###########################



##### Model 1 (null model) #####
# ------------------------------

Cand.mod[[1]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~1 + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[1]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok (if
# # random effects are included, otherwise Moran's I test is significant)!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[1]])
# performance::check_autocorrelation(Cand.mod[[1]])
# performance::check_collinearity(Cand.mod[[1]])
# performance::check_singularity(Cand.mod[[1]])
#
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[1]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[1]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[1]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[1]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[1]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[1]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=1, R2=R2[[1]]))



##### Model 2 #####
# -----------------

Cand.mod[[2]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[2]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[2]])
# performance::check_autocorrelation(Cand.mod[[2]])
# performance::check_collinearity(Cand.mod[[2]])
# performance::check_singularity(Cand.mod[[2]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[2]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[2]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[2]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[2]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[2]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[2]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=2, R2=R2[[1]]))



##### Model 3 #####
# -----------------

Cand.mod[[3]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + geomem
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[3]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[3]])
# performance::check_autocorrelation(Cand.mod[[3]])
# performance::check_collinearity(Cand.mod[[3]])
# performance::check_singularity(Cand.mod[[3]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[3]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[3]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[3]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[3]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[3]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[3]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=3, R2=R2[[1]]))



##### Model 4 #####
# -----------------

Cand.mod[[4]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + levelling
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[4]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[4]])
# performance::check_autocorrelation(Cand.mod[[4]])
# performance::check_collinearity(Cand.mod[[4]])
# performance::check_singularity(Cand.mod[[4]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[4]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[4]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[4]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[4]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[4]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[4]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=4, R2=R2[[1]]))



##### Model 5 #####
# -----------------

Cand.mod[[5]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + log2(stand_surface)
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[5]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[5]])
# performance::check_autocorrelation(Cand.mod[[5]])
# performance::check_collinearity(Cand.mod[[5]])
# performance::check_singularity(Cand.mod[[5]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[5]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[5]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[5]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[5]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[5]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[5]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=5, R2=R2[[1]]))



##### Model 6 #####
# -----------------

Cand.mod[[6]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + stripsoverlap_ok
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[6]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[6]])
# performance::check_autocorrelation(Cand.mod[[6]])
# performance::check_collinearity(Cand.mod[[6]])
# performance::check_singularity(Cand.mod[[6]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[6]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[6]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[6]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[6]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[6]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[6]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=6, R2=R2[[1]]))




##### Model 7 #####
# -----------------

Cand.mod[[7]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + pb_fixation
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[7]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[7]])
# performance::check_autocorrelation(Cand.mod[[7]]) # Close!
# performance::check_collinearity(Cand.mod[[7]])
# performance::check_singularity(Cand.mod[[7]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[7]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[7]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[7]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[7]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[7]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[7]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=7, R2=R2[[1]]))



##### Model 8 #####
# -----------------

Cand.mod[[8]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + tarpfix_pierced
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[8]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[8]])
# performance::check_autocorrelation(Cand.mod[[8]])
# performance::check_collinearity(Cand.mod[[8]])
# performance::check_singularity(Cand.mod[[8]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[8]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[8]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[8]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[8]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[8]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[8]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=8, R2=R2[[1]]))



##### Model 9 #####
# -----------------

Cand.mod[[9]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + woven_geotex
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[9]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[9]])
# performance::check_autocorrelation(Cand.mod[[9]])
# performance::check_collinearity(Cand.mod[[9]])
# performance::check_singularity(Cand.mod[[9]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[9]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[9]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[9]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[9]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[9]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[9]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=9, R2=R2[[1]]))



##### Model 10 #####
# -----------------

Cand.mod[[10]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + plantation
                                  + (1|xp_id), data = rtarped,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[10]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[10]])
# performance::check_autocorrelation(Cand.mod[[10]]) # Nope! The cause of this is not clear. However,
# # adding another random effect only makes things worse and this problem does not seem to affect inferences.
# # Therefore, we keep the model as is to avoid overfitting.
# performance::check_collinearity(Cand.mod[[10]])
# performance::check_singularity(Cand.mod[[10]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[10]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[10]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[10]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[10]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[10]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[10]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=10, R2=R2[[1]]))



##### Model 11 #####
# -----------------

Cand.mod[[11]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarpfix_pierced
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[11]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[11]])
# performance::check_autocorrelation(Cand.mod[[11]]) # Not ok! The cause of this is not clear. However,
# # adding another random effect only makes things worse and this problem does not seem to affect inferences.
# # Therefore, we keep the model as is to avoid overfitting.
# performance::check_collinearity(Cand.mod[[11]])
# performance::check_singularity(Cand.mod[[11]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[11]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[11]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[11]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[11]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[11]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[11]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=11, R2=R2[[1]]))



##### Model 12 #####
# -----------------

Cand.mod[[12]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[12]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Not ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[12]])
# performance::check_autocorrelation(Cand.mod[[12]]) # Not ok! The cause of this is not clear. However,
# # adding another random effect only makes things worse and this problem does not seem to affect inferences.
# # Therefore, we keep the model as is to avoid overfitting.
# performance::check_collinearity(Cand.mod[[12]])
# performance::check_singularity(Cand.mod[[12]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[12]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[12]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[12]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[12]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[12]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[12]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=12, R2=R2[[1]]))



##### Model 13 #####
# -----------------

Cand.mod[[13]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + woven_geotex
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[13]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[13]])
# performance::check_autocorrelation(Cand.mod[[13]]) # Same as before!
# performance::check_collinearity(Cand.mod[[13]])
# performance::check_singularity(Cand.mod[[13]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[13]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[13]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[13]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[13]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[13]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[13]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=13, R2=R2[[1]]))



##### Model 14 #####
# -----------------

Cand.mod[[14]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + levelling
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[14]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[14]])
# performance::check_autocorrelation(Cand.mod[[14]]) # Close!
# performance::check_collinearity(Cand.mod[[14]])
# performance::check_singularity(Cand.mod[[14]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[14]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[14]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[14]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[14]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[14]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[14]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=14, R2=R2[[1]]))



##### Model 15 #####
# -----------------

Cand.mod[[15]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + log2(sedicover_height + 1)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[15]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Nope! This time, including another
# # random effect factor (i.e. manager_id) solves most problems. However, as it creates singularity and
# # does not improve inference or predictive power, we chose to keep the most parsimonious model!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[15]])
# performance::check_autocorrelation(Cand.mod[[15]])
# performance::check_collinearity(Cand.mod[[15]])
# performance::check_singularity(Cand.mod[[15]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[15]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[15]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[15]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[15]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[15]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[15]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=15, R2=R2[[1]]))



##### Model 16 #####
# -----------------

Cand.mod[[16]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + pb_fixation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[16]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[16]])
# performance::check_autocorrelation(Cand.mod[[16]]) # Nope (but close)! Same as model 15.
# performance::check_collinearity(Cand.mod[[16]])
# performance::check_singularity(Cand.mod[[16]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[16]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[16]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[16]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[16]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[16]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[16]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=16, R2=R2[[1]]))



##### Model 17 #####
# -----------------

Cand.mod[[17]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + tarpfix_pierced
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[17]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[17]])
# performance::check_autocorrelation(Cand.mod[[17]]) # Close!
# performance::check_collinearity(Cand.mod[[17]])
# performance::check_singularity(Cand.mod[[17]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[17]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[17]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[17]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[17]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[17]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[17]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=17, R2=R2[[1]]))



##### Model 18 #####
# -----------------

Cand.mod[[18]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[18]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope! Same as before.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[18]])
# performance::check_autocorrelation(Cand.mod[[18]]) # Close!
# performance::check_collinearity(Cand.mod[[18]])
# performance::check_singularity(Cand.mod[[18]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[18]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[18]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[18]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[18]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[18]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[18]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=18, R2=R2[[1]]))



##### Model 19 #####
# -----------------

Cand.mod[[19]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~plantation + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[19]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[19]])
# performance::check_autocorrelation(Cand.mod[[19]]) # Close!
# performance::check_collinearity(Cand.mod[[19]])
# performance::check_singularity(Cand.mod[[19]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[19]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[19]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[19]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[19]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[19]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[19]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=19, R2=R2[[1]]))



##### Model 20 #####
# -----------------

Cand.mod[[20]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~plantation + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[20]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[20]])
# performance::check_autocorrelation(Cand.mod[[20]]) # Close!
# performance::check_collinearity(Cand.mod[[20]])
# performance::check_singularity(Cand.mod[[20]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[20]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[20]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[20]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[20]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[20]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[20]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=20, R2=R2[[1]]))



##### Model 21 #####
# -----------------

Cand.mod[[21]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[21]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[21]])
# performance::check_autocorrelation(Cand.mod[[21]])
# performance::check_collinearity(Cand.mod[[21]])
# performance::check_singularity(Cand.mod[[21]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[21]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[21]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[21]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[21]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[21]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[21]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=21, R2=R2[[1]]))



##### Model 22 #####
# -----------------

Cand.mod[[22]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + levelling
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[22]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[22]])
# performance::check_autocorrelation(Cand.mod[[22]])
# performance::check_collinearity(Cand.mod[[22]])
# performance::check_singularity(Cand.mod[[22]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[22]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[22]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[22]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[22]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[22]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[22]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=22, R2=R2[[1]]))



##### Model 23 #####
# -----------------

Cand.mod[[23]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[23]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before (including another
# # random effect factor improves residuals but creates a singular model without affecting coefficients or
# # standard errors, AIC, etc., so I decided to keep the model as is)
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[23]])
# performance::check_autocorrelation(Cand.mod[[23]]) # Close!
# performance::check_collinearity(Cand.mod[[23]])
# performance::check_singularity(Cand.mod[[23]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[23]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[23]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[23]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[23]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[23]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[23]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=23, R2=R2[[1]]))



##### Model 24 #####
# ------------------

Cand.mod[[24]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + stripsoverlap_ok
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[24]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[24]])
# performance::check_autocorrelation(Cand.mod[[24]])
# performance::check_collinearity(Cand.mod[[24]])
# performance::check_singularity(Cand.mod[[24]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[24]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[24]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[24]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[24]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[24]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[24]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=24, R2=R2[[1]]))



##### Model 25 #####
# -----------------

Cand.mod[[25]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + pb_fixation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[25]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[25]])
# performance::check_autocorrelation(Cand.mod[[25]]) # Close!
# performance::check_collinearity(Cand.mod[[25]])
# performance::check_singularity(Cand.mod[[25]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[25]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[25]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[25]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[25]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[25]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[25]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=25, R2=R2[[1]]))



##### Model 26 #####
# -----------------

Cand.mod[[26]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~tarpfix_pierced + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[26]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[26]])
# performance::check_autocorrelation(Cand.mod[[26]]) # Close!
# performance::check_collinearity(Cand.mod[[26]])
# performance::check_singularity(Cand.mod[[26]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[26]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[26]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[26]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[26]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[26]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[26]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=26, R2=R2[[1]]))



##### Model 27 #####
# -----------------

Cand.mod[[27]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(sedicover_height + 1) + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[27]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[27]])
# performance::check_autocorrelation(Cand.mod[[27]])
# performance::check_collinearity(Cand.mod[[27]])
# performance::check_singularity(Cand.mod[[27]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[27]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[27]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[27]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[27]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[27]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[27]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=27, R2=R2[[1]]))



##### Model 28 #####
# -----------------

Cand.mod[[28]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(sedicover_height + 1) + stripsoverlap_ok
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[28]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[28]])
# performance::check_autocorrelation(Cand.mod[[28]])
# performance::check_collinearity(Cand.mod[[28]])
# performance::check_singularity(Cand.mod[[28]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[28]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[28]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[28]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[28]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[28]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[28]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=28, R2=R2[[1]]))



##### Model 29 #####
# -----------------

Cand.mod[[29]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(sedicover_height + 1) + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[29]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope! Same as before.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[29]])
# performance::check_autocorrelation(Cand.mod[[29]])
# performance::check_collinearity(Cand.mod[[29]])
# performance::check_singularity(Cand.mod[[29]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[29]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[29]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[29]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[29]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[29]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[29]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=29, R2=R2[[1]]))



##### Model 30 #####
# -----------------

Cand.mod[[30]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~stripsoverlap_ok + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[30]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[30]])
# performance::check_autocorrelation(Cand.mod[[30]])
# performance::check_collinearity(Cand.mod[[30]])
# performance::check_singularity(Cand.mod[[30]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[30]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[30]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[30]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[30]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[30]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[30]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=30, R2=R2[[1]]))



##### Model 31 #####
# -----------------

Cand.mod[[31]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~stripsoverlap_ok + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[31]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[31]])
# performance::check_autocorrelation(Cand.mod[[31]])
# performance::check_collinearity(Cand.mod[[31]])
# performance::check_singularity(Cand.mod[[31]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[31]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[31]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[31]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[31]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[31]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[31]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=31, R2=R2[[1]]))



##### Model 32 #####
# -----------------

Cand.mod[[32]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[32]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[32]])
# performance::check_autocorrelation(Cand.mod[[32]]) # Same as before.
# performance::check_collinearity(Cand.mod[[32]])
# performance::check_singularity(Cand.mod[[32]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[32]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[32]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[32]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[32]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[32]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[32]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=32, R2=R2[[1]]))



##### Model 33 #####
# -----------------

Cand.mod[[33]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[33]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[33]])
# performance::check_autocorrelation(Cand.mod[[33]])
# performance::check_collinearity(Cand.mod[[33]])
# performance::check_singularity(Cand.mod[[33]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[33]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[33]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[33]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[33]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[33]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[33]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=33, R2=R2[[1]]))



##### Model 34 #####
# -----------------

Cand.mod[[34]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control + tarpfix_pierced
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[34]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[34]])
# performance::check_autocorrelation(Cand.mod[[34]])
# performance::check_collinearity(Cand.mod[[34]])
# performance::check_singularity(Cand.mod[[34]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[34]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[34]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[34]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[34]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[34]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[34]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=34, R2=R2[[1]]))



##### Model 35 #####
# -----------------

Cand.mod[[35]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[35]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[35]])
# performance::check_autocorrelation(Cand.mod[[35]])
# performance::check_collinearity(Cand.mod[[35]])
# performance::check_singularity(Cand.mod[[35]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[35]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[35]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[35]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[35]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[35]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[35]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=35, R2=R2[[1]]))



##### Model 36 #####
# -----------------

Cand.mod[[36]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + add_control + woven_geotex
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[36]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[36]])
# performance::check_autocorrelation(Cand.mod[[36]])
# performance::check_collinearity(Cand.mod[[36]])
# performance::check_singularity(Cand.mod[[36]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[36]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[36]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[36]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[36]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[36]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[36]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=36, R2=R2[[1]]))



##### Model 37 #####
# -----------------

Cand.mod[[37]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + woven_geotex + levelling
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[37]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[37]])
# performance::check_autocorrelation(Cand.mod[[37]])
# performance::check_collinearity(Cand.mod[[37]])
# performance::check_singularity(Cand.mod[[37]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[37]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[37]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[37]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[37]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[37]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[37]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=37, R2=R2[[1]]))



##### Model 38 #####
# -----------------

Cand.mod[[38]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + woven_geotex + log(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[38]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[38]])
# performance::check_autocorrelation(Cand.mod[[38]])
# performance::check_collinearity(Cand.mod[[38]])
# performance::check_singularity(Cand.mod[[38]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[38]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[38]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[38]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[38]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[38]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[38]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=38, R2=R2[[1]]))



##### Model 39 #####
# ------------------

Cand.mod[[39]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + woven_geotex + tarpfix_pierced
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[39]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[39]])
# performance::check_autocorrelation(Cand.mod[[39]])
# performance::check_collinearity(Cand.mod[[39]])
# performance::check_singularity(Cand.mod[[39]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[39]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[39]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[39]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[39]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[39]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[39]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=39, R2=R2[[1]]))



##### Model 40 #####
# -----------------

Cand.mod[[40]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + woven_geotex + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[40]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[40]])
# performance::check_autocorrelation(Cand.mod[[40]])
# performance::check_collinearity(Cand.mod[[40]])
# performance::check_singularity(Cand.mod[[40]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[40]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[40]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[40]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[40]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[40]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[40]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=40, R2=R2[[1]]))



##### Model 41 #####
# -----------------

Cand.mod[[41]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + tarpfix_pierced + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[41]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[41]])
# performance::check_autocorrelation(Cand.mod[[41]])
# performance::check_collinearity(Cand.mod[[41]])
# performance::check_singularity(Cand.mod[[41]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[41]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[41]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[41]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[41]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[41]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[41]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=41, R2=R2[[1]]))



##### Model 42 #####
# -----------------

Cand.mod[[42]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + tarpfix_pierced + levelling
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[42]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[42]])
# performance::check_autocorrelation(Cand.mod[[42]])
# performance::check_collinearity(Cand.mod[[42]])
# performance::check_singularity(Cand.mod[[42]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[42]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[42]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[42]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[42]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[42]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[42]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=42, R2=R2[[1]]))



##### Model 43 #####
# -----------------

Cand.mod[[43]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + tarpfix_pierced + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[43]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[43]])
# performance::check_autocorrelation(Cand.mod[[43]])
# performance::check_collinearity(Cand.mod[[43]])
# performance::check_singularity(Cand.mod[[43]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[43]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[43]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[43]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[43]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[43]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[43]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=43, R2=R2[[1]]))



##### Model 44 #####
# -----------------

Cand.mod[[44]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + tarpfix_pierced + pb_fixation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[44]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[44]])
# performance::check_autocorrelation(Cand.mod[[44]]) # Close!
# performance::check_collinearity(Cand.mod[[44]])
# performance::check_singularity(Cand.mod[[44]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[44]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[44]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[44]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[44]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[44]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[44]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=44, R2=R2[[1]]))



##### Model 45 #####
# -----------------

Cand.mod[[45]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + stripsoverlap_ok + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[45]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[45]])
# performance::check_autocorrelation(Cand.mod[[45]])
# performance::check_collinearity(Cand.mod[[45]])
# performance::check_singularity(Cand.mod[[45]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[45]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[45]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[45]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[45]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[45]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[45]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=45, R2=R2[[1]]))



##### Model 46 #####
# -----------------

Cand.mod[[46]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + stripsoverlap_ok + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[46]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[46]])
# performance::check_autocorrelation(Cand.mod[[46]])
# performance::check_collinearity(Cand.mod[[46]])
# performance::check_singularity(Cand.mod[[46]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[46]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[46]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[46]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[46]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[46]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[46]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=46, R2=R2[[1]]))



##### Model 47 #####
# -----------------

Cand.mod[[47]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~obstacles + log2(sedicover_height + 1) + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[47]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[47]])
# performance::check_autocorrelation(Cand.mod[[47]])
# performance::check_collinearity(Cand.mod[[47]])
# performance::check_singularity(Cand.mod[[47]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[47]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[47]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[47]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[47]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[47]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[47]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=47, R2=R2[[1]]))



##### Model 48 #####
# -----------------

Cand.mod[[48]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + woven_geotex + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[48]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[48]])
# performance::check_autocorrelation(Cand.mod[[48]]) # Close!
# performance::check_collinearity(Cand.mod[[48]])
# performance::check_singularity(Cand.mod[[48]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[48]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[48]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[48]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[48]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[48]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[48]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=48, R2=R2[[1]]))



##### Model 49 #####
# -----------------

Cand.mod[[49]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + woven_geotex + tarpfix_pierced
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[49]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[49]])
# performance::check_autocorrelation(Cand.mod[[49]]) # close!
# performance::check_collinearity(Cand.mod[[49]])
# performance::check_singularity(Cand.mod[[49]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[49]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[49]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[49]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[49]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[49]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[49]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=49, R2=R2[[1]]))



##### Model 50 #####
# -----------------

Cand.mod[[50]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + woven_geotex + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[50]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[50]])
# performance::check_autocorrelation(Cand.mod[[50]]) # Nope!
# performance::check_collinearity(Cand.mod[[50]])
# performance::check_singularity(Cand.mod[[50]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[50]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[50]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[50]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[50]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[50]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[50]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=50, R2=R2[[1]]))



##### Model 51 #####
# -----------------

Cand.mod[[51]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarpfix_pierced + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[51]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[51]])
# performance::check_autocorrelation(Cand.mod[[51]]) # Close!
# performance::check_collinearity(Cand.mod[[51]])
# performance::check_singularity(Cand.mod[[51]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[51]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[51]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[51]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[51]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[51]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[51]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=51, R2=R2[[1]]))



##### Model 52 #####
# -----------------

Cand.mod[[52]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarpfix_pierced + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[52]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[52]])
# performance::check_autocorrelation(Cand.mod[[52]]) # Close!
# performance::check_collinearity(Cand.mod[[52]])
# performance::check_singularity(Cand.mod[[52]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[52]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[52]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[52]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[52]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[52]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[52]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=52, R2=R2[[1]]))



##### Model 53 #####
# -----------------

Cand.mod[[53]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarpfix_pierced + pb_fixation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[53]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[53]])
# performance::check_autocorrelation(Cand.mod[[53]]) # Nope!
# performance::check_collinearity(Cand.mod[[53]])
# performance::check_singularity(Cand.mod[[53]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[53]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[53]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[53]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[53]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[53]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[53]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=53, R2=R2[[1]]))



##### Model 54 #####
# -----------------

Cand.mod[[54]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + tarpfix_pierced + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[54]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[54]])
# performance::check_autocorrelation(Cand.mod[[54]]) # Close!
# performance::check_collinearity(Cand.mod[[54]])
# performance::check_singularity(Cand.mod[[54]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[54]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[54]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[54]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[54]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[54]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[54]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=54, R2=R2[[1]]))



##### Model 55 #####
# -----------------

Cand.mod[[55]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + plantation + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[55]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[55]])
# performance::check_autocorrelation(Cand.mod[[55]]) # Close!
# performance::check_collinearity(Cand.mod[[55]])
# performance::check_singularity(Cand.mod[[55]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[55]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[55]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[55]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[55]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[55]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[55]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=55, R2=R2[[1]]))



##### Model 56 #####
# -----------------

Cand.mod[[56]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + plantation + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[56]], n = 1000, plot = FALSE)
# plot(simu.resid) # NOPE! Same as before (except here, including an additional random effect worsen things)
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[56]])
# performance::check_autocorrelation(Cand.mod[[56]]) # Close!
# performance::check_collinearity(Cand.mod[[56]])
# performance::check_singularity(Cand.mod[[56]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[56]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[56]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[56]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[56]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[56]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[56]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=56, R2=R2[[1]]))



##### Model 57 #####
# -----------------

Cand.mod[[57]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + log2(stand_surface) + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[57]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[57]])
# performance::check_autocorrelation(Cand.mod[[57]]) # Close!
# performance::check_collinearity(Cand.mod[[57]])
# performance::check_singularity(Cand.mod[[57]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[57]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[57]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[57]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[57]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[57]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[57]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=57, R2=R2[[1]]))



##### Model 58 #####
# -----------------

Cand.mod[[58]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + log2(stand_surface) + pb_fixation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[58]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[58]])
# performance::check_autocorrelation(Cand.mod[[58]]) # Nope! Same as before.
# performance::check_collinearity(Cand.mod[[58]])
# performance::check_singularity(Cand.mod[[58]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[58]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[58]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[58]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[58]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[58]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[58]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=58, R2=R2[[1]]))



##### Model 59 #####
# -----------------

Cand.mod[[59]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~add_control + log2(stand_surface) + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[59]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$add_control) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[59]])
# performance::check_autocorrelation(Cand.mod[[59]]) # Nope!
# performance::check_collinearity(Cand.mod[[59]])
# performance::check_singularity(Cand.mod[[59]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[59]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[59]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[59]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[59]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[59]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[59]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=59, R2=R2[[1]]))



##### Model 60 #####
# -----------------

Cand.mod[[60]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + tarpfix_pierced + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[60]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[60]])
# performance::check_autocorrelation(Cand.mod[[60]]) # Close!
# performance::check_collinearity(Cand.mod[[60]])
# performance::check_singularity(Cand.mod[[60]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[60]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[60]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[60]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[60]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[60]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[60]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=60, R2=R2[[1]]))



##### Model 60 #####
# -----------------

Cand.mod[[61]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + tarpfix_pierced + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[61]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[61]])
# performance::check_autocorrelation(Cand.mod[[61]]) # Close!
# performance::check_collinearity(Cand.mod[[61]])
# performance::check_singularity(Cand.mod[[61]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[61]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[61]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[61]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[61]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[61]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[61]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=61, R2=R2[[1]]))



##### Model 62 #####
# -----------------

Cand.mod[[62]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(stand_surface) + tarpfix_pierced + plantation
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[62]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$plantation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[62]])
# performance::check_autocorrelation(Cand.mod[[62]]) # Close!
# performance::check_collinearity(Cand.mod[[62]])
# performance::check_singularity(Cand.mod[[62]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[62]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[62]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[62]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[62]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[62]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[62]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=62, R2=R2[[1]]))



##### Model 63 #####
# -----------------

Cand.mod[[63]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(stand_surface) + tarpfix_pierced + stripsoverlap_ok
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[63]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[63]])
# performance::check_autocorrelation(Cand.mod[[63]])
# performance::check_collinearity(Cand.mod[[63]])
# performance::check_singularity(Cand.mod[[63]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[63]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[63]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[63]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[63]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[63]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[63]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=63, R2=R2[[1]]))



##### Model 64 #####
# -----------------

Cand.mod[[64]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(stand_surface) + tarpfix_pierced + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[64]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarpfix_pierced) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[64]])
# performance::check_autocorrelation(Cand.mod[[64]]) # Close!
# performance::check_collinearity(Cand.mod[[64]])
# performance::check_singularity(Cand.mod[[64]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[64]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[64]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[64]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[64]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[64]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[64]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=64, R2=R2[[1]]))



##### Model 65 #####
# -----------------

Cand.mod[[65]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + log2(sedicover_height + 1) + log2(stand_surface)
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[65]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = rtarped$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[65]])
# performance::check_autocorrelation(Cand.mod[[65]]) # Close!
# performance::check_collinearity(Cand.mod[[65]])
# performance::check_singularity(Cand.mod[[65]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[65]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[65]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[65]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[65]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[65]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[65]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=65, R2=R2[[1]]))



##### Model 66 #####
# -----------------

Cand.mod[[66]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + log2(sedicover_height + 1) + stripsoverlap_ok
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[66]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$stripsoverlap_ok) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[66]])
# performance::check_autocorrelation(Cand.mod[[66]])
# performance::check_collinearity(Cand.mod[[66]])
# performance::check_singularity(Cand.mod[[66]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[66]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[66]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[66]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[66]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[66]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[66]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=66, R2=R2[[1]]))



##### Model 67 #####
# -----------------

Cand.mod[[67]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~woven_geotex + log2(sedicover_height + 1) + tarping_duration
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[67]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok!
# DHARMa::plotResiduals(simu.resid, form = rtarped$sedicover_height) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope! Same as before.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[67]])
# performance::check_autocorrelation(Cand.mod[[67]])
# performance::check_collinearity(Cand.mod[[67]])
# performance::check_singularity(Cand.mod[[67]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[67]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[67]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[67]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[67]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[67]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[67]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=67, R2=R2[[1]]))



##### Model 68 #####
# -----------------

Cand.mod[[68]] <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~log2(sedicover_height + 1) + tarping_duration + geomem
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[68]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$woven_geotex) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = rtarped$tarping_duration) # Nope!
# DHARMa::plotResiduals(simu.resid, form = rtarped$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = rtarped$longitude, y = rtarped$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[68]])
# performance::check_autocorrelation(Cand.mod[[68]]) # Close!
# performance::check_collinearity(Cand.mod[[68]])
# performance::check_singularity(Cand.mod[[68]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[68]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[68]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[68]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[68]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[68]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_nakagawa(Cand.mod[[68]]) # Nakagawa's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=68, R2=R2[[1]]))




###### hahahahahahahahaha #######


# Variations de 41
rrr <- glmmTMB::glmmTMB(formula = lreg_tarpedarea~stripsoverlap_ok + tarpfix_multimethod + slope
                                   + (1|xp_id), data = rtarped,
                                   family = stats::binomial(link = "logit"), REML = FALSE)
summary(rrr)
summary(Cand.mod[[51]])


"levelling"  "obstacles" "pb_durability"  "pb_fixation"   "repairs"
"slope"  "stripsoverlap_ok"
[25] "tarpfix_multimethod"
 "woven_geotex"

###################################
###################################
###################################
###################################




# ------------------------------------- #
##### Model selection and averaging #####
# ------------------------------------- #

Model <- (1:50)

Candidate <- c("null",
               "log2(distance + 1) + add_control", "log2(distance + 1) + obstacles",
               "log2(distance + 1) + repairs", "log2(distance + 1) + pb_fixation", "log2(distance + 1) + log2(stand_surface)",
               "log2(distance + 1) + geomem", "log2(distance + 1) + tarping_duration", "log2(distance + 1) + log2(trench_depth + 1)",
               "log2(distance + 1) * log2(trench_depth + 1)", "log2(trench_depth + 1) + add_control", "log2(trench_depth + 1) + obstacles",
               "log2(trench_depth + 1) + log2(stand_surface)", "log2(stand_surface) + add_control", "log2(stand_surface) + obstacles",
               "log2(stand_surface) + pb_fixation", "log2(stand_surface) + tarping_duration", "log2(trench_depth + 1) + repairs",
               "log2(trench_depth + 1) + slope", "add_control + obstacles", "add_control + pb_fixation",
               "add_control + slope", "add_control + tarping_duration", "obstacles + pb_fixation",
               "obstacles + geomem", "obstacles + repairs", "slope + geomem", "slope + pb_fixation",
               "slope + tarping_duration", "log2(distance + 1) + log2(trench_depth + 1) + add_control",
               "log2(distance + 1) + log2(trench_depth + 1) + obstacles", "log2(distance + 1) + log2(trench_depth + 1) + geomem",
               "log2(distance + 1) + log2(trench_depth + 1) + log2(stand_surface)", "log2(distance + 1) + log2(trench_depth + 1) + slope",
               "log2(distance + 1) + log2(stand_surface) + add_control", "log2(distance + 1) + log2(stand_surface) + obstacles",
               "log2(distance + 1) + log2(stand_surface) + pb_fixation", "log2(distance + 1) + geomem + tarping_duration",
               "log2(stand_surface) + log2(trench_depth + 1) + obstacles", "log2(stand_surface) + log2(trench_depth + 1) + slope",
               "tarping_duration + log2(distance + 1) + add_control", "obstacles + repairs + add_control",
               "obstacles + add_control + geomem", "obstacles + add_control + tarping_duration",
               "obstacles + add_control + log2(trench_depth + 1)", "obstacles + add_control + log2(stand_surface)",
               "obstacles + add_control + slope", "log2(trench_depth + 1) + log2(distance + 1) + repairs",
               "log2(trench_depth + 1) + obstacles + repairs", "log2(trench_depth + 1) + obstacles + geomem")

Cand.model <- data.frame(Model, Candidate)



##### Model selection #####
# -------------------------

### Rank models based on AICc:
AICc <- MuMIn::model.sel(object = Cand.mod)

### Improved formating:
AICc.model <- as.data.frame(AICc)
AICc.model$csweigth <- cumsum(AICc.model$weight) # Add a column containing the cumulative sum of AICc weights
AICc.model$Model <- row.names(AICc.model) # Add a column containing the n° of each candidate model

AICc.model <- merge(AICc.model, Cand.model, by="Model") # CAUTION: this merger reorder the rows of the
# the table based on Model, and it removes the NULL model from the table (as it doesn't exist in Cand.model).
# Consequently, the new ordering begins at 10 (10, 11, 12, ..., 20, 21, ...)!

AICc.model <- merge(AICc.model, R.ajust, by="Model") # This merger adds the computed R2 (here: pseudo-R2)

AICc.model <- AICc.model[order(AICc.model$delta),] # To reorder rows according to delta AICc
AICc.model$Response <- "lreg_tarpedarea"
AICc.model$Rank <- 1:nrow(AICc.model)
AICc.model <- AICc.model[,c("Rank", "Model", "Response", "Candidate", "df", "AICc", "delta", "weight", "R2")]
AICc.model$AICc <- format(round(AICc.model$AICc, digits=1))
AICc.model$delta <- format(round(AICc.model$delta, digits=3))
AICc.model$weight <- format(round(AICc.model$weight, digits=3))
AICc.model$R2 <- format(round(AICc.model$R2, digits=3))
colnames(AICc.model)[colnames(AICc.model) == 'Candidate'] <- 'Candidate model'
colnames(AICc.model)[colnames(AICc.model) == 'df'] <- 'k'
colnames(AICc.model)[colnames(AICc.model) == 'delta'] <- 'delta AICc'
colnames(AICc.model)[colnames(AICc.model) == 'weight'] <- 'W'
colnames(AICc.model)[colnames(AICc.model) == 'R2'] <- 'R2glmm'

### Table export:
readr::write_csv2(x = AICc.model, file = here::here("output", "tables", "Models_lreg_tarpedarea.csv"))



##### Multimodel Inference (averaging) #####
# ------------------------------------------

Parameters <- c("Intercept", "add_control", "log2(distance + 1)", "geomem", "obstacles",
                "pb_fixation", "slope", "log2(stand_surface)", "tarping_duration",
                "log2(trench_depth + 1)", "repairs")
Var <- c("cond((Int))", "cond(add_control1)", "cond(log2(distance + 1))", "cond(geomem1)",
         "cond(obstacles)", "cond(pb_fixation1)", "cond(slope)", "cond(log2(stand_surface))",
         "cond(tarping_duration)", "cond(log2(trench_depth + 1))","cond(repairs1)")
Var.Imp <- c("cond((Int))", "cond(add_control)", "cond(log2(distance + 1))", "cond(geomem)",
             "cond(obstacles)", "cond(pb_fixation)", "cond(slope)", "cond(log2(stand_surface))",
             "cond(tarping_duration)", "cond(log2(trench_depth + 1))","cond(repairs)")

Para.model <- data.frame(Parameters, Var, Var.Imp)

### Select the top models:
# Removing the model including an interaction before averaging:
AICc.2 <- AICc[-1,]
#top.models <- MuMIn::get.models(AICc.2, cumsum(weight) <= 0.95) # To take those with a cumulated sum of
# AICc weights <= 0.95
top.models <- MuMIn::get.models(AICc.2, cumsum(weight) <= 1) # To take them all
# We could also select models according to their delta AICc (see Burnham & Anderson, 2002)!

### Actual model parameters averaging:
Parameter <- MuMIn::model.avg(top.models, revised.var=T, adjusted=T, fit=T)
Parameter.model <- as.data.frame(cbind(MuMIn::coefTable(Parameter), stats::confint(Parameter)))

### Improved formating:
Parameter.model$Var <- row.names(Parameter.model)
Parameter.model <- merge(Parameter.model, Para.model, by="Var", all=TRUE)
Imp <- as.data.frame(format(round(MuMIn::importance(Parameter), digits=2)))
colnames(Imp) <- "Imp."
Imp$Var.Imp <- row.names(Imp)
Parameter.model <- merge(Parameter.model, Imp, by="Var.Imp", all=TRUE)

Parameter.model$'Estimate (±SE)' <- paste0(format(round(Parameter.model$Estimate, digits=3), trim=T),
                                           " (±", format(round(Parameter.model$'Std. Error', digits=3),
                                                         trim=T), ")")
Parameter.model$'(95% CI)' <- paste0("(", format(round(Parameter.model$'2.5 %', digits=3), trim=T),
                                     "; ", format(round(Parameter.model$'97.5 %', digits=3), trim=T), ")")
Parameter.model$Response <- "lreg_tarpedarea"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Response", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)",
                                      'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Imp.), ]

### Table export:
readr::write_csv2(x = Parameter.model, file = here::here("output", "tables", "Parameters_lreg_tarpedarea.csv"))
