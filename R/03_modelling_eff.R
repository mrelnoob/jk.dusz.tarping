################################
##### MODELLING efficiency #####
################################

### IMPORTANT NOTE: to avoid creating notes for unquoted variables, I must add the following code at
# the beginning of every source file (e.g. R/myscript.R) that uses unquoted variables, so in front of
# all my scripts doing any kind of analyses (otherwise, I should always assign each variable, e.g.
# mydata$myvariable, which is quite time consuming and wearisome).
# if(getRversion() >= "2.15.1")  utils::globalVariables(c(
#   "manager_id", "xp_id", "country", "efficiency", "eff_eradication", "high_eff",
#   "latitude", "elevation", "goals", "slope", "coarse_env", "obstacles", "flood",
#   "geomem", "maxveg", "uprootexcav", "stand_surface", "fully_tarped", "distance", "tarping_duration",
#   "stripsoverlap_ok", "tarpfix_pierced", "tarpfix_multimethod", "sedicover_height",
#   "trench_depth", "plantation", "followups",
#   "pb_fixation", "pb_durability"))
# -------------------------------------------------------------------------------------------------------- #





# --------------------------------------------------- #
##### Data preparation for modelling 'efficiency' #####
# --------------------------------------------------- #

# List of used packages (for publication or package building): here, readr, (gamlss, gamlss.dist), MuMIn,
# emmeans, glmmTMB, ggplot2, broom.mixed

.pardefault <- par() # To save the default graphical parameters (in case I want to restore them).

library(jk.dusz.tarping)
readr::read_csv(here::here("mydata", "erad.csv"), col_names = TRUE, col_types =
                         readr::cols(
                           manager_id = readr::col_factor(),
                           xp_id = readr::col_factor(),
                           eff_eradication = readr::col_factor(c("0", "1")),
                           high_eff = readr::col_factor(c("0", "1")),
                           goals = readr::col_factor(),
                           geomem = readr::col_factor(c("0", "1")),
                           maxveg = readr::col_factor(c("0", "1")),
                           uprootexcav = readr::col_factor(c("0", "1")),
                           fully_tarped = readr::col_factor(c("0", "1")),
                           stripsoverlap_ok = readr::col_factor(c("0", "1")),
                           tarpfix_multimethod = readr::col_factor(c("0", "1")),
                           tarpfix_pierced = readr::col_factor(c("0", "1")),
                           plantation = readr::col_factor(c("0", "1")),
                           pb_fixation = readr::col_factor(c("0", "1")),
                           pb_durability = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(efficiency = efficiency/10) %>% # For 'efficiency' to be look like a percentage that could
  # be modelled using a Beta Regression Model.
  dplyr::mutate(efficiency = ifelse(efficiency == 1, 0.999, efficiency)) -> eff
summary(eff)



### Custom functions for modelling:
# To help creating QQplots for beta regression models:

qq.line = function(x) {
  # following four lines from base R's qqline()
  y <- quantile(x[!is.na(x)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  return(c(int = int, slope = slope))
}

# Since the MuMIn::r.squaredGLMM() function does not work for beta models with logit links, here is a way
# to compute (marginal and conditional) Pseudo-R2 based on: Johnson, P.C.D. (2014) Extension of Nakagawa &
# Schielzeth’s R_GLMM² to random slopes models. Methods in Ecology and Evolution 5: 44-946:

pseudo.R2glmm <- function(model){
  X <- model.matrix(model)
  n <- nrow(X)
  Beta <- glmmTMB::fixef(model)$cond
  Sf <- var(X %*% Beta)
  Sigma.list <- glmmTMB::VarCorr(model)
  Sl <-
    sum(
      sapply(Sigma.list$cond,
             function(Sigma)
             {
               Z <-X[,rownames(Sigma)]
               sum(diag(Z %*% Sigma %*% t(Z)))/n
             }))
  Se <- attr(Sigma.list$cond, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  Marginal_R2 <- Sf / total.var # Marginal R_GLMM² represents the variance explained by the fixed effects
  Conditional_R2 <- (Sf + Sl) / total.var # Conditional R_GLMM² represents the variance explained by the entire model
  R2 <- data.frame(Marginal_R2, Conditional_R2)
  print(R2)
}





# ---------------------------------------- #
##### Building of the candidate models #####
# ---------------------------------------- #

Cand.mod <- list()
Cand.mod[[1]] <- gamlss::gamlss(efficiency~distance, data = eff, family = BEOI())



##### Test model with GAMLSS #####
# --------------------------------

# library(gamlss)
# mod1 <- gamlss(formula = efficiency~distance+re(random = ~1|manager_id),
#                        nu.formula = ~distance+re(random = ~1|manager_id),
#                        data = eff, family = BEINF())
#
# # Checking the residuals:
# plot(mod1) # Obs. n°67 is probably problematic?
# plot(fitted(mod1)~eff$efficiency)
# # Checking the goodness-of-fit and overdispersion:
# dat.resid <- sum(resid(mod1, type = "weighted", what = "mu")^2)
# 1 - pchisq(dat.resid, mod1$df.resid) # Conclusions: Pearson residuals do indicate a lack of fit (p values
# # far lower than 0.05)
# dat.resid/df.residual(mod1) # Conclusions: the data are definitely not overdispersed.
# # Exploring parameters and hypotheses:
# summary(mod1)
# binomial()$linkinv(coef(mod1))  # logit inverse: to get interpretable parameter estimates
# # Further trend exploration:
# gamlss::Rsq(mod1) # R^2 seems too good to be true!



##### Model 1 (null model) #####
# ------------------------------

Cand.mod[[1]] <- glmmTMB::glmmTMB(formula = efficiency~(1|manager_id), data = eff,
                         family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[1]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[1]])))
QQline <- qq.line(resid(Cand.mod[[1]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[1]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[1]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[1]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[1]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[1]])) # Ok
dat.resid/df.residual(Cand.mod[[1]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[1]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[1]])
family(Cand.mod[[1]])$linkinv(glmmTMB::fixef(Cand.mod[[1]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[1]])
broom.mixed::glance(Cand.mod[[1]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[1]], ~1))[1] - logLik(Cand.mod[[1]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[1]])



##### Model 2 #####
# -----------------

Cand.mod[[2]] <- glmmTMB::glmmTMB(formula = efficiency~distance + (1|manager_id), data = eff,
                         family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[2]],
                                                   type = "pearson"), x = fitted(Cand.mod[[2]])))
QQline <- qq.line(resid(Cand.mod[[2]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[2]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[2]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[2]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[2]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[2]])) # Ok
dat.resid/df.residual(Cand.mod[[2]]) # Ok
# To plot the observed vs. fitted values:
par(.pardefault) # To restore defaults graphical parameters
plot(x = fitted(Cand.mod[[2]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[2]])
family(Cand.mod[[2]])$linkinv(glmmTMB::fixef(Cand.mod[[2]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[2]])
broom.mixed::glance(Cand.mod[[2]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[2]], ~1))[1] - logLik(Cand.mod[[2]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[2]])



##### Model 3 #####
# -----------------

Cand.mod[[3]] <- glmmTMB::glmmTMB(formula = efficiency~followups + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[3]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[3]])))
QQline <- qq.line(resid(Cand.mod[[3]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[3]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[3]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[3]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[3]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[3]])) # Ok
dat.resid/df.residual(Cand.mod[[3]]) # Ok
# To plot the observed vs. fitted values:
par(.pardefault) # To restore defaults graphical parameters
plot(x = fitted(Cand.mod[[3]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[3]])
family(Cand.mod[[3]])$linkinv(glmmTMB::fixef(Cand.mod[[3]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[3]])
broom.mixed::glance(Cand.mod[[3]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[3]], ~1))[1] - logLik(Cand.mod[[3]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[3]])



##### Model 4 #####
# -----------------

Cand.mod[[4]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[4]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[4]])))
QQline <- qq.line(resid(Cand.mod[[4]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[4]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[4]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[4]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[4]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[4]])) # Ok
dat.resid/df.residual(Cand.mod[[4]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[4]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[4]])
family(Cand.mod[[4]])$linkinv(glmmTMB::fixef(Cand.mod[[4]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[4]])
broom.mixed::glance(Cand.mod[[4]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[4]], ~1))[1] - logLik(Cand.mod[[4]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[4]])



##### Model 5 #####
# -----------------

Cand.mod[[5]] <- glmmTMB::glmmTMB(formula = efficiency~geomem + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[5]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[5]])))
QQline <- qq.line(resid(Cand.mod[[5]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[5]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[5]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[5]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[5]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[5]])) # Ok
dat.resid/df.residual(Cand.mod[[5]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[5]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[5]])
family(Cand.mod[[5]])$linkinv(glmmTMB::fixef(Cand.mod[[5]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[5]])
broom.mixed::glance(Cand.mod[[5]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[5]], ~1))[1] - logLik(Cand.mod[[5]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[5]])



##### Model 6 #####
# -----------------

Cand.mod[[6]] <- glmmTMB::glmmTMB(formula = efficiency~obstacles + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[6]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[6]])))
QQline <- qq.line(resid(Cand.mod[[6]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[6]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[6]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[6]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[6]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[6]])) # Ok
dat.resid/df.residual(Cand.mod[[6]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[6]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[6]])
family(Cand.mod[[6]])$linkinv(glmmTMB::fixef(Cand.mod[[6]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[6]])
broom.mixed::glance(Cand.mod[[6]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[6]], ~1))[1] - logLik(Cand.mod[[6]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[6]])



##### Model 7 #####
# -----------------

Cand.mod[[7]] <- glmmTMB::glmmTMB(formula = efficiency~plantation + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[7]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[7]])))
QQline <- qq.line(resid(Cand.mod[[7]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[7]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[7]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[7]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[7]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[7]])) # Ok
dat.resid/df.residual(Cand.mod[[7]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[7]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[7]])
family(Cand.mod[[7]])$linkinv(glmmTMB::fixef(Cand.mod[[7]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[7]])
broom.mixed::glance(Cand.mod[[7]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[7]], ~1))[1] - logLik(Cand.mod[[7]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[7]])



##### Model 8 #####
# -----------------

Cand.mod[[8]] <- glmmTMB::glmmTMB(formula = efficiency~pb_fixation + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[8]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[8]])))
QQline <- qq.line(resid(Cand.mod[[8]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[8]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[8]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[8]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[8]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[8]])) # Ok
dat.resid/df.residual(Cand.mod[[8]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[8]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[8]])
family(Cand.mod[[8]])$linkinv(glmmTMB::fixef(Cand.mod[[8]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[8]])
broom.mixed::glance(Cand.mod[[8]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[8]], ~1))[1] - logLik(Cand.mod[[8]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[8]])



##### Model 9 #####
# -----------------

Cand.mod[[9]] <- glmmTMB::glmmTMB(formula = efficiency~sedicover_height + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[9]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[9]])))
QQline <- qq.line(resid(Cand.mod[[9]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[9]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[9]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[9]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[9]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[9]])) # Ok
dat.resid/df.residual(Cand.mod[[9]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[9]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[9]])
family(Cand.mod[[9]])$linkinv(glmmTMB::fixef(Cand.mod[[9]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[9]])
broom.mixed::glance(Cand.mod[[9]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[9]], ~1))[1] - logLik(Cand.mod[[9]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[9]])



##### Model 10 #####
# -----------------

Cand.mod[[10]] <- glmmTMB::glmmTMB(formula = efficiency~slope + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[10]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[10]])))
QQline <- qq.line(resid(Cand.mod[[10]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[10]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[10]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[10]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[10]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[10]])) # Ok
dat.resid/df.residual(Cand.mod[[10]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[10]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[10]])
family(Cand.mod[[10]])$linkinv(glmmTMB::fixef(Cand.mod[[10]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[10]])
broom.mixed::glance(Cand.mod[[10]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[10]], ~1))[1] - logLik(Cand.mod[[10]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[10]])



##### Model 11 #####
# -----------------

Cand.mod[[11]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + (1|manager_id), data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[11]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[11]])))
QQline <- qq.line(resid(Cand.mod[[11]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[11]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[11]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[11]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[11]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[11]])) # Ok
dat.resid/df.residual(Cand.mod[[11]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[11]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[11]])
family(Cand.mod[[11]])$linkinv(glmmTMB::fixef(Cand.mod[[11]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[11]])
broom.mixed::glance(Cand.mod[[11]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[11]], ~1))[1] - logLik(Cand.mod[[11]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[11]])




##### Model 12 #####
# -----------------

Cand.mod[[12]] <- glmmTMB::glmmTMB(formula = efficiency~tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[12]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[12]])))
QQline <- qq.line(resid(Cand.mod[[12]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[12]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[12]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[12]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[12]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[12]])) # Ok
dat.resid/df.residual(Cand.mod[[12]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[12]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[12]])
family(Cand.mod[[12]])$linkinv(glmmTMB::fixef(Cand.mod[[12]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[12]])
broom.mixed::glance(Cand.mod[[12]])

### Computing a quasiR^2:
1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[12]], ~1))[1] - logLik(Cand.mod[[12]])[1])) # Methods from flutterbys.com
pseudo.R2glmm(model = Cand.mod[[12]])




##### Model 13 #####
# -----------------

Cand.mod[[13]] <- glmmTMB::glmmTMB(formula = efficiency~uprootexcav + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[13]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[13]])))
QQline <- qq.line(resid(Cand.mod[[13]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[13]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[13]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[13]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[13]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[13]])) # Ok
dat.resid/df.residual(Cand.mod[[13]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[13]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[13]])
family(Cand.mod[[13]])$linkinv(glmmTMB::fixef(Cand.mod[[13]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[13]])
broom.mixed::glance(Cand.mod[[13]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[13]])




##### Model 14 #####
# -----------------

Cand.mod[[14]] <- glmmTMB::glmmTMB(formula = efficiency~distance + followups + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[14]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[14]])))
QQline <- qq.line(resid(Cand.mod[[14]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[14]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[14]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[14]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[14]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[14]])) # Ok
dat.resid/df.residual(Cand.mod[[14]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[14]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[14]])
family(Cand.mod[[14]])$linkinv(glmmTMB::fixef(Cand.mod[[14]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[14]])
broom.mixed::glance(Cand.mod[[14]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[14]])




##### Model 15 #####
# -----------------

Cand.mod[[15]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[15]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[15]])))
QQline <- qq.line(resid(Cand.mod[[15]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[15]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[15]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[15]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[15]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[15]])) # Ok
dat.resid/df.residual(Cand.mod[[15]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[15]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[15]])
family(Cand.mod[[15]])$linkinv(glmmTMB::fixef(Cand.mod[[15]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[15]])
broom.mixed::glance(Cand.mod[[15]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[15]])



##### Model 16 #####
# -----------------

Cand.mod[[16]] <- glmmTMB::glmmTMB(formula = efficiency~distance + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[16]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[16]])))
QQline <- qq.line(resid(Cand.mod[[16]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[16]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[16]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[16]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[16]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[16]])) # Ok
dat.resid/df.residual(Cand.mod[[16]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[16]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[16]])
family(Cand.mod[[16]])$linkinv(glmmTMB::fixef(Cand.mod[[16]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[16]])
broom.mixed::glance(Cand.mod[[16]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[16]])



##### Model 17 #####
# -----------------

Cand.mod[[17]] <- glmmTMB::glmmTMB(formula = efficiency~distance + log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[17]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[17]])))
QQline <- qq.line(resid(Cand.mod[[17]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[17]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[17]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[17]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[17]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[17]])) # Ok
dat.resid/df.residual(Cand.mod[[17]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[17]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[17]])
family(Cand.mod[[17]])$linkinv(glmmTMB::fixef(Cand.mod[[17]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[17]])
broom.mixed::glance(Cand.mod[[17]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[17]])



##### Model 18 #####
# -----------------

Cand.mod[[18]] <- glmmTMB::glmmTMB(formula = efficiency~distance * log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[18]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[18]])))
QQline <- qq.line(resid(Cand.mod[[18]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[18]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[18]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[18]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[18]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[18]])) # Ok
dat.resid/df.residual(Cand.mod[[18]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[18]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[18]])
family(Cand.mod[[18]])$linkinv(glmmTMB::fixef(Cand.mod[[18]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[18]])
broom.mixed::glance(Cand.mod[[18]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[18]])



