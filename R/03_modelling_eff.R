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
                           repairs = readr::col_factor(c("0", "1")),
                           add_control = readr::col_factor(c("0", "1")),
                           pb_fixation = readr::col_factor(c("0", "1")),
                           pb_durability = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(efficiency = efficiency/10) %>% # For 'efficiency' to be look like a percentage that could
  # be modelled using a Beta Regression Model.
  dplyr::mutate(efficiency = ifelse(efficiency == 1, 0.999, efficiency)) -> eff
summary(eff)



### Custom functions for modelling:
# To help creating QQ-plots for beta regression models:

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


§§§§§§§ TESTER SI ça change quand stand_surface/10!!!!!

  §§§§§§§§
++++ centering pour interactions!!!
  ++ tester R2 de nagelkerke
++ VIF
++ autres checks de performance::!!!!!!
  ++ tester avec repairs/add_control§§§
++ Huber-White SE!!!!!!!
---- de stand_surface, followups et ++++ de repairs/add_control et d'obstacles!'



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
performance::r2_nakagawa(model = Cand.mod[[2]])



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



##### Model 19 #####
# -----------------

Cand.mod[[19]] <- glmmTMB::glmmTMB(formula = efficiency~distance + uprootexcav + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[19]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[19]])))
QQline <- qq.line(resid(Cand.mod[[19]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[19]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[19]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[19]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[19]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[19]])) # Ok
dat.resid/df.residual(Cand.mod[[19]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[19]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[19]])
family(Cand.mod[[19]])$linkinv(glmmTMB::fixef(Cand.mod[[19]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[19]])
broom.mixed::glance(Cand.mod[[19]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[19]])



##### Model 20 #####
# -----------------

Cand.mod[[20]] <- glmmTMB::glmmTMB(formula = efficiency~followups + fully_tarped + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[20]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[20]])))
QQline <- qq.line(resid(Cand.mod[[20]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[20]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[20]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[20]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[20]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[20]])) # Ok
dat.resid/df.residual(Cand.mod[[20]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[20]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[20]])
family(Cand.mod[[20]])$linkinv(glmmTMB::fixef(Cand.mod[[20]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[20]])
broom.mixed::glance(Cand.mod[[20]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[20]])



##### Model 21 #####
# -----------------

Cand.mod[[21]] <- glmmTMB::glmmTMB(formula = efficiency~followups + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[21]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[21]])))
QQline <- qq.line(resid(Cand.mod[[21]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[21]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[21]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[21]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[21]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[21]])) # Ok
dat.resid/df.residual(Cand.mod[[21]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[21]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[21]])
family(Cand.mod[[21]])$linkinv(glmmTMB::fixef(Cand.mod[[21]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[21]])
broom.mixed::glance(Cand.mod[[21]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[21]])



##### Model 22 #####
# -----------------

Cand.mod[[22]] <- glmmTMB::glmmTMB(formula = efficiency~followups * pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[22]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[22]])))
QQline <- qq.line(resid(Cand.mod[[22]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[22]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[22]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[22]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[22]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[22]])) # Ok
dat.resid/df.residual(Cand.mod[[22]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[22]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[22]])
family(Cand.mod[[22]])$linkinv(glmmTMB::fixef(Cand.mod[[22]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[22]])
broom.mixed::glance(Cand.mod[[22]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[22]])



##### Model 23 #####
# -----------------

Cand.mod[[23]] <- glmmTMB::glmmTMB(formula = efficiency~followups + log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[23]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[23]])))
QQline <- qq.line(resid(Cand.mod[[23]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[23]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[23]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[23]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[23]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[23]])) # Ok
dat.resid/df.residual(Cand.mod[[23]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[23]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[23]])
family(Cand.mod[[23]])$linkinv(glmmTMB::fixef(Cand.mod[[23]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[23]])
broom.mixed::glance(Cand.mod[[23]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[23]])



##### Model 24 #####
# -----------------

Cand.mod[[24]] <- glmmTMB::glmmTMB(formula = efficiency~followups + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[24]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[24]])))
QQline <- qq.line(resid(Cand.mod[[24]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[24]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[24]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[24]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[24]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[24]])) # Ok
dat.resid/df.residual(Cand.mod[[24]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[24]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[24]])
family(Cand.mod[[24]])$linkinv(glmmTMB::fixef(Cand.mod[[24]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[24]])
broom.mixed::glance(Cand.mod[[24]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[24]])



##### Model 25 #####
# -----------------

Cand.mod[[25]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[25]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[25]])))
QQline <- qq.line(resid(Cand.mod[[25]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[25]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[25]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[25]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[25]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[25]])) # Ok
dat.resid/df.residual(Cand.mod[[25]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[25]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[25]])
family(Cand.mod[[25]])$linkinv(glmmTMB::fixef(Cand.mod[[25]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[25]])
broom.mixed::glance(Cand.mod[[25]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[25]])



##### Model 26 #####
# -----------------

Cand.mod[[26]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[26]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[26]])))
QQline <- qq.line(resid(Cand.mod[[26]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[26]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[26]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[26]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[26]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[26]])) # Ok
dat.resid/df.residual(Cand.mod[[26]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[26]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[26]])
family(Cand.mod[[26]])$linkinv(glmmTMB::fixef(Cand.mod[[26]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[26]])
broom.mixed::glance(Cand.mod[[26]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[26]])



##### Model 27 #####
# -----------------

Cand.mod[[27]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + sedicover_height + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[27]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[27]])))
QQline <- qq.line(resid(Cand.mod[[27]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[27]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[27]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[27]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[27]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[27]])) # Ok
dat.resid/df.residual(Cand.mod[[27]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[27]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[27]])
family(Cand.mod[[27]])$linkinv(glmmTMB::fixef(Cand.mod[[27]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[27]])
broom.mixed::glance(Cand.mod[[27]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[27]])



##### Model 28 #####
# -----------------

Cand.mod[[28]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[28]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[28]])))
QQline <- qq.line(resid(Cand.mod[[28]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[28]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[28]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[28]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[28]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[28]])) # Ok
dat.resid/df.residual(Cand.mod[[28]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[28]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[28]])
family(Cand.mod[[28]])$linkinv(glmmTMB::fixef(Cand.mod[[28]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[28]])
broom.mixed::glance(Cand.mod[[28]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[28]])



##### Model 29 #####
# -----------------

Cand.mod[[29]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped * log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[29]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[29]])))
QQline <- qq.line(resid(Cand.mod[[29]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[29]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[29]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[29]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[29]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[29]])) # Ok
dat.resid/df.residual(Cand.mod[[29]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[29]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[29]])
family(Cand.mod[[29]])$linkinv(glmmTMB::fixef(Cand.mod[[29]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[29]])
broom.mixed::glance(Cand.mod[[29]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[29]])



##### Model 30 #####
# -----------------

Cand.mod[[30]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[30]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[30]])))
QQline <- qq.line(resid(Cand.mod[[30]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[30]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[30]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[30]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[30]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[30]])) # Ok
dat.resid/df.residual(Cand.mod[[30]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[30]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[30]])
family(Cand.mod[[30]])$linkinv(glmmTMB::fixef(Cand.mod[[30]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[30]])
broom.mixed::glance(Cand.mod[[30]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[30]])



##### Model 31 #####
# -----------------

Cand.mod[[31]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[31]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[31]])))
QQline <- qq.line(resid(Cand.mod[[31]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[31]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[31]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[31]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[31]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[31]])) # Ok
dat.resid/df.residual(Cand.mod[[31]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[31]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[31]])
family(Cand.mod[[31]])$linkinv(glmmTMB::fixef(Cand.mod[[31]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[31]])
broom.mixed::glance(Cand.mod[[31]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[31]])



##### Model 32 #####
# -----------------

Cand.mod[[32]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[32]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[32]])))
QQline <- qq.line(resid(Cand.mod[[32]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[32]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[32]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[32]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[32]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[32]])) # Ok
dat.resid/df.residual(Cand.mod[[32]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[32]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[32]])
family(Cand.mod[[32]])$linkinv(glmmTMB::fixef(Cand.mod[[32]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[32]])
broom.mixed::glance(Cand.mod[[32]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[32]])



##### Model 33 #####
# -----------------

Cand.mod[[33]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + plantation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[33]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[33]])))
QQline <- qq.line(resid(Cand.mod[[33]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[33]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[33]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[33]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[33]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[33]])) # Ok
dat.resid/df.residual(Cand.mod[[33]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[33]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[33]])
family(Cand.mod[[33]])$linkinv(glmmTMB::fixef(Cand.mod[[33]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[33]])
broom.mixed::glance(Cand.mod[[33]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[33]])



##### Model 34 #####
# -----------------

Cand.mod[[34]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + slope + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[34]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[34]])))
QQline <- qq.line(resid(Cand.mod[[34]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[34]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[34]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[34]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[34]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[34]])) # Ok
dat.resid/df.residual(Cand.mod[[34]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[34]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[34]])
family(Cand.mod[[34]])$linkinv(glmmTMB::fixef(Cand.mod[[34]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[34]])
broom.mixed::glance(Cand.mod[[34]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[34]])



##### Model 35 #####
# -----------------

Cand.mod[[35]] <- glmmTMB::glmmTMB(formula = efficiency~log(stand_surface) + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[35]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[35]])))
QQline <- qq.line(resid(Cand.mod[[35]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[35]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[35]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[35]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[35]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[35]])) # Ok
dat.resid/df.residual(Cand.mod[[35]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[35]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[35]])
family(Cand.mod[[35]])$linkinv(glmmTMB::fixef(Cand.mod[[35]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[35]])
broom.mixed::glance(Cand.mod[[35]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[35]])



##### Model 36 #####
# -----------------

Cand.mod[[36]] <- glmmTMB::glmmTMB(formula = efficiency~pb_fixation + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[36]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[36]])))
QQline <- qq.line(resid(Cand.mod[[36]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[36]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[36]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[36]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[36]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[36]])) # Ok
dat.resid/df.residual(Cand.mod[[36]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[36]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[36]])
family(Cand.mod[[36]])$linkinv(glmmTMB::fixef(Cand.mod[[36]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[36]])
broom.mixed::glance(Cand.mod[[36]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[36]])



##### Model 37 #####
# -----------------

Cand.mod[[37]] <- glmmTMB::glmmTMB(formula = efficiency~uprootexcav + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[37]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[37]])))
QQline <- qq.line(resid(Cand.mod[[37]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[37]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[37]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[37]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[37]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[37]])) # Ok
dat.resid/df.residual(Cand.mod[[37]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[37]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[37]])
family(Cand.mod[[37]])$linkinv(glmmTMB::fixef(Cand.mod[[37]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[37]])
broom.mixed::glance(Cand.mod[[37]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[37]])



##### Model 38 #####
# -----------------

Cand.mod[[38]] <- glmmTMB::glmmTMB(formula = efficiency~uprootexcav + plantation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[38]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[38]])))
QQline <- qq.line(resid(Cand.mod[[38]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[38]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[38]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[38]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[38]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[38]])) # Ok
dat.resid/df.residual(Cand.mod[[38]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[38]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[38]])
family(Cand.mod[[38]])$linkinv(glmmTMB::fixef(Cand.mod[[38]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[38]])
broom.mixed::glance(Cand.mod[[38]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[38]])



##### Model 39 #####
# -----------------

Cand.mod[[39]] <- glmmTMB::glmmTMB(formula = efficiency~obstacles + plantation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[39]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[39]])))
QQline <- qq.line(resid(Cand.mod[[39]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[39]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[39]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[39]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[39]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[39]])) # Ok
dat.resid/df.residual(Cand.mod[[39]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[39]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[39]])
family(Cand.mod[[39]])$linkinv(glmmTMB::fixef(Cand.mod[[39]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[39]])
broom.mixed::glance(Cand.mod[[39]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[39]])



##### Model 40 #####
# -----------------

Cand.mod[[40]] <- glmmTMB::glmmTMB(formula = efficiency~obstacles + slope + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[40]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[40]])))
QQline <- qq.line(resid(Cand.mod[[40]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[40]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[40]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[40]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[40]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[40]])) # Ok
dat.resid/df.residual(Cand.mod[[40]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[40]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[40]])
family(Cand.mod[[40]])$linkinv(glmmTMB::fixef(Cand.mod[[40]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[40]])
broom.mixed::glance(Cand.mod[[40]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[40]])



##### Model 41 #####
# -----------------

Cand.mod[[41]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + followups + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[41]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[41]])))
QQline <- qq.line(resid(Cand.mod[[41]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[41]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[41]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[41]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[41]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[41]])) # Ok
dat.resid/df.residual(Cand.mod[[41]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[41]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[41]])
family(Cand.mod[[41]])$linkinv(glmmTMB::fixef(Cand.mod[[41]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[41]])
broom.mixed::glance(Cand.mod[[41]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[41]])



##### Model 42 #####
# -----------------

Cand.mod[[42]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[42]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[42]])))
QQline <- qq.line(resid(Cand.mod[[42]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[42]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[42]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[42]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[42]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[42]])) # Ok
dat.resid/df.residual(Cand.mod[[42]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[42]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[42]])
family(Cand.mod[[42]])$linkinv(glmmTMB::fixef(Cand.mod[[42]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[42]])
broom.mixed::glance(Cand.mod[[42]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[42]])



##### Model 43 #####
# -----------------

Cand.mod[[43]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[43]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[43]])))
QQline <- qq.line(resid(Cand.mod[[43]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[43]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[43]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[43]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[43]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[43]])) # Ok
dat.resid/df.residual(Cand.mod[[43]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[43]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[43]])
family(Cand.mod[[43]])$linkinv(glmmTMB::fixef(Cand.mod[[43]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[43]])
broom.mixed::glance(Cand.mod[[43]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[43]])



##### Model 44 #####
# -----------------

Cand.mod[[44]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + plantation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[44]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[44]])))
QQline <- qq.line(resid(Cand.mod[[44]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[44]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[44]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[44]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[44]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[44]])) # Ok
dat.resid/df.residual(Cand.mod[[44]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[44]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[44]])
family(Cand.mod[[44]])$linkinv(glmmTMB::fixef(Cand.mod[[44]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[44]])
broom.mixed::glance(Cand.mod[[44]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[44]])



##### Model 45 #####
# -----------------

Cand.mod[[45]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[45]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[45]])))
QQline <- qq.line(resid(Cand.mod[[45]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[45]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[45]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[45]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[45]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[45]])) # Ok
dat.resid/df.residual(Cand.mod[[45]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[45]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[45]])
family(Cand.mod[[45]])$linkinv(glmmTMB::fixef(Cand.mod[[45]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[45]])
broom.mixed::glance(Cand.mod[[45]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[45]])



##### Model 46 #####
# -----------------

Cand.mod[[46]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + slope + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[46]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[46]])))
QQline <- qq.line(resid(Cand.mod[[46]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[46]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[46]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[46]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[46]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[46]])) # Ok
dat.resid/df.residual(Cand.mod[[46]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[46]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[46]])
family(Cand.mod[[46]])$linkinv(glmmTMB::fixef(Cand.mod[[46]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[46]])
broom.mixed::glance(Cand.mod[[46]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[46]])



##### Model 47 #####
# -----------------

Cand.mod[[47]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + sedicover_height + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[47]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[47]])))
QQline <- qq.line(resid(Cand.mod[[47]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[47]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[47]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[47]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[47]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[47]])) # Ok
dat.resid/df.residual(Cand.mod[[47]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[47]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[47]])
family(Cand.mod[[47]])$linkinv(glmmTMB::fixef(Cand.mod[[47]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[47]])
broom.mixed::glance(Cand.mod[[47]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[47]])



##### Model 48 #####
# -----------------

Cand.mod[[48]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[48]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[48]])))
QQline <- qq.line(resid(Cand.mod[[48]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[48]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[48]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[48]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[48]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[48]])) # Ok
dat.resid/df.residual(Cand.mod[[48]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[48]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[48]])
family(Cand.mod[[48]])$linkinv(glmmTMB::fixef(Cand.mod[[48]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[48]])
broom.mixed::glance(Cand.mod[[48]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[48]])



##### Model 49 #####
# -----------------

Cand.mod[[49]] <- glmmTMB::glmmTMB(formula = efficiency~distance + fully_tarped + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[49]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[49]])))
QQline <- qq.line(resid(Cand.mod[[49]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[49]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[49]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[49]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[49]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[49]])) # Ok
dat.resid/df.residual(Cand.mod[[49]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[49]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[49]])
family(Cand.mod[[49]])$linkinv(glmmTMB::fixef(Cand.mod[[49]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[49]])
broom.mixed::glance(Cand.mod[[49]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[49]])



##### Model 50 #####
# -----------------

Cand.mod[[50]] <- glmmTMB::glmmTMB(formula = efficiency~distance + log(stand_surface) + followups + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[50]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[50]])))
QQline <- qq.line(resid(Cand.mod[[50]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[50]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[50]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[50]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[50]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[50]])) # Ok
dat.resid/df.residual(Cand.mod[[50]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[50]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[50]])
family(Cand.mod[[50]])$linkinv(glmmTMB::fixef(Cand.mod[[50]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[50]])
broom.mixed::glance(Cand.mod[[50]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[50]])



##### Model 51 #####
# -----------------

Cand.mod[[51]] <- glmmTMB::glmmTMB(formula = efficiency~distance + log(stand_surface) + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[51]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[51]])))
QQline <- qq.line(resid(Cand.mod[[51]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[51]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[51]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[51]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[51]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[51]])) # Ok
dat.resid/df.residual(Cand.mod[[51]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[51]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[51]])
family(Cand.mod[[51]])$linkinv(glmmTMB::fixef(Cand.mod[[51]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[51]])
broom.mixed::glance(Cand.mod[[51]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[51]])



##### Model 52 #####
# -----------------

Cand.mod[[52]] <- glmmTMB::glmmTMB(formula = efficiency~distance + log(stand_surface) + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[52]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[52]])))
QQline <- qq.line(resid(Cand.mod[[52]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[52]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[52]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[52]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[52]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[52]])) # Ok
dat.resid/df.residual(Cand.mod[[52]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[52]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[52]])
family(Cand.mod[[52]])$linkinv(glmmTMB::fixef(Cand.mod[[52]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[52]])
broom.mixed::glance(Cand.mod[[52]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[52]])



##### Model 53 #####
# -----------------

Cand.mod[[53]] <- glmmTMB::glmmTMB(formula = efficiency~distance + log(stand_surface) + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[53]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[53]])))
QQline <- qq.line(resid(Cand.mod[[53]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[53]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[53]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[53]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[53]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[53]])) # Ok
dat.resid/df.residual(Cand.mod[[53]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[53]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[53]])
family(Cand.mod[[53]])$linkinv(glmmTMB::fixef(Cand.mod[[53]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[53]])
broom.mixed::glance(Cand.mod[[53]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[53]])



##### Model 54 #####
# -----------------

Cand.mod[[54]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + followups + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[54]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[54]])))
QQline <- qq.line(resid(Cand.mod[[54]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[54]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[54]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[54]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[54]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[54]])) # Ok
dat.resid/df.residual(Cand.mod[[54]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[54]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[54]])
family(Cand.mod[[54]])$linkinv(glmmTMB::fixef(Cand.mod[[54]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[54]])
broom.mixed::glance(Cand.mod[[54]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[54]])



##### Model 55 #####
# -----------------

Cand.mod[[55]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[55]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[55]])))
QQline <- qq.line(resid(Cand.mod[[55]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[55]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[55]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[55]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[55]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[55]])) # Ok
dat.resid/df.residual(Cand.mod[[55]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[55]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[55]])
family(Cand.mod[[55]])$linkinv(glmmTMB::fixef(Cand.mod[[55]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[55]])
broom.mixed::glance(Cand.mod[[55]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[55]])



##### Model 56 #####
# -----------------

Cand.mod[[56]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + plantation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[56]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[56]])))
QQline <- qq.line(resid(Cand.mod[[56]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[56]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[56]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[56]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[56]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[56]])) # Ok
dat.resid/df.residual(Cand.mod[[56]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[56]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[56]])
family(Cand.mod[[56]])$linkinv(glmmTMB::fixef(Cand.mod[[56]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[56]])
broom.mixed::glance(Cand.mod[[56]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[56]])



##### Model 57 #####
# -----------------

Cand.mod[[57]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[57]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[57]])))
QQline <- qq.line(resid(Cand.mod[[57]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[57]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[57]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[57]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[57]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[57]])) # Ok
dat.resid/df.residual(Cand.mod[[57]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[57]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[57]])
family(Cand.mod[[57]])$linkinv(glmmTMB::fixef(Cand.mod[[57]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[57]])
broom.mixed::glance(Cand.mod[[57]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[57]])



##### Model 58 #####
# -----------------

Cand.mod[[58]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + sedicover_height + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[58]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[58]])))
QQline <- qq.line(resid(Cand.mod[[58]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[58]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[58]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[58]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[58]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[58]])) # Ok
dat.resid/df.residual(Cand.mod[[58]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[58]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[58]])
family(Cand.mod[[58]])$linkinv(glmmTMB::fixef(Cand.mod[[58]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[58]])
broom.mixed::glance(Cand.mod[[58]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[58]])



##### Model 59 #####
# -----------------

Cand.mod[[59]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + slope + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[59]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[59]])))
QQline <- qq.line(resid(Cand.mod[[59]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[59]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[59]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[59]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[59]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[59]])) # Ok
dat.resid/df.residual(Cand.mod[[59]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[59]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[59]])
family(Cand.mod[[59]])$linkinv(glmmTMB::fixef(Cand.mod[[59]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[59]])
broom.mixed::glance(Cand.mod[[59]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[59]])



##### Model 60 #####
# -----------------

Cand.mod[[60]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + uprootexcav + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[60]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[60]])))
QQline <- qq.line(resid(Cand.mod[[60]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[60]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[60]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[60]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[60]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[60]])) # Ok
dat.resid/df.residual(Cand.mod[[60]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[60]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[60]])
family(Cand.mod[[60]])$linkinv(glmmTMB::fixef(Cand.mod[[60]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[60]])
broom.mixed::glance(Cand.mod[[60]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[60]])



##### Model 61 #####
# -----------------

Cand.mod[[61]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + log(stand_surface) + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[61]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[61]])))
QQline <- qq.line(resid(Cand.mod[[61]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[61]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[61]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[61]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[61]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[61]])) # Ok
dat.resid/df.residual(Cand.mod[[61]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[61]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[61]])
family(Cand.mod[[61]])$linkinv(glmmTMB::fixef(Cand.mod[[61]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[61]])
broom.mixed::glance(Cand.mod[[61]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[61]])



##### Model 62 #####
# -----------------

Cand.mod[[62]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + followups + geomem + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[62]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[62]])))
QQline <- qq.line(resid(Cand.mod[[62]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[62]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[62]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[62]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[62]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[62]])) # Ok
dat.resid/df.residual(Cand.mod[[62]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[62]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[62]])
family(Cand.mod[[62]])$linkinv(glmmTMB::fixef(Cand.mod[[62]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[62]])
broom.mixed::glance(Cand.mod[[62]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[62]])



##### Model 63 #####
# -----------------

Cand.mod[[63]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + followups + obstacles + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[63]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[63]])))
QQline <- qq.line(resid(Cand.mod[[63]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[63]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[63]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[63]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[63]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[63]])) # Ok
dat.resid/df.residual(Cand.mod[[63]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[63]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[63]])
family(Cand.mod[[63]])$linkinv(glmmTMB::fixef(Cand.mod[[63]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[63]])
broom.mixed::glance(Cand.mod[[63]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[63]])



##### Model 64 #####
# -----------------

Cand.mod[[64]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + followups + tarping_duration + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[64]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[64]])))
QQline <- qq.line(resid(Cand.mod[[64]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[64]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[64]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  plot(e, main = i, las = 1)
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[64]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[64]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[64]])) # Ok
dat.resid/df.residual(Cand.mod[[64]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[64]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[64]])
family(Cand.mod[[64]])$linkinv(glmmTMB::fixef(Cand.mod[[64]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[64]])
broom.mixed::glance(Cand.mod[[64]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[64]])



##### Model 65 #####
# -----------------

Cand.mod[[65]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + pb_fixation + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[65]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[65]])))
QQline <- qq.line(resid(Cand.mod[[65]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[65]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[65]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[65]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[65]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[65]])) # Ok
dat.resid/df.residual(Cand.mod[[65]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[65]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[65]])
family(Cand.mod[[65]])$linkinv(glmmTMB::fixef(Cand.mod[[65]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[65]])
broom.mixed::glance(Cand.mod[[65]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[65]])



##### Model 66 #####
# -----------------

Cand.mod[[66]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + slope + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[66]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[66]])))
QQline <- qq.line(resid(Cand.mod[[66]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[66]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[66]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[66]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[66]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[66]])) # Ok
dat.resid/df.residual(Cand.mod[[66]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[66]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[66]])
family(Cand.mod[[66]])$linkinv(glmmTMB::fixef(Cand.mod[[66]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[66]])
broom.mixed::glance(Cand.mod[[66]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[66]])



##### Model 67 #####
# -----------------

Cand.mod[[67]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + uprootexcav + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[67]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[67]])))
QQline <- qq.line(resid(Cand.mod[[67]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[67]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[67]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[67]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[67]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[67]])) # Ok
dat.resid/df.residual(Cand.mod[[67]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[67]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[67]])
family(Cand.mod[[67]])$linkinv(glmmTMB::fixef(Cand.mod[[67]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[67]])
broom.mixed::glance(Cand.mod[[67]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[67]])



##### Model 68 #####
# -----------------

Cand.mod[[68]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + distance + followups + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[68]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[68]])))
QQline <- qq.line(resid(Cand.mod[[68]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[68]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[68]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[68]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[68]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[68]])) # Ok
dat.resid/df.residual(Cand.mod[[68]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[68]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[68]])
family(Cand.mod[[68]])$linkinv(glmmTMB::fixef(Cand.mod[[68]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[68]])
broom.mixed::glance(Cand.mod[[68]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[68]])



##### Model 69 #####
# -----------------

Cand.mod[[69]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + distance + log(stand_surface) + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[69]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[69]])))
QQline <- qq.line(resid(Cand.mod[[69]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[69]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[69]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[69]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[69]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[69]])) # Ok
dat.resid/df.residual(Cand.mod[[69]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[69]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[69]])
family(Cand.mod[[69]])$linkinv(glmmTMB::fixef(Cand.mod[[69]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[69]])
broom.mixed::glance(Cand.mod[[69]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[69]])



##### Model 70 #####
# -----------------

Cand.mod[[70]] <- glmmTMB::glmmTMB(formula = efficiency~fully_tarped + obstacles + distance * uprootexcav + (1|manager_id), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation:
ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[70]],
                                                                              type = "pearson"), x = fitted(Cand.mod[[70]])))
QQline <- qq.line(resid(Cand.mod[[70]], type = "pearson"))
ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[70]], type = "pearson"))) +
  ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_15):
dat.sim <- simulate(Cand.mod[[70]], n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
  e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
  resid[i] <- e(eff$efficiency[i] + runif(250, -0.5, 0.5))
}
par(.pardefault) # To restore defaults graphical parameters
plot(resid ~ fitted(Cand.mod[[70]]))

### Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(Cand.mod[[70]], type = "pearson")^2)
1 - pchisq(dat.resid, df.residual(Cand.mod[[70]])) # Ok
dat.resid/df.residual(Cand.mod[[70]]) # Ok
# To plot the observed vs. fitted values:
plot(x = fitted(Cand.mod[[70]]), y = eff$efficiency, xlab = "Fitted values", ylab = "Observed values")

### Exploring the model parameters and test hypotheses:
summary(Cand.mod[[70]])
family(Cand.mod[[70]])$linkinv(glmmTMB::fixef(Cand.mod[[70]])$cond) # To get interpretable coefficients.
broom.mixed::tidy(Cand.mod[[70]])
broom.mixed::glance(Cand.mod[[70]])

### Computing a Pseudo-R2:
pseudo.R2glmm(model = Cand.mod[[70]])



