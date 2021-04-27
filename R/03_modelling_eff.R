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



##### Data preparation for modelling 'efficiency' #####
# --------------------------------------------------- #

# List of used packages (for publication or package building): here, readr, gamlss, gamlss.dist, MuMIn,
# emmeans,

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
  dplyr::mutate(efficiency = efficiency/10) -> eff # For 'efficiency' to be look like a percentage that could
# be modelled using a Beta Regression Model.
summary(eff)


##### Test model #####
# --------------------
#Cand.mod <- list()

#Cand.mod[[1]] <- gamlss::gamlss(efficiency~distance, data = eff, family = BEOI())
#mod1 <- gamlss::gamlss(efficiency~distance, data = eff, family = gamlss.dist::BEOI())
library(gamlss)
mod1 <- gamlss::gamlss(formula = efficiency~distance + re(random = ~1|manager_id, method = "ML"),
                       nu.formula = ~distance + re(random = ~1|manager_id, method = "ML"),
                       data = eff, family = gamlss.dist::BEOI())
#mod2 <- gamlss::gamlss(efficiency~distance + gamlss::random(x = manager_id),
#                      data = eff, family = gamlss.dist::BEOI())
# Checking the residuals:
plot(mod1)
#plot(mod2)
plot(fitted(mod1)~eff$efficiency)

# Checking the goodness-of-fit and overdispersion:
dat.resid <- sum(resid(mod1, type = "weighted", what = "mu")^2)
1 - pchisq(dat.resid, mod1$df.resid) # Conclusions: Pearson residuals do indicate a lack of fit (p values
# far lower than 0.05)
dat.resid/df.residual(mod1) # Conclusions: the data are definitely not overdispersed.
# Exploring parameters and hypotheses:
summary(mod1)
binomial()$linkinv(coef(mod1))  # logit inverse: to get interpretable parameter estimates
# Further trend exploration:
gamlss::Rsq(mod1) # R^2 seems very good!


## Dans gamlss, je dois utiliser re() pour fitter un random effect (s'utilise globalement comme lme() - voir
# ressources en ligne, dont l'aide dans mes favoris), car random() convient lorsque Y est normal,
# (enfin je crois).
## Mes modèles mixtes doivent être ajustés avec du ML et pas du REML (sinon, ils ne sont pas comparables) !
# Mais s'il me reste qu'un "best model", je peux le refitter en REML car ses estimations de SE seront plus
# précises.

# Donc je CENTRE mes variables impliquées dans des interactions. Mais si à la fin, je vire les modèles à
# interaction, alors pas besoin, non ? ASK CV§§§§§§§§§





####### Code de Philippe 2021
Model <- (1:17)

Candidate <- c('Manag + Year',
               'Manag + Year + Elevation', 'Manag + Year + Distance', 'Manag + Year + Prop_Sedim',
               'Manag + Year + Depth_Sedim', "Manag + Year + Dura_Inund", "Manag + Year + Inte_Inund",
               "Manag + Year + Mean_Temp", "Manag + Year + Sum_Prec",
               'Manag + Year * Elevation', 'Manag + Year * Distance', 'Manag + Year * Prop_Sedim',
               'Manag + Year * Depth_Sedim', "Manag + Year * Dura_Inund", "Manag + Year * Inte_Inund",
               "Manag + Year * Mean_Temp", "Manag + Year * Sum_Prec")

Cand.model <- data.frame(Model, Candidate)

AICc <- MuMIn::model.sel(object=Cand.mod) # Rank models (but where does Cand.mod come from?)
AICc.model <- as.data.frame(AICc)
AICc.model$csweigth <- cumsum(AICc.model$weight)
AICc.model$Model <- row.names(AICc.model)
AICc.model <- merge(AICc.model, Cand.model, by="Model")
AICc.model <- merge(AICc.model, R.ajust, by="Model")
AICc.model <- AICc.model[order(AICc.model$delta),]
AICc.model$Variables <- "FD_sla"
AICc.model$Rank <- 1:nrow(AICc.model)
AICc.model <- AICc.model[,c("Rank", "Model", "Variables", "Candidate", "df", "AICc", "delta", "weight", "R2")]
AICc.model$AICc <- format(round(AICc.model$AICc, digits=1))
AICc.model$delta <- format(round(AICc.model$delta, digits=3))
AICc.model$weight <- format(round(AICc.model$weight, digits=3))
AICc.model$R2 <- format(round(AICc.model$R2, digits=3))
colnames(AICc.model)[colnames(AICc.model) == 'Candidate'] <- 'Candidate model'
colnames(AICc.model)[colnames(AICc.model) == 'df'] <- 'k'
colnames(AICc.model)[colnames(AICc.model) == 'delta'] <- 'delta AICc'
colnames(AICc.model)[colnames(AICc.model) == 'weight'] <- 'W'
write.table(AICc.model,"output/Model/Models_FD_sla.csv", sep=";", row.names=F, col.names=T)


Parameters <- c('Manag', 'Year', 'Elevation', 'Distance', 'Prop_Sedim', 'Depth_Sedim', 'Dura_Inund',
                'Inte_Inund', 'Mean_Temp', 'Sum_Prec',
                'Year * Elevation', 'Year * Distance', 'Year * Prop_Sedim', 'Year * Depth_Sedim',
                'Year * Dura_Inund', 'Year * Inte_Inund', 'Year * Mean_Temp', 'Year * Sum_Prec')
Var <- c("managplowing", "annee19", "topographie.2", "distance_eau.2", "sub_limons.2", "h_alluvions.2",
         "Dura_inund.2", "Inte_inund.2", "Mean_Temp.2", "Sum_Prec.2",
         "annee19:topographie.2", "annee19:distance_eau.2", "annee19:sub_limons.2", "annee19:h_alluvions.2",
         "annee19ura_inund.2", "annee19:Inte_inund.2", "annee19:Mean_Temp.2", "annee19:Sum_Prec.2")
Var.Imp <- c("manag", "annee", "topographie.2", "distance_eau.2", "sub_limons.2", "h_alluvions.2",
             "Dura_inund.2", "Inte_inund.2", "Mean_Temp.2", "Sum_Prec.2",
             "annee:topographie.2", "annee:distance_eau.2", "annee:sub_limons.2", "annee:h_alluvions.2",
             "anneeura_inund.2", "annee:Inte_inund.2", "annee:Mean_Temp.2", "annee:Sum_Prec.2")

Para.model <- data.frame(Parameters, Var, Var.Imp)
top.models <- MuMIn::get.models(AICc, cumsum(weight)<=1) # Pour prendre l'ensemble des modèles a priori
# (sinon 0.95)
# top.models <- get.models(AICc, delta<=7)
Parameter <- MuMIn::model.avg(top.models, revised.var=T, adjusted=T, fit=T) # Model averaging
Parameter.model <- as.data.frame(cbind(coefTable(Parameter), confint(Parameter)))
Parameter.model <- Parameter.model[rownames(Parameter.model) != "X.Intercept.", ]
Parameter.model$Var <- row.names(Parameter.model)
Parameter.model <- merge(Parameter.model, Para.model, by="Var", all=TRUE)
Imp <- as.data.frame(format(round(importance(Parameter), digits=2)))
colnames(Imp) <- "Imp."
Imp$Var.Imp <- row.names(Imp)
Parameter.model <- merge(Parameter.model, Imp, by="Var.Imp", all=TRUE)
Parameter.model$'Estimate (±SE)' <- paste0(format(round(Parameter.model$Estimate, digits=3), trim=T),
                                           " (±", format(round(Parameter.model$'Std. Error', digits=3),
                                                         trim=T), ")")
Parameter.model$'(95% CI)' <- paste0("(", format(round(Parameter.model$'2.5 %', digits=3), trim=T),
                                     "; ", format(round(Parameter.model$'97.5 %', digits=3), trim=T), ")")
Parameter.model$Variables <- "FD_sla"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Variables", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)",
                                      'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Parameters), ]
write.table(Parameter.model,"output/Model/Parameters_FD_sla.csv", sep=";", row.names=F, col.names=T)








######## Code from https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples:
# download the file directly from the website
dater<-read.csv("https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples/Example%20data.csv?attredirects=0&d=1")
head(dater)
summary(dater)

require(MuMIn)
#First, fit 4 candidate linear models to explain variation in density
mod1<-lm(density~distance+elev, data = dater)
mod2<-lm(density~slope+pct.cover, data = dater)
mod3<-lm(density~slope+distance, data = dater)
mod4<-lm(density~slope+distance+elev, data = dater)
# use the mod.sel function to conduct model selection
# and put output into object out.put
out.put<-model.sel(mod1,mod2,mod3,mod4)
out.put

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!! (weights add to 1)
subset(out.put, delta <5)
# select models 95% cumulative weight criteria
# IMPORTANT: Weights have been renormalized!!
subset(out.put, cumsum(out.put$weight) <= .95)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want
sel.table<-as.data.frame(out.put)[6:11]
# a little clean-up, lets round things a bit
sel.table[,3:4]<- round(sel.table[,3:4],2)
sel.table[,5:6]<- round(sel.table[,5:6],3)
sel.table
# that’s better

# how about a little renaming columns to fit proper conventions
# number of parameters (df) should be K
names(sel.table)[1] = "K"
## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
# replace Model name with formulas little tricky so be careful
for(i in 1:nrow(sel.table)) sel.table$Model[i]<- as.character(formula(paste(sel.table$Model[i])))[3]
# let's see what is in there
sel.table

# Importance weights for individual predictor variables
# calculated using the importance function
importance(out.put)
# The number of candidate models that a parameter occurs can have a big effect of the importance weight.
# For example, the intercept is included in all models, so the importance weight is 1 (hence it is never
# shown). In the above output, pct.cover is in only one model so interpret weights with caution.

# There are two methods for model-averaging defined by Burnham and Anderson as , where parameter estimates
# are averaged over all models in which predictor xj occurs and  where parameter estimates are averaged over
# all models not just those in which predictor xj occurs. MuMIn function model.avg conducts both types of
# model averaging and reports the first type of model averaging as “subset” and the second type as “full.
# Model average using all candidate models, always use revised.var = TRUE
MA.ests<-model.avg(out.put, revised.var = TRUE)
MA.ests
# Le reste du code ne marche pas (le package MuMin a changé depuis), donc je m'arrête-là.
