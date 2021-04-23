eff <- readr::read_csv(here::here("mydata", "erad.csv"), col_names = TRUE, col_types =
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
                           pb_durability = readr::col_factor(c("0", "1"))))


## Dans gamlss, je dois utiliser re() pour fitter un random effect (s'utilise globalement comme lme() - voir
# ressources en ligne, dont l'aide dans mes favoris),
# car random() convient lorsque Y est normal, (enfin je crois).
## Mes modèles mixtes doivent être ajustés avec du ML et pas du REML (sinon, ils ne sont pas comparables) !


Model <- (1:17)

Candidate <- c('Manag + Year',
               'Manag + Year + Elevation', 'Manag + Year + Distance', 'Manag + Year + Prop_Sedim', 'Manag + Year + Depth_Sedim', "Manag + Year + Dura_Inund", "Manag + Year + Inte_Inund", "Manag + Year + Mean_Temp", "Manag + Year + Sum_Prec",
               'Manag + Year * Elevation', 'Manag + Year * Distance', 'Manag + Year * Prop_Sedim', 'Manag + Year * Depth_Sedim', "Manag + Year * Dura_Inund", "Manag + Year * Inte_Inund", "Manag + Year * Mean_Temp", "Manag + Year * Sum_Prec")

Cand.model <- data.frame(Model, Candidate)

AICc <- model.sel(object=Cand.mod)
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


Parameters <- c('Manag', 'Year', 'Elevation', 'Distance', 'Prop_Sedim', 'Depth_Sedim', 'Dura_Inund', 'Inte_Inund', 'Mean_Temp', 'Sum_Prec',
                'Year * Elevation', 'Year * Distance', 'Year * Prop_Sedim', 'Year * Depth_Sedim', 'Year * Dura_Inund', 'Year * Inte_Inund', 'Year * Mean_Temp', 'Year * Sum_Prec')
Var <- c("managplowing", "annee19", "topographie.2", "distance_eau.2", "sub_limons.2", "h_alluvions.2", "Dura_inund.2", "Inte_inund.2", "Mean_Temp.2", "Sum_Prec.2",
         "annee19:topographie.2", "annee19:distance_eau.2", "annee19:sub_limons.2", "annee19:h_alluvions.2", "annee19ura_inund.2", "annee19:Inte_inund.2", "annee19:Mean_Temp.2", "annee19:Sum_Prec.2")
Var.Imp <- c("manag", "annee", "topographie.2", "distance_eau.2", "sub_limons.2", "h_alluvions.2", "Dura_inund.2", "Inte_inund.2", "Mean_Temp.2", "Sum_Prec.2",
             "annee:topographie.2", "annee:distance_eau.2", "annee:sub_limons.2", "annee:h_alluvions.2", "anneeura_inund.2", "annee:Inte_inund.2", "annee:Mean_Temp.2", "annee:Sum_Prec.2")

Para.model <- data.frame(Parameters, Var, Var.Imp)
top.models <- get.models(AICc, cumsum(weight)<=1) # Pour prendre l'ensemble des modèles a priori (sinon 0.95)
# top.models <- get.models(AICc, delta<=7)
Parameter <- model.avg(top.models, revised.var=T, adjusted=T, fit=T) # Model averaging
Parameter.model <- as.data.frame(cbind(coefTable(Parameter), confint(Parameter)))
Parameter.model <- Parameter.model[rownames(Parameter.model) != "X.Intercept.", ]
Parameter.model$Var <- row.names(Parameter.model)
Parameter.model <- merge(Parameter.model, Para.model, by="Var", all=TRUE)
Imp <- as.data.frame(format(round(importance(Parameter), digits=2)))
colnames(Imp) <- "Imp."
Imp$Var.Imp <- row.names(Imp)
Parameter.model <- merge(Parameter.model, Imp, by="Var.Imp", all=TRUE)
Parameter.model$'Estimate (±SE)' <- paste0(format(round(Parameter.model$Estimate, digits=3), trim=T), " (±", format(round(Parameter.model$'Std. Error', digits=3), trim=T), ")")
Parameter.model$'(95% CI)' <- paste0("(", format(round(Parameter.model$'2.5 %', digits=3), trim=T), "; ", format(round(Parameter.model$'97.5 %', digits=3), trim=T), ")")
Parameter.model$Variables <- "FD_sla"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Variables", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)", 'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Parameters), ]
write.table(Parameter.model,"output/Model/Parameters_FD_sla.csv", sep=";", row.names=F, col.names=T)
