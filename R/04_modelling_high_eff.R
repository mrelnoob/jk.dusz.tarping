# ******************************************************************************************** #
# ******************************************************************************************** #
# ************************************ MODELLING high_eff ************************************ #
# ******************************************************************************************** #
# ******************************************************************************************** #

### IMPORTANT NOTE (only when in 'package development' mod): to avoid creating notes for unquoted variables,
# I must add the following code at the beginning of every source file (e.g. R/myscript.R) that uses
# unquoted variables, so in front of all my scripts doing any kind of analyses (otherwise, I should
# always assign each variable, e.g. mydata$myvariable, which is quite time consuming and wearisome).
# if(getRversion() >= "2.29.1")  utils::globalVariables(c(
#   "manager_id", "xp_id", "country", "efficiency", "eff_eradication", "high_eff",
#   "latitude", "elevation", "goals", "slope", "coarse_env", "obstacles", "flood",
#   "geomem", "maxveg", "uprootexcav", "stand_surface", "fully_tarped", "distance", "tarping_duration",
#   "stripsoverlap_ok", "tarpfix_pierced", "tarpfix_multimethod", "sedicover_height",
#   "trench_depth", "plantation", "followups",
#   "pb_fixation", "pb_durability"))
# NOTE: This variable list has been copy-pasted from another R script and should thus be updated to
# contain the actual list of variables used in this script!
# -------------------------------------------------------------------------------------------------------- #





# ================================================= #
##### Data preparation for modelling 'high_eff' #####
# ================================================= #

# List of used packages (for publication or package building): here, readr, MuMIn, lme4, ggplot2,
# (broom.mixed), stats

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
                    pb_durability = readr::col_factor(c("0", "1")))) -> eff
summary(eff)



##### Final pre-modelling assumption checks ##### (run only when required)
# --------------------------------------------- #

# ### Testing the relevance of the random effect structure:
# m0.glm <- stats::glm(high_eff ~ fully_tarped, data = eff, family = binomial)
# m0.glmer <- lme4::glmer(high_eff ~ fully_tarped + (1|manager_id), data = eff, family = binomial)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glmer <- AIC(logLik(m0.glmer))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
# pchisq(as.numeric(null.id), df=1, lower.tail=F)
# rm(m0.glm, m0.glmer, aic.glm, aic.glmer, null.id)
# # The Likelihood Ratio Tests are NOT significant so the use of the random effect structure may not be
# # necessary! However, further tests on the model residuals may indicate otherwise.


# IMPORTANT NOTE:
# Since we'll use a regularization modelling method (i.e. penalized regression) to avoid overfitting
# and deal with our "complete separation" problem, we will not use model selection and multimodel inference
# for this response variable.
# Consequently, we will assess further logistic regression assumptions using the most parsimonious full model
# we could think of, i.e. including all the available predictors we thought should be important to explain
# our outcome (i.e. the eradication or near-eradication of knotweeds after tarping) based on our knowledge.



# ### (Re-)Assessing the linearity assumption:
# eff2 <- eff[,c("high_eff", "add_control", "distance", "elevation", "fully_tarped", "geomem", "obstacles",
#                "plantation", "pb_fixation", "pb_durability", "repairs", "slope", "stand_surface",
#                "stripsoverlap_ok", "sedicover_height", "tarping_duration", "uprootexcav")]
# model <- glm(high_eff ~., data = eff2, family = binomial)
# # Predict the probability (p) of high efficacy:
# probabilities <- predict(model, type = "response")
# # Transforming predictors
# mydata <- eff2 %>%
#   dplyr::select_if(is.numeric) %>%
#   dplyr::mutate("stand_surface (log)" = log(stand_surface)) %>%
#   dplyr::mutate("distance (log+1)" = log(distance+1)) # These transformations were made to linearize the relationships
# predictors <- colnames(mydata)
# # Bind the logit and tidying the data for plot (ggplot2, so long format)
# mydata <- mydata %>%
#   dplyr::mutate(logit = log(probabilities/(1-probabilities))) %>%
#   tidyr::gather(key = "predictors", value = "predictor.value", -logit)
# # Create scatterplot
# ggplot2::ggplot(mydata, ggplot2::aes(y = logit, x = predictor.value))+
#   ggplot2::geom_point(size = 0.5, alpha = 0.5) +
#   ggplot2::geom_smooth(method = "loess") +
#   ggplot2::theme_bw() +
#   ggplot2::facet_wrap(~predictors, scales = "free_x")
# # It appears that including binary variables in the model impedes the linearity assessment as logit values
# # become bimodal. However, penalized regression methods such as Ridge or LASSO do not explicitly assume
# # the same assumptions as un-penalized models. Moreover, all variables are standardized before performing
# # Ridge or LASSO regressions and it would be quite odd to transform (e.g. log) prior to standardize them and
# # it would surely highly complicate interpretations. Therefore, we will fit the model with non-transformed
# # variables.



# ### Assessing multicollinearity:
# car::vif(mod = model) # There is no signs of strong multicollinearity as all GVIF value are under 2.1
# rm(eff2, model, probabilities, predictors, mydata)





# =========================================================================== #
##### Cross-validated Ridge logistic regression with bootstrap validation #####
# =========================================================================== #


1) Faire un autre modèle ridge TOUTES les interactions (voir Favoris SO).
2) Implémenter la procédure d'enhanced bootstrap' (bootstrapping optimism) avec les deux types de modèles
ridge, i.e. avec et sans interactions (voir Favoris RStudio & Freerangestats).
3) Comparer les résultats !

4) Refaire tourner ça, mais en comparant si la standardisation manuelle marche mieux en utilisant le meilleur
modèle (voir Favoris CV).
5) Faire des prédictions avec des valeurs moyennes pour tous les prédicteurs et des plages de valeurs
pour distance et stand_surface!
6) Plotter les résultats.
7) Enjoy!


##### Comparing Ridge model with and without interactions #####
# ----------------------------------------------------------- #

### Data preparation:

mydata <- eff[,c("high_eff", "add_control", "distance", "fully_tarped", "geomem", "obstacles",
                 "plantation", "pb_fixation", "pb_durability", "repairs", "slope", "stand_surface",
                 "stripsoverlap_ok", "sedicover_height", "tarping_duration", "uprootexcav")]
# Create data matrices (as accepted by glmnet):
x <- stats::model.matrix(high_eff~., mydata)[,-1] # Matrix of potential predictors (WITHOUT interactions)
f <- as.formula(high_eff~add_control+distance+fully_tarped+geomem+obstacles+plantation+pb_fixation
                +pb_durability+repairs+slope+stand_surface+stripsoverlap_ok+sedicover_height+tarping_duration
                +uprootexcav+distance*stand_surface+pb_fixation*repairs+obstacles*add_control+0)
x_int <- stats::model.matrix(f, mydata)[,-1] # Matrix of potential predictors (WITH interactions)
y <- jk.dusz.tarping::as.numfactor(x = mydata$high_eff) %>% as.matrix() # Same for Y



### Ridge regression:
# Find the optimal value of lambda that minimizes the cross-validation error:
set.seed(653)
cv.ridge <- glmnet::cv.glmnet(x = x, y = y, alpha = 0, family = "binomial",
                              type.measure = "deviance", nfolds = 10)
cv.ridge_int <- glmnet::cv.glmnet(x = x_int, y = y, alpha = 0, family = "binomial",
                              type.measure = "deviance", nfolds = 10)
# plot(cv.ridge)
# plot(cv.ridge_int)
# glmnet::coef.glmnet(object = cv.ridge, s = cv.ridge$lambda.min) # To have look at the coefficients
# glmnet::coef.glmnet(object = cv.ridge_int, s = cv.ridge$lambda.min) # To have look at the coefficients

# Compute the full Ridge model:
ridge.model <- glmnet::glmnet(x = x, y = y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min)
ridge.model_int <- glmnet::glmnet(x = x_int, y = y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min)



### Evaluate the performance of the full Ridge model (WITHOUT interactions):
pred.class <- predict(object = ridge.model, newx = x,
                         s = cv.ridge$lambda.min, type = "class") %>% as.factor()
pred.prob  <- predict(ridge.model, newx = x, s = cv.ridge$lambda.min, type = "response")

MLmetrics::LogLoss(y_pred = pred.prob, y_true = y) # 0.472
MLmetrics::AUC(y_pred = pred.prob, y_true = y) # 0.859
MLmetrics::R2_Score(y_pred = pred.prob, y_true = y) # 0.316
stats::deviance(ridge.model) # 80.268

conf.mat <- MLmetrics::ConfusionMatrix(y_pred = pred.class, y_true = y) # Or table(y, pred.class)
error_rate <- (conf.mat[2]+conf.mat[3])/length(y)
error_rate # Mean error rate = 0.223
tpr <- (conf.mat[4])/(conf.mat[2]+conf.mat[4])
tpr # True positive rate = 0.567 (not so good)

# # To plot the ROC curve:
# pred <- ROCR::prediction(predictions = pred.prob, labels = mydata$high_eff)
# auc <- ROCR::performance(prediction.obj = pred, measure = "auc")@y.values[[1]][1]
# perf <- ROCR::performance(prediction.obj = pred, measure = "tpr","fpr")
# plot(perf, col="navyblue", cex.main=1,
#      main= paste("Logistic Regression ROC Curve: AUC =", round(auc,3)))
# abline(a=0, b = 1, col='darkorange1')



### Evaluate the performance of the full Ridge model (WITH interactions):
pred.class_int <- predict(object = ridge.model_int, newx = x_int,
                         s = cv.ridge_int$lambda.min, type = "class") %>% as.factor()
pred.prob_int  <- predict(ridge.model_int, newx = x_int, s = cv.ridge$lambda.min, type = "response")

MLmetrics::LogLoss(y_pred = pred.prob_int, y_true = y) # 0.469
MLmetrics::AUC(y_pred = pred.prob_int, y_true = y) # 0.86
MLmetrics::R2_Score(y_pred = pred.prob_int, y_true = y) # 0.321
stats::deviance(ridge.model_int) # 79.678

conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = pred.class_int, y_true = y) # Or table(y, pred.class_int)
error_rate_int <- (conf.mat_int[2]+conf.mat_int[3])/length(y)
error_rate_int # Mean error rate = 0.212
tpr_int <- (conf.mat_int[4])/(conf.mat_int[2]+conf.mat_int[4])
tpr_int # True positive rate = 0.567 (not so good)

# # To plot the ROC curve:
# pred <- ROCR::prediction(predictions = pred.prob_int, labels = mydata$high_eff)
# auc <- ROCR::performance(prediction.obj = pred, measure = "auc")@y.values[[1]][1]
# perf <- ROCR::performance(prediction.obj = pred, measure = "tpr","fpr")
# plot(perf, col="navyblue", cex.main=1,
#      main= paste("Logistic Regression ROC Curve: AUC =", round(auc,3)))
# abline(a=0, b = 1, col='darkorange1')



### Enhanced (optimism) bootstrap comparison:
# Create a function suitable for boot that will return the optimism estimates for statistics testing
# models against the full original sample:
compare_opt <- function(orig_data, i){
  # Create the resampled data
  train_data <- orig_data[i, ]

  x.train <- stats::model.matrix(high_eff~., train_data)[,-1] # Matrix of potential predictors (WITHOUT interactions)
  f <- as.formula(high_eff~add_control+distance+fully_tarped+geomem+obstacles+plantation+pb_fixation
                  +pb_durability+repairs+slope+stand_surface+stripsoverlap_ok+sedicover_height+tarping_duration
                  +uprootexcav+distance*stand_surface+pb_fixation*repairs+obstacles*add_control+0)
  x.train_int <- stats::model.matrix(f, train_data)[,-1] # Matrix of potential predictors (WITH interactions)
  y.train <- jk.dusz.tarping::as.numfactor(x = train_data$high_eff) %>% as.matrix() # Same for Y

  # Run the entire modelling process:
  cv.lambda <- glmnet::cv.glmnet(x = x.train, y = y.train, alpha = 0, family = "binomial",
                                type.measure = "deviance", nfolds = 10)
  model_full <- glmnet::glmnet(x = x.train, y = y.train, alpha = 0, family = "binomial", lambda = cv.lambda$lambda.min)
  cv.lambda_int <- glmnet::cv.glmnet(x = x.train_int, y = y.train, alpha = 0, family = "binomial",
                                type.measure = "deviance", nfolds = 10)
  model_full_int <- glmnet::glmnet(x = x.train_int, y = y.train, alpha = 0, family = "binomial", lambda = cv.lambda_int$lambda.min)

  # Predict the values on the trained, resampled data
  train_pred.class <- predict(object = model_full, newx = x.train,
                              s = cv.lambda$lambda.min, type = "class") %>% as.factor()
  train_pred.prob  <- predict(object = model_full, newx = x.train,
                              s = cv.lambda$lambda.min, type = "response")
  train_pred.class_int <- predict(object = model_full_int, newx = x.train_int,
                              s = cv.lambda_int$lambda.min, type = "class") %>% as.factor()
  train_pred.prob_int  <- predict(object = model_full_int, newx = x.train_int,
                              s = cv.lambda_int$lambda.min, type = "response")

  # Predict the values on the original, unresampled data
  full_pred.class <- predict(object = model_full, newx = x,
                             s = cv.lambda$lambda.min, type = "class") %>% as.factor()
  full_pred.prob  <- predict(object = model_full, newx = x,
                             s = cv.lambda$lambda.min, type = "response")
  full_pred.class_int <- predict(object = model_full_int, newx = x_int,
                             s = cv.lambda_int$lambda.min, type = "class") %>% as.factor()
  full_pred.prob_int  <- predict(object = model_full_int, newx = x_int,
                             s = cv.lambda_int$lambda.min, type = "response")

  train_conf.mat <- MLmetrics::ConfusionMatrix(y_pred = train_pred.class, y_true = y.train)
  train_error_rate <- (train_conf.mat[2]+train_conf.mat[3])/length(y.train)
  train_tpr <- (train_conf.mat[4])/(train_conf.mat[2]+train_conf.mat[4])
  full_conf.mat <- MLmetrics::ConfusionMatrix(y_pred = full_pred.class, y_true = y)
  full_error_rate <- (full_conf.mat[2]+full_conf.mat[3])/length(y)
  full_tpr <- (full_conf.mat[4])/(full_conf.mat[2]+full_conf.mat[4])

  train_conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = train_pred.class_int, y_true = y.train)
  train_error_rate_int <- (train_conf.mat_int[2]+train_conf.mat_int[3])/length(y.train)
  train_tpr_int <- (train_conf.mat_int[4])/(train_conf.mat_int[2]+train_conf.mat_int[4])
  full_conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = full_pred.class_int, y_true = y)
  full_error_rate_int <- (full_conf.mat_int[2]+full_conf.mat_int[3])/length(y)
  full_tpr_int <- (full_conf.mat_int[4])/(full_conf.mat_int[2]+full_conf.mat_int[4])


  # Return a vector of summary optimism results
  results <- c(
    boot_LogLoss = MLmetrics::LogLoss(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::LogLoss(y_pred = full_pred.prob, y_true = y),
    boot_AUC = MLmetrics::AUC(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::AUC(y_pred = full_pred.prob, y_true = y),
    boot_R2 = MLmetrics::R2_Score(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::R2_Score(y_pred = full_pred.prob, y_true = y),
    boot_m.error = train_error_rate - full_error_rate,
    boot_tpr = train_tpr - full_tpr,
    boot_LogLoss_int = MLmetrics::LogLoss(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::LogLoss(y_pred = full_pred.prob_int, y_true = y),
    boot_AUC_int = MLmetrics::AUC(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::AUC(y_pred = full_pred.prob_int, y_true = y),
    boot_R2_int = MLmetrics::R2_Score(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::R2_Score(y_pred = full_pred.prob_int, y_true = y),
    boot_m.error_int = train_error_rate_int - full_error_rate_int,
    boot_tpr_int = train_tpr_int - full_tpr_int
  )
  return(results)
}



### Perform bootstrapping and return optimism-corrected results:
res_opt <- boot::boot(data = mydata, statistic = compare_opt, R = 1000) # For large datasets, use a parallel
# cluster computing (it's quite straightforward)

# Calculate the results for the original model:
original <- c(
  MLmetrics::LogLoss(y_pred = pred.prob, y_true = y),
  MLmetrics::AUC(y_pred = pred.prob, y_true = y),
  MLmetrics::R2_Score(y_pred = pred.prob, y_true = y),
  error_rate,
  tpr,
  MLmetrics::LogLoss(y_pred = pred.prob_int, y_true = y),
  MLmetrics::AUC(y_pred = pred.prob_int, y_true = y),
  MLmetrics::R2_Score(y_pred = pred.prob_int, y_true = y),
  error_rate_int,
  tpr_int
)

# Compute the mean bootstrapped results and the enhanced results:
optimism <- apply(na.omit(res_opt$t), 2, mean)
corrected_results <- original - optimism %>% as.data.frame(row.names = c("LogLoss", "AUC",
                                                          "R2", "Mean Error Rate", "True Positive Rate",
                                                          "LogLoss (w/ inter.)", "AUC (w/ inter.)",
                                                          "R2 (w/ inter.)", "Mean Error Rate (w/ inter.)",
                                                          "True Positive Rate (w/ inter.)"))
print(corrected_results <- dplyr::rename(.data = corrected_results, "Optimism-corrected value" = .))


### NUL!











# ############## Elastic Net avec caret tuning ! #################
# ################################################################
# # To prepare for a 10-fold cross validation with 5 repeats:
# train.control <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 5)
#
# # To convert the data to a matrix format accepted by glmnet:
# x <- stats::model.matrix(high_eff~., mydata)[,-1] # Creates a matrix of the potential predictors
# y <- jk.dusz.tarping::as.numfactor(x = mydata$high_eff) %>% as.matrix() # Same for Y
#
#
# ### Train a model using the repeated 10-fold cross validation to find the best hyperparameters:
# set.seed(21)
# trained.mod <- caret::train(high_eff~., data = mydata, method = "glmnet", family = "binomial",
#                                trControl = train.control,
#                                tuneLength = 15)
#
# # To extract the best hyperparameters:
# get_best_result <- function(caret_fit) {
#   best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
#   best_result = caret_fit$results[best, ]
#   rownames(best_result) = NULL
#   best_result
# }
# best_param <- get_best_result(caret_fit = trained.mod)
# best_param
# #      alpha       lambda  Accuracy     Kappa AccuracySD   KappaSD
# #1 0.9357143 0.0008972655 0.7420556 0.2757236  0.1260169 0.3631528 # Some other tuning gave slightly better
# # results but I'm having trouble with reproducibility here...
#
# # Final elastic-net model (to observe what variables are kept and their associated coefficients):
# elasticnet.model <- glmnet::glmnet(x = x, y = y, alpha = best_param$alpha,
#                                    lambda = best_param$lambda,
#                                    family = "binomial")
# glmnet::coef.glmnet(object = elasticnet.model, s = best_param$lambda)
#
#
# ### Make predictions using the final elastic net model (glmnet):
# predic.glmnet <- predict(object = elasticnet.model, newx = x,
#                        s = best_param$lambda, type = "class") %>% as.factor()
# mean(predic.glmnet == mydata$high_eff) # 0.776 ±= Accuracy?
# pred.obs <- data.frame(cbind(predic.glmnet, mydata[,1]))
# sum(pred.obs$predic.glmnet == pred.obs$high_eff) # 66 out of 85 observations
# sum(pred.obs$predic.glmnet == 1 & pred.obs$high_eff == 1) # Correctly predicted eradication 10 times
# # out of 23 eradication events (not very good)!
#
# #foldid <- sample(1:10, size = length(y), replace = TRUE)


