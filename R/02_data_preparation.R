############################
##### DATA PREPARATION #####
############################

### IMPORTANT NOTE: to avoid creating notes for unquoted variables, I must add the following code at
# the beginning of every source file (e.g. R/myscript.R) that uses unquoted variables, so in front of
# all my scripts doing any kind of analyses (otherwise, I should always assign each variable, e.g.
# mydata$myvariable, which is quite time consuming and wearisome).
if(getRversion() >= "2.15.1")  utils::globalVariables(c(
  "xp_id", "country", "latitude", "longitude", "elevation", "tarping_date", "planned_duration", "goals",
  "restoration", "operation_type", "multiple_ops", "freq_monitoring", "slope", "difficulty_access",
  "shade", "forest", "ruggedness", "granulometry", "obstacles", "flood", "environment", "fabric_type",
  "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
  "other_unknown", "grammage", "thickness", "resi_punc", "resi_trac", "season", "maxveg", "preparation",
  "levelling", "stand_surface", "age", "fully_tarped", "distance", "multi_strips", "strips_overlap",
  "strips_fixation", "staples_distance", "fabric_fixation", "tarpfix_multimethod", "sedicover_height",
  "trench", "trench_depth", "pierced_tarpinstall", "plantation", "repairs", "add_control",
  "add_control_type", "degradation", "pb_fixation", "pb_durability", "pb_trampiercing", "pb_vandalism",
  "regrowth_during", "reg_staples", "reg_stripoverlaps", "reg_obstacles", "reg_holes", "reg_plantations",
  "reg_pierced", "reg_edges", "reg_nearby", "untarped_regrowth", "tarping_abandoned", "tarping_completed",
  "tarping_ongoing", "tarping_duration", "latest_condition", "latest_regrowth", "latest_months",
  "eff_expansion", "eff_dispersal", "eff_vigour", "eff_eradication"))
# Alternatively, here's another solution (to insert within a given function) when the number of unquoted
# variables is low:
# planned_duration <- age <- plantation <- NULL
# Avoids potentially problematic NOTE during the R CMD check (i.e. devtools::check()) due to the
# fact that these 3 variables are un-quoted (also known as non-standard evaluation (NSE)). If you want
# to find again how I know that, just run a Google search with some of the NOTE message "package: no
# visible binding for global variable" (I read a r-blogger.com post about it).


# _________________________________________
### Creation of a function that will create separate datasets for each response variable (depending on the
# explanatory variables their models will potentially include):

#' Datasets for Models Building
#'
#' @description The `model_datasets` function creates the respective reduced datasets that should be used
#' to model each response variables. For instance, if `response.var = "fixation"`, the function will
#' produce a dataset containing \emph{pb_fixation} as the \strong{response} variable and a subset of
#' variables to be used as \strong{predictors/explanatory variables} (removing all the variables that
#' should not be used to model \emph{pb_fixation}'s variations).
#' @param response.var A character string specifying which modelling dataset should be produced (either:
#' efficiency", "edges", "overlaps", "latest_condition", or "fixation"):
#' * "efficiency" will produce the dataset for the 3 efficiency evaluation variables (namely \emph{
#' eff_eradication}, \emph{eff_expansion} and \emph{eff_vigour});
#' * "edges" will produce the dataset having for response variable the tarping operations that observed
#' regrowth at the edge of the tarped area;
#' * "overlaps" produce the dataset having for response variable the tarping operations that observed
#' regrowth at at strip overlaps;
#' * "latest_condition" will produce the dataset having for response variable the condition (good, bad, etc.)
#' of the fabric at the date of the latest observation;
#' * "fixation" will produce the dataset having for response variable the tarping operations that reported
#' fixation problems during the operation.
#'
#' @return A tibble.
#' @export
#' @import dplyr stringr
#'
#' @examples
#' \dontrun{
#' eff_model <- model_datasets(response.var = "efficiency")
#' }
model_datasets <- function(response.var = c("efficiency", "edges", "overlaps",
                                            "latest_condition", "fixation")){

  ppp <- jk.dusz.tarping::clean_my_data()
  ppp %>%
    dplyr::filter(!stringr::str_detect(operation_type, "crushing_tarping_trial")) -> qqq # Exclude rows which
  # belong to crushing-tarping trials!

  ### For the 3 "efficiency" models
  if (response.var == "efficiency") {
    tapioca <- qqq[,c("xp_id", "latitude", "longitude", "elevation", "goals",
                  "eff_expansion", "eff_vigour", "eff_eradication",
                  "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
                  "obstacles", "flood",
                  "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
                  "other_unknown", "grammage", "thickness",
                  "season", "maxveg", "preparation", "stand_surface", "age", "fully_tarped", "distance",
                  "strips_overlap", "tarpfix_multimethod", "sedicover_height", "trench_depth",
                  "pierced_tarpinstall", "plantation", "repairs", "add_control", "add_control_type",
                  "degradation", "pb_fixation", "pb_durability",
                  "tarping_duration")]
  }

  ### For the "latest_reg_edges" model
  if (response.var == "edges") {
    qqq <- dplyr::mutate(.data = qqq, latest_reg_edges =
                           ifelse(latest_regrowth == "edges" | latest_regrowth == "edges_nearby",
                                  yes = 1, no = 0)) %>%
      dplyr::mutate(.data = qqq, reg_elsewhere =
                      ifelse(reg_staples == 1 | reg_stripoverlaps == 1 | reg_obstacles == 1 |
                               reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0))

    tapioca <- qqq[,c("xp_id", "latitude", "longitude", "elevation", "latest_reg_edges",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "season", "preparation",
      "stand_surface", "age", "fully_tarped", "distance", "strips_overlap",
      "fabric_fixation", "tarpfix_multimethod", "sedicover_height", "trench_depth", "plantation",
      "repairs", "add_control", "add_control_type",
      "degradation", "pb_fixation", "pb_durability",
      "reg_elsewhere",
      "tarping_duration")]
  }

  ### For the "overlaps" model
  if (response.var == "overlaps") {
    qqq <- dplyr::mutate(.data = qqq, reg_elsewhere =
                           ifelse(reg_staples == 1 | reg_edges == 1 | reg_obstacles == 1 |
                                    reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0))

    tapioca <- qqq[,c("xp_id", "latitude", "longitude", "elevation", "reg_stripoverlaps",
      "freq_monitoring", "slope", "difficulty_access","shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "season", "preparation",
      "stand_surface", "age", "fully_tarped", "distance",
      "strips_overlap", "strips_fixation", "staples_distance", "fabric_fixation", "sedicover_height", "plantation",
      "reg_elsewhere","tarping_duration")]
  }

  ### For the "latest_condition" model
  if (response.var == "latest_condition") {
    tapioca <- qqq[,c("xp_id", "latitude", "longitude", "elevation", "latest_condition",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
      "other_unknown", "grammage", "thickness", "resi_punc", "resi_trac",
      "preparation", "levelling", "stand_surface", "age", "fully_tarped", "distance", "strips_overlap",
      "strips_fixation", "staples_distance", "fabric_fixation", "sedicover_height",
      "pierced_tarpinstall", "plantation", "repairs", "add_control",
      "degradation", "regrowth_during", "tarping_duration")]
  }

  ### For the "fixation" model
  if (response.var == "fixation") {
    tapioca <- qqq[,c("xp_id", "latitude", "longitude", "elevation", "pb_fixation",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
      "other_unknown", "grammage", "thickness",
      "preparation", "levelling", "stand_surface", "age", "fully_tarped", "distance", "multi_strips",
      "strips_overlap", "strips_fixation", "staples_distance", "fabric_fixation", "tarpfix_multimethod",
      "sedicover_height", "trench_depth", "plantation", "repairs", "pb_durability", "pb_trampiercing",
      "regrowth_during", "tarping_duration")]
  }
  return(tapioca)
}





# _________________________________________
### Creation of a function that draws univariate boxplots for all numeric variables in a given dataset:

#' Univariate boxplots
#'
#' @description The `uni.boxplot` function draws, within a single panel, an independent boxplot for each
#' numeric (continuous or discrete) variable in a given dataset. It is particularly useful for data
#' exploration (e.g. Zuur \emph{et al.}, 2010). For instance, to simultaneously observe the
#' distributions of all numeric variables or \strong{to detect their univariate outliers}.
#'
#' @details The `uni.boxplot` function only modifies the graphical parameters of the
#' \code{\link[graphics:boxplot]{boxplot}} function in the `graphics` package to match some predefined
#' preferences. Therefore, default values of `uni.boxplot` create nice looking boxplots but retain
#' default \emph{heuristic} aspects of `boxplot` (such as the length of whiskers or the plotting of
#' outliers). These aspects can however be changed as in \code{\link[graphics:boxplot]{boxplot}}. \cr
#' On the other hand, panel parameters are internally controlled using `par`. However, to avoid unforeseen
#' conflicts with other internal parameters, it is not possible to tune panel parameters as we would
#' do with `par`. Instead, parametrization is only possible with the given subset of parameters.
#'
#' @param dataset
#' @param mar
#' @param cex.lab
#' @param font.lab
#' @param bty
#' @param fg
#' @param col.axis
#' @param col.lab
#' @param cex.par
#' @param tcl
#' @param mgp
#' @param oma
#' @param type
#' @param border
#' @param col
#' @param lty
#' @param staplewex
#' @param whisklwd
#' @param boxwex
#' @param boxlwd
#' @param medlwd
#' @param pch
#' @param cex
#' @param ...
#'
#' @return
#' @export
#' @import graphics
#'
#' @examples
#' \dontrun{
#' data(TRUC)
#' etc.
#' }
uni.boxplot <- function(dataset, mar=c(0.5,4.1,1.1,1.5), cex.lab=1, font.lab=2, bty = "n", fg = "gray35",
                        col.axis = "gray35", col.lab = "gray20", cex.par = 0.8, tcl = -0.3,
                        mgp = c(2.4, 0.6, 0), oma = c(1, 0, 0, 0),
                        type = "n", border = "lightcoral", col = "moccasin",  lty = 1, staplewex = 0,
                        whisklwd = 2, boxwex = 0.7, boxlwd = 0.1, medlwd = 2.6, pch = 19, cex = cex, ...){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!

  graphics::par(mfrow= c(ncol.adjust,4), mar=mar, cex.lab=cex.lab, font.lab=font.lab, bty=bty, fg=fg,
      col.axis=col.axis, col.lab=col.lab, cex=cex.par, tcl=tcl, mgp=mgp, oma=oma)
  for (i in c(1:ncol.data)) {
    graphics::boxplot(num.data[,i],ylab =(nam[i]), type=type, border=border, col = col,
            lty=lty, staplewex=staplewex, whisklwd=whisklwd, boxwex=boxwex, boxlwd=boxlwd,
            medlwd=medlwd, pch=pch, cex=cex, ...) }
  }

