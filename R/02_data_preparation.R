############################
##### DATA PREPARATION #####
############################

### IMPORTANT NOTE: to avoid creating notes for unquoted variables, I must add the following code at
# the beginning of every source file (e.g. R/myscript.R) that uses unquoted variables, so in front of
# all my scripts doing any kind of analyses (otherwise, I should always assign each variable, e.g.
# mydata$myvariable, which is quite time consuming and wearisome).
if(getRversion() >= "2.15.1")  utils::globalVariables(c(
  "manager_id", "xp_id", "country", "latitude", "longitude", "elevation", "tarping_date",
  "planned_duration", "goals",
  "restoration", "operation_type", "multiple_ops", "freq_monitoring", "slope", "difficulty_access",
  "shade", "forest", "ruggedness", "granulometry", "obstacles", "flood", "environment", "fabric_type",
  "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
  "other_unknown", "grammage", "thickness", "resi_punc", "resi_trac", "season", "maxveg", "preparation",
  "levelling", "stand_surface", "age", "fully_tarped", "distance", "multi_strips", "strips_overlap",
  "stripsoverlap_ok", "tarpfix_pierced", "geomem", "geotex", "uprootexcav", "efficiency", "high_eff",
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
    dplyr::filter(!stringr::str_detect(operation_type, "crushing_tarping_trial")) %>%
    dplyr::mutate(uprootexcav = ifelse(preparation == "excavation" |
                                         preparation == "uprooting", 1, 0)) -> qqq # Exclude rows which
  # belong to crushing-tarping trials, and creates a new variable called "uprootexcav"!

  ### For the 3 "efficiency" models
  if (response.var == "efficiency") {
    qqq <- dplyr::mutate(.data = qqq,
                         high_eff = ifelse(stringr::str_detect(latest_regrowth, stringr::regex("few.", dotall = TRUE)) |
                                                          eff_eradication == "1", 1, 0)) %>%
      dplyr::mutate_if(is_binary, factor)
    # With str_detect and regex, I matched all obs. with latest_regrowth that begins with "few"
    # (quite similar to grep)!

    qqq <- dplyr::mutate(.data = qqq,
                         efficiency = rowMeans(x = qqq[,80:82], na.rm = FALSE))


    tapioca <- qqq[,c("manager_id", "xp_id", "latitude", "longitude", "elevation", "goals",
                  "efficiency", "eff_eradication", "high_eff",
                  "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness",
                  "granulometry", "obstacles", "flood",
                  "geomem", "geotex", "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex",
                  "pla_geotex", "weedsp_geotex", "other_unknown", "grammage", "thickness",
                  "maxveg", "uprootexcav", "stand_surface", "age", "fully_tarped", "distance", "tarping_duration",
                  "stripsoverlap_ok", "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "trench_depth",
                  "pierced_tarpinstall", "plantation", "repairs", "add_control", "add_control_type",
                  "degradation", "pb_fixation", "pb_durability")]
  }

  ### For the "latest_reg_edges" model
  if (response.var == "edges") {
    qqq <- dplyr::mutate(.data = qqq, latest_reg_edges =
                           ifelse(latest_regrowth == "edges" | latest_regrowth == "edges_nearby",
                                  yes = 1, no = 0)) %>%
      dplyr::mutate(.data = qqq, reg_elsewhere =
                      ifelse(reg_staples == 1 | reg_stripoverlaps == 1 | reg_obstacles == 1 |
                               reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0))

    tapioca <- qqq[,c("manager_id", "xp_id", "latitude", "longitude", "elevation", "latest_reg_edges",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "geomem", "geotex",
      "season", "uprootexcav",
      "stand_surface", "age", "fully_tarped", "distance", "tarping_duration", "stripsoverlap_ok","strips_overlap",
      "fabric_fixation", "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "trench_depth", "plantation",
      "repairs", "add_control", "add_control_type",
      "degradation", "pb_fixation", "pb_durability",
      "reg_elsewhere")]
  }

  ### For the "overlaps" model
  if (response.var == "overlaps") {
    qqq <- dplyr::mutate(.data = qqq, reg_elsewhere =
                           ifelse(reg_staples == 1 | reg_edges == 1 | reg_obstacles == 1 |
                                    reg_holes == 1 | reg_plantations == 1 | reg_pierced == 1, yes = 1, no = 0))

    tapioca <- qqq[,c("manager_id", "xp_id", "latitude", "longitude", "elevation", "reg_stripoverlaps",
      "freq_monitoring", "slope", "difficulty_access","shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "geomem", "geotex",
      "season", "uprootexcav",
      "stand_surface", "age", "fully_tarped", "distance","tarping_duration",
      "stripsoverlap_ok", "strips_overlap", "strips_fixation", "staples_distance", "fabric_fixation", "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "plantation",
      "reg_elsewhere")]
  }

  ### For the "latest_condition" model
  if (response.var == "latest_condition") {
    tapioca <- qqq[,c("manager_id", "xp_id", "latitude", "longitude", "elevation", "latest_condition",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "geomem", "geotex", "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
      "other_unknown", "grammage", "thickness", "resi_punc", "resi_trac",
      "uprootexcav", "levelling", "stand_surface", "age", "fully_tarped", "distance", "tarping_duration",
      "stripsoverlap_ok", "strips_overlap", "strips_fixation", "staples_distance", "fabric_fixation", "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height",
      "pierced_tarpinstall", "plantation", "repairs", "add_control",
      "degradation", "regrowth_during")]
  }

  ### For the "fixation" model
  if (response.var == "fixation") {
    tapioca <- qqq[,c("manager_id", "xp_id", "latitude", "longitude", "elevation", "pb_fixation",
      "freq_monitoring", "slope", "difficulty_access", "shade", "forest", "ruggedness", "granulometry",
      "obstacles", "flood",
      "geomem", "geotex", "liner_geomem", "agri_geomem", "woven_geotex", "mulching_geotex", "pla_geotex", "weedsp_geotex",
      "other_unknown", "grammage", "thickness",
      "uprootexcav", "levelling", "stand_surface", "age", "fully_tarped", "distance", "tarping_duration",
      "multi_strips", "stripsoverlap_ok", "strips_overlap", "strips_fixation", "staples_distance", "fabric_fixation",
      "tarpfix_multimethod", "tarpfix_pierced", "sedicover_height", "trench_depth", "plantation",
      "repairs", "pb_durability", "pb_trampiercing",
      "regrowth_during")]
  }
  return(tapioca)
}





# _________________________________________
### Creation of a function that draws univariate boxplots for all numeric variables in a given dataset:

#' Univariate boxplots
#'
#' @description The `uni.boxplots` function draws, within a single panel, an independent boxplot for each
#' numeric (continuous or discrete) variable in a given dataset. It is particularly useful for data
#' exploration (e.g. Zuur \emph{et al.}, 2010). For instance, to simultaneously observe the
#' distributions of all numeric variables or \strong{to detect their univariate outliers}.
#'
#' @details The `uni.boxplots` function only modifies the graphical parameters of the
#' \code{\link[graphics:boxplot]{boxplot}} function in the `graphics` package to match some predefined
#' preferences. Therefore, default values of `uni.boxplots` create nice looking boxplots but retain
#' default \emph{heuristic} aspects of `boxplot` (such as the length of whiskers or the plotting of
#' outliers). These aspects can however be changed as in \code{\link[graphics:boxplot]{boxplot}}. \cr
#' On the other hand, panel parameters are internally controlled using `par`. However, to avoid unforeseen
#' conflicts with other internal parameters, it is not possible to tune panel parameters as we would
#' do with `par`. Instead, parametrization is only possible with the given subset of parameters.
#'
#' @note To avoid \emph{recursive argument errors}, internal arguments should be called using upper case
#' letters (e.g. CEX.LAB = 0.9) whereas other arguments from the `boxplot` function should be called with
#' their normal case writing (e.g. outline = FALSE)!
#'
#' @param dataset The input dataset containing all variables to be plotted. It may contain all kinds of
#' variables, the `uni.boxplots` function will automatically detect and plot numeric variables (columns).
#' @param ... Any other parameter that can be incorporated in \code{\link[graphics:boxplot]{boxplot}}.
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting of `cex`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param TYPE The type of boxplot to draw. Default is "n".
#' @param BORDER An optional vector of colors for the outlines of the boxplots. The values in border
#' are recycled if the length of border is less than the number of plots. Default is "lightcoral".
#' @param COL If col is non-null it is assumed to contain colors to be used to colour the bodies of
#' the boxplots. Default is "moccasin".
#' @param LTY The line type. Line types can either be specified as an integer (0=blank, 1=solid (default),
#' 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses ‘invisible lines’
#' (i.e., does not draw them).
#' @param STAPLEWEX Staple line width expansion, proportional to box width. Default is 0.
#' @param WHISKLWD Whisker line width expansion. Default is 2.
#' @param BOXWEX A scale factor to be applied to all boxes. When there are only a few groups, the
#' appearance of the plot can be improved by making the boxes narrower. Default is 0.7.
#' @param BOXLWD Width of boxplot outer lines. Default is 0.1.
#' @param MEDLWD Width of the median line. Default is 2.6.
#' @param PCH The type of points to be drawn for outliers. Default is 19. See \code{\link[graphics:points]{points}}
#' for possible values and their interpretation.
#'
#' @return A panel of univariate boxplots.
#' @export
#' @import graphics
#'
#' @examples
#' data("mtcars")
#' uni.boxplots(dataset = mtcars)
uni.boxplots <- function(dataset, MAR=c(0.5,4.1,1.1,1.5), CEX.LAB=1, FONT.LAB=2, BTY = "n", FG = "gray35",
                        COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.8, TCL = -0.3,
                        MGP = c(2.4, 0.6, 0), OMA = c(1, 0, 0, 0),
                        TYPE = "n", BORDER = "lightcoral", COL = "moccasin",  LTY = 1, STAPLEWEX = 0,
                        WHISKLWD = 2, BOXWEX = 0.7, BOXLWD = 0.1, MEDLWD = 2.6, PCH = 19, ...){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!

  graphics::par(mfrow= c(ncol.adjust,4), mar=MAR, cex.lab=CEX.LAB, font.lab=FONT.LAB, bty=BTY, fg=FG,
      col.axis=COL.AXIS, col.lab=COL.LAB, cex=CEX.PAR, tcl=TCL, mgp=MGP, oma=OMA)
  for (i in c(1:ncol.data)) {
    graphics::boxplot(num.data[,i],ylab =(nam[i]), type=TYPE, border=BORDER, col = COL,
            lty=LTY, staplewex=STAPLEWEX, whisklwd=WHISKLWD, boxwex=BOXWEX, boxlwd=BOXLWD,
            medlwd=MEDLWD, pch=PCH, cex=0.7, ...) }
}





# _________________________________________
### Creation of a function that draws univariate Cleveland dotplots for all numeric variables in a
### given dataset:

#' Univariate Cleveland dotplots
#'
#' @description The `uni.dotplots` function draws, within a single panel, an independent Cleveland dotplot
#' (i.e. plotting value against rank) for each numeric (continuous or discrete) variable in a given
#' dataset. It is particularly useful for data exploration (e.g. Zuur \emph{et al.}, 2010). For instance,
#' to simultaneously observe the distributions of all numeric variables or \strong{to detect their
#' univariate outliers}.
#'
#' @details The `uni.dotplots` function only modifies the graphical parameters of the
#' \code{\link[graphics:plot.default]{plot}} function in the `graphics` package to match some predefined
#' preferences. Therefore, default values of `uni.dotplots` create nice looking dotplots but retain
#' some aspects of the original `plot` function. These aspects can however be changed as in
#' \code{\link[graphics:plot.default]{plot}}. \cr
#' On the other hand, panel parameters are internally controlled using `par`. However, to avoid unforeseen
#' conflicts with other internal parameters, it is not possible to tune panel parameters as we would
#' do with `par`. Instead, parametrization is only possible with the given subset of parameters.
#'
#' @note To avoid \emph{recursive argument errors}, internal arguments should be called using upper case
#' letters (e.g. CEX.LAB = 0.9) whereas other arguments from the `plot` function should be called with
#' their normal case writing (e.g. sub = "My subtitle")!
#'
#' @param dataset The input dataset containing all variables to be plotted (must be a `data.frame` with
#' at least 2 variables). It may contain all kinds of columns, the `uni.dotplots` function will
#' automatically detect and plot numeric variables (columns).
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting
#' of `CEX.PAR`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param LAB A numerical vector of the form `c(x, y, len)` which modifies the default way that axes
#' are annotated. The values of `x` and `y` give the (approximate) number of tickmarks on the x and y
#' axes and len specifies the label length. The default is `c(5, 5, 7)`. Note that this only affects
#' the way the parameters xaxp and yaxp are set when the user coordinate system is set up, and is
#' not consulted when axes are drawn. `len` is unimplemented in `R`.
#' @param COL.PCH The color to be used for points. Default is "lightcoral".
#' @param PCH The type of points to be drawn. Default is 19. See \code{\link[graphics:points]{points}}
#' for possible values and their interpretation.
#' @param COL.GRID The color of the background grid. Default is "lavender".
#' @param NX The number of lines for the grid in x. Default value is 5.
#' @param NY The number of lines for the grid in Y. Default value is 9.
#' @param LTY The type of lines to be drawn in the background grid. Default value is 6.
#' @param ... Any other parameter that can be incorporated in \code{\link[graphics:plot.default]{plot}}.
#'
#' @return A panel of univariate dotplots.
#' @export
#' @import graphics
#'
#' @examples
#' data("mtcars")
#' uni.dotplots(dataset = mtcars, COL.GRID = "lightblue", LTY = 1, NX = 10, NY = 20)
uni.dotplots <- function(dataset, MAR=c(3,2,0.5,1.5), CEX.LAB = 1.2, FONT.LAB = 2, BTY = "n",
                         FG = "gray35", COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.6,
                         TCL = -0.3, MGP = c(1.7, 0.6, 0.1), OMA = c(1, 0, 1, 0), LAB = c(5, 10, 7),
                         COL.PCH = "lightcoral", PCH = 19, COL.GRID = "lavender", NX = 5, NY = 9, LTY = 6,
                         ...){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!
  num.data <- as.matrix(num.data)

  graphics::par(mfrow= c (ncol.adjust,4), mar=MAR, cex.lab = CEX.LAB, font.lab=FONT.LAB, bty = BTY, fg = FG,
      col.axis = COL.AXIS, col.lab = COL.LAB, cex = CEX.PAR, tcl = TCL,
      mgp = MGP, oma = OMA, lab = LAB)
  for (i in c(1:ncol(num.data))) {
    graphics::plot(x = num.data[,i], y = 1:length(num.data[,i]), type = "p", xlab = nam[i], ylab = "",
         col = COL.PCH, pch = PCH, panel.first = {
           grid(col=COL.GRID,nx = NX,ny = NY, lty = LTY)
         }, ...) }
  # Here, the argument panel.first={} is used to draw the grid first, so behind the points!
}





# _________________________________________
### Creation of a function that generate random samples from various distributions based on my variables
### parameters:

#' Random sample simulation
#'
#' @description The `uni.simudistrib` function **automatically generates 5 Cleveland dotplots** of random
#' samples from different distributions (either `normal`, `log-normal`, or `poisson`) based on the
#' parameters of the variables in the input `data.frame` or `matrix` (see \emph{Details}).  \cr
#' This function is useful to see whether variables' extreme values are actual outliers or whether they
#' lie in a range of values possible for a random sample drawn from a `normal`, `log-normal`, or a `poisson`
#' distribution. \emph{In fine}, it may help determine if the original variable can be approximated by
#' these distribution with or without a transformation.
#'
#' @details The `uni.simudistrib` function extracts some key parameters from the input variables (sample
#' size, mean and standard deviation) and generates random samples based on these parameters. For instance,
#' if `simu.var` contains \emph{i} variables `X1`, `X2`, ... `Xi` and if `distribution = "normal"`, the
#' function will return a panel of \emph{i}x5 plots:
#'
#' * The 1st row will contain five dotplots for five random samples with \emph{n} = `length(X1)` and drawn
#' from a Normal distribution with the same mean and standard deviation as `X1`.
#' * The 2nd row will contain five dotplots for five random samples with \emph{n} = `length(X2)` and drawn
#' from a Normal distribution with the same mean and standard deviation as `X2`.
#' * Etc.
#'
#' \strong{Warning}: the function may fail for `log-normal` and `poisson` distributions if input
#' variables contain negative values (because these distributions are by definition positive). Additionally,
#' if `distribution = "poisson"`, the resulting plots will return \strong{integer} values as Poisson
#' is a discrete probability distribution.
#'
#' @param simu.var A `data.frame` or a `matrix`. For obvious layout and readability reasons, `simu.var`
#' should not include too many variables (`p` < 12 is advised) if the plot is to be printed in a page
#' or window of limited dimensions. For any use in HTML documents (e.g. with `RMarkdown`), the number of
#' input variables should not be a problem.
#' @param distribution Either `"normal"`, `"log-normal"` or `"poisson"` (case sensitive).
#'
#' @return A panel of `p`*5 plots, where `p` is the number of variables in `simu.var`.
#' @export
#' @import graphics
#' @importFrom stats na.omit sd rnorm rlnorm rpois
#'
#' @examples
#' uni.simudistrib(simu.var = iris[,1:4], distribution = "normal")
uni.simudistrib <- function(simu.var, distribution){
  if (!is.data.frame(simu.var) && !is.matrix(simu.var)) {
    warning("The input dataset (i.e. simu.var) must either be a data.frame or a matrix!")
  }
  if(missing(distribution)){
    stop("Please specify the desired distribution. See ?uni.simudistrib for details.")
  }

  ncol.data <- ncol(simu.var)
  nam <- names(simu.var)
  num.data <- as.matrix(simu.var)

  graphics::par(mfrow = c(ncol.data,5), mar = c(3,2,0.1,1.5), cex.lab = 1, font.lab=2, bty = "n",
                fg = "gray35", col.axis = "gray35", col.lab = "gray20", cex = 0.6, tcl = -0.3,
                mgp = c(1.8, 0.6, 0.1), oma = c(0.2, 0.2, 0, 0), lab = c(5, 10, 7))
  for (i in c(1:ncol(num.data))) {
    mu <- mean(stats::na.omit(num.data[,i]))
    std <- stats::sd(stats::na.omit(num.data[,i]))
    n <- length(num.data[,i])
    rrr <- NULL

    if (distribution == "normal") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rnorm(n = n, mean = mu, sd = std)
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rnorm based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }

    if (distribution == "log-normal") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rlnorm(n = n,
                           meanlog = log(mu^2 / sqrt(std^2 + mu^2)),
                           sdlog = sqrt(log(1 + (std^2 / mu^2)))) # That's how you use rlnorm!!!!!!
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rlnorm based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }

    if (distribution == "poisson") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rpois(n = n, lambda = mu)
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rpois based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }
  }
}





# _________________________________________
### Creation of a function that generate histograms for each numeric variable:

#' Histogram Panel
#'
#' @description The `uni.histograms` function draws, within a single panel, an independent histogram
#' or each numeric (continuous or discrete) variable in a given dataset. It is particularly useful
#' for data exploration (e.g. Zuur \emph{et al.}, 2010). For instance, to simultaneously observe
#' the distributions of all numeric variables and determine which one will require transformation.
#'
#' @param dataset The input dataset containing all variables to be plotted (must be a `data.frame` with
#' at least 2 variables). It may contain all kinds of columns, the `uni.histograms` function will
#' automatically detect and plot numeric variables (columns).
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting
#' of `CEX.PAR`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param LAB A numerical vector of the form `c(x, y, len)` which modifies the default way that axes
#' are annotated. The values of `x` and `y` give the (approximate) number of tickmarks on the x and y
#' axes and len specifies the label length. The default is `c(5, 5, 7)`. Note that this only affects
#' the way the parameters xaxp and yaxp are set when the user coordinate system is set up, and is
#' not consulted when axes are drawn. `len` is unimplemented in `R`.
#' @param BREAKS One of:
#' * a vector giving the breakpoints between histogram cells,
#' * a function to compute the vector of breakpoints,
#' * a single number giving the number of cells for the histogram,
#' * a character string naming an algorithm to compute the number of cells (see
#' \code{\link[graphics:hist]{hist}}),
#' * a function to compute the number of cells.
#'
#' In the last three cases the number is a suggestion only; as the breakpoints will be set to pretty
#' values, the number is limited to 1e6 (with a warning if it was larger). If breaks is a function,
#' the x vector is supplied to it as the only argument (and the number of breaks is only limited by
#' the amount of available memory).
#' @param COL The color of the bar of the histograms (bins).
#' @param BORDER The color of the border of the bars.
#'
#' @return A panel of histograms.
#' @import graphics
#' @export
#'
#' @examples
#' uni.histograms(dataset = iris[,1:4])
uni.histograms <- function(dataset, MAR=c(3,2,0.5,1.5), CEX.LAB = 1.2, FONT.LAB = 2, BTY = "n",
                     FG = "gray35", COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.6,
                     TCL = -0.3, MGP = c(1.7, 0.6, 0.1), OMA = c(1, 0, 1, 0), LAB = c(5, 10, 7),
                     BREAKS = 10, COL = "moccasin", BORDER = "white"){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!
  num.data <- as.matrix(num.data)

  graphics::par(mfrow= c (ncol.adjust,4), mar=MAR, cex.lab = CEX.LAB, font.lab=FONT.LAB, bty = BTY, fg = FG,
                col.axis = COL.AXIS, col.lab = COL.LAB, cex = CEX.PAR, tcl = TCL,
                mgp = MGP, oma = OMA, lab = LAB)
  for (i in c(1:ncol(num.data))) {
    graphics::hist(num.data[,i], breaks = BREAKS, col = COL, border = BORDER,
                   main = "", xlab = nam[i], ylab = "")
  }
}
