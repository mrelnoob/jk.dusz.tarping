#########################
##### DATA CLEANING #####
#########################

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


# _____________________________________________________________________
### Creation of function to automatically detect boolean (binary) data:
#' Find Binary Variables
#'
#' @description Automatically detect binary (boolean) variables, even if there is NAs in the data.
#' @param x Any variable (vector, columns etc.) of any length.
#'
#' @return Logical (i.e. a TRUE or FALSE answer).
#' @export
#' @importFrom stats na.omit
#'
#' @examples
#' xx <- c(0, 0, 1, 0, NA)
#' is_binary(xx) # TRUE
is_binary <- function(x) {
  x0 <- na.omit(x)
  length(unique(x0)) %in% 1:2 && all(x0 %in% 0:1)
} # Asks: is the length of the unique elements of x0 (excluding NAs) is 1 or 2 AND all either 0 or 1?





# ________________________________________________________________________________________
### Creation of a function to automatically convert a factor into its exact numeric value:
#' Convert factor values as exact numeric values
#'
#' @description Automatically convert the input `x (factor)` into a `numeric` vector while keeping the
#' exact value of `x`.
#' @details  In \strong{R}, \emph{numeric factors} should be coded as `character` but it is rarely the
#' case because most functions actually prefer `factors`. It's one of the shortcomings of \strong{R}.
#' Consequently, there is no \strong{R} function that converts a \emph{numeric factor} (e.g. a binary
#' variable with only zeros and ones) into a `numeric` vector while keeping the correct value (if you use
#' `as.numeric(myfactor)`, you will get the underlying level codes, not the values as numbers), hence the
#' usefulness of `as.numfactor`.
#'
#' @param x A `factor` variable (may contain NAs).
#'
#' @return A `numeric` vector.
#' @export
#'
#' @examples
#' f <- as.factor(cars$speed)
#' as.numeric(f) # Converts the values into level ranks.
#' as.numfactor(f) # Works!
as.numfactor <- function(x){as.numeric(levels(x))[x]}





###________________
#' Clean Data Types
#'
#' @description The function `clean_my_data` automatically imports the \emph{raw_data} of the knotweed's
#' tarping survey and cleans by transforming character variables into factors, ordinal variables into
#' ordered factors, and boolean/binary variables into factors.
#' @return A cleaned tibble.
#' @note The variables "plantation" and "age" are ordinal variables and I have thus coded them
#' as such (as an ordered factor). However, it might be preferable, from a statistical point of view,
#' to consider it as a numeric variable. The 2nd solution would be more parsimonious (less levels and
#' thus lighter models) but would assume that intervals between each level (between 0 and 1,
#' between 1 and 2, etc.) are equals when it's not a necessarily true assumption (for "plantation",
#' it's probably not). The 1st solution will cause models to add polynomial terms to my factor levels
#' (see marked pages in my web browsers). \cr
#' I also coded "planned_duration" as an ordinal variable but it there's no problem here because I will
#' probably not use it in statistical analyses.
#'
#' @export
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' my_cleaned_data <- clean_my_data()
#' }
clean_my_data <- function(){
  raw_data <- jk.dusz.tarping::import_raw_data()
  # Transform character variables into factors, ordinal variables into ordered factors, and boolean/binary
  # variables into factors:
  raw_data %>%
    dplyr::mutate(planned_duration = factor(x = planned_duration, ordered = TRUE),
                  age = factor(x = age, ordered = TRUE),
                  plantation = factor(x = plantation, ordered = TRUE)) %>%
    dplyr::mutate(geomem = ifelse(grepl("geomem$", fabric_type) | grepl("tarp$", fabric_type) |
                                    grepl("mixed", fabric_type) | grepl("unknown", fabric_type), 1, 0)) %>%
    dplyr::mutate(geotex = ifelse(grepl("geotex$", fabric_type) | grepl("mixed", fabric_type), 1, 0)) %>%
    dplyr::mutate(tarpfix_pierced = ifelse(c(grepl("*staples*", fabric_fixation) | fabric_fixation == "wired_stakes") &
                                             fabric_fixation != "staples_and_taped_patches", 1, 0)) %>%
    dplyr::mutate_if(is.character, factor) %>%
    dplyr::mutate_if(is_binary, factor) -> cleaned_data

  return(cleaned_data)
}





### __________________________
#' Lighter Version of the Data
#'
#' @description This function creates a slightly lighter version of the cleaned knotweed's
#' tarping survey dataset by removing the variables that will not be used in any of the future models
#' (such as purely descriptive variables).
#'
#' @return A tibble.
#' @export
#' @importFrom dplyr select
#'
#' @examples
#' \dontrun{
#' my_light_dataset <- lighten_my_data()
#' }
lighten_my_data <- function(){
  # Importation and cleaning
  toto <- clean_my_data()

  # Lighten my data
  toto %>%
    dplyr::select(xp_id, latitude, longitude, elevation, tarping_date, goals, operation_type, freq_monitoring,
                  slope, difficulty_access, shade, forest, ruggedness, granulometry, obstacles, flood,
                  liner_geomem, agri_geomem, woven_geotex, mulching_geotex, pla_geotex, weedsp_geotex,
                  other_unknown, grammage, thickness, resi_punc, resi_trac,
                  season, maxveg, preparation, levelling, stand_surface, age, fully_tarped, distance,
                  multi_strips, strips_overlap, strips_fixation, staples_distance, fabric_fixation,
                  tarpfix_multimethod, sedicover_height, trench_depth, pierced_tarpinstall, plantation,
                  repairs, add_control, add_control_type,
                  degradation, pb_fixation, pb_durability, pb_trampiercing, pb_vandalism,
                  reg_staples, reg_stripoverlaps, reg_obstacles, reg_holes, reg_plantations, reg_pierced,
                  reg_edges, reg_nearby,
                  tarping_abandoned, tarping_completed, tarping_ongoing, tarping_duration, latest_condition,
                  latest_regrowth, latest_months,
                  eff_expansion, eff_dispersal, eff_vigour, eff_eradication) -> tata
  return(tata)
}
