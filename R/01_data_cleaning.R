#########################
##### DATA CLEANING #####
#########################

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
    dplyr::mutate_if(is.character, factor) %>%
    dplyr::mutate_if(is_binary, factor) -> cleaned_data
  return(cleaned_data)
}
