
#' Import the Raw Dataset
#'
#' Import the raw .csv file of the 2020 knotweed's tarping survey data. \cr To avoid errors, please ensure
#' that there are no ";" in the cells of the CSV (e.g. if comments have been made in the table cells) and
#' that strings within the table are written in lower cases, without spaces or any special characters
#' (you can write in English for instance).
#'
#' @return A tibble (i.e. a kind of improved data.frame). For further information on tibbles, please refer to
#' the `tidyverse` or \link[readr]{readr} documentation.
#' @export
#' @importFrom readr read_delim
#' @importFrom here here
#'
#' @examples
#' \dontrun{
#' mydata <- import_raw_data()
#' }
import_raw_data <- function(){
  x <- readr::read_delim(here::here("mydata", "data_tarping_x.csv"), delim = ";", col_names = TRUE)
  return(x)
}
