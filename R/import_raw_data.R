
#' Import the Raw Dataset
#'
#' Import the raw CSV file of the 2020 tarping survey data.
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' mydata <- import_raw_data()
#' }
import_raw_data <- function(){
  readr::read_csv2(file = here::here("data", "raw_data", "data_tarping.csv"), col_names = TRUE)
}
