
import_raw_data <- function(){
  readr::read_csv2(file = here::here("data", "raw_data", "data_tarping.csv"), col_names = TRUE)
}
