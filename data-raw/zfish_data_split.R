library(dplyr)


zfishdev <- readRDS("./data-raw/zfishdev.rds")

zfishdev <- zfishdev %>%
  filter(stringr::str_detect(well__substance, 'Caffeine|Saccharin|Dieldrin')) %>%
  mutate(directionality = as.integer(directionality), n_in = as.integer(n_in), N = as.integer(N)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality)

devtools::use_data(zfishdev, overwrite =  TRUE)



zfishbeh <- readRDS("./data-raw/zfishbeh.rds")

zfishbeh <-  zfishbeh %>%
  filter(stringr::str_detect(well__substance, 'Caffeine|Saccharin|Dieldrin')) %>%
  mutate(directionality = as.integer(directionality)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality)

devtools::use_data(zfishbeh, overwrite =  TRUE)
