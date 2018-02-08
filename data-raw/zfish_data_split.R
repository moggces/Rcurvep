library(dplyr)
library(Rcurvep)
## the zfish development data
zfishdev <- readRDS("./data-raw/zfishdev.rds")

zfishdev <- zfishdev %>%
  filter(stringr::str_detect(well__substance, 'Caffeine|Saccharin|Dieldrin')) %>%
  mutate(directionality = as.integer(directionality), n_in = as.integer(n_in), N = as.integer(N)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality)

devtools::use_data(zfishdev, overwrite =  TRUE)


## the zfish behavior data
zfishbeh <- readRDS("./data-raw/zfishbeh.rds")

zfishbeh <-  zfishbeh %>%
  filter(stringr::str_detect(well__substance, 'Caffeine|Saccharin|Dieldrin')) %>%
  mutate(directionality = as.integer(directionality)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality)


devtools::use_data(zfishbeh, overwrite =  TRUE)

## zfish development data results

percentd <- readRDS("./data-raw/biobide32_percent_boot_1000.rds")
percentd <- percentd %>% bind_rows(.id = "endpoint")
act <- extract_curvep_data(percentd, "act")
hl <- extract_curvep_data(percentd, "concs_hl")
zfishdev_act <- list(act, hl) %>%
  purrr::reduce(inner_join) %>%
  dplyr::filter(endpoint == "percent_affected_96" | endpoint == "percent_mortality_96") %>%
  dplyr::select(-endpoint)

devtools::use_data(zfishdev_act, overwrite =  TRUE)
