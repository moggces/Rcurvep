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

## the zfish development complete
zfishdev <- readRDS("./data-raw/zfishdev.rds")
zfishdev_all <- zfishdev %>%
  mutate(directionality = as.integer(directionality), n_in = as.integer(n_in), N = as.integer(N)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality) %>%
  filter(stringr::str_detect(endpoint, "96"))
devtools::use_data(zfishdev_all, overwrite =  TRUE)
#saveRDS(zfishdev_all, here::here("data-raw", "zfishdev_complete.rds"))

## the zfish behavior data
zfishbeh <- readRDS("./data-raw/zfishbeh.rds")

zfishbeh <-  zfishbeh %>%
  filter(stringr::str_detect(well__substance, 'Caffeine|Saccharin|Dieldrin')) %>%
  mutate(directionality = as.integer(directionality)) %>%
  rename(endpoint = id, chemical = well__substance) %>%
  select(-directionality)


devtools::use_data(zfishbeh, overwrite =  TRUE)

## zfish development data results

# percentd <- readRDS("./data-raw/biobide32_percent_boot_1000.rds")
# percentd <- percentd %>% bind_rows(.id = "endpoint")
# act <- extract_curvep_data(percentd, "act")
# hl <- extract_curvep_data(percentd, "concs_hl")
# zfishdev_act <- list(act, hl) %>%
#   purrr::reduce(inner_join) %>%
#   dplyr::filter(endpoint == "percent_affected_96" | endpoint == "percent_mortality_96") %>%
#   tidyr::separate(.data$dduid, c("endpoint", "chemical", "directionality"), sep = "#") %>%
#   dplyr::select(-directionality) %>%
#   dplyr::mutate(direction = 1) %>%
#   dplyr::rename(threshold = thres)
#
# class(zfishdev_act) <- c("rcurvep_out", class(zfishdev_act))
# devtools::use_data(zfishdev_act, overwrite =  TRUE)

set.seed(300)
acts <- run_curvep_batch(zfishdev_all,
                         directionality = 1, n_sample = 100,
                         threshold = seq(5, 95, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE),
                         simplify_output = TRUE)

zfishdev_act <- acts
devtools::use_data(zfishdev_act, overwrite =  TRUE)
