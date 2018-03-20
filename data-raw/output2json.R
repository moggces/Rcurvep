# prepare the dataset
data("zfishbeh")
threshold_range <- seq(5, 10, by = 5) # usually seq(5, 95, by = 5)
n_sample_test <- 2 # usually for 1000
outd <- run_curvep_job(zfishbeh,
      directionality = 1, n_sample = n_sample_test,
      threshold = threshold_range, other_paras = list(CARR = 20, TrustHi = TRUE))


# summary data

## simple version (works)
sum_act <- extract_curvep_data(outd, "summary")
sum_json <- sum_act %>% jsonlite::toJSON()
sum_act_re <- sum_json %>% jsonlite::fromJSON()
dplyr::all_equal(sum_act, sum_act_re)

## verbose version (works)
sum_act <- extract_curvep_data(outd, "summary")
sum_json <- sum_act %>% jsonlite::serializeJSON()
sum_act_re <- sum_json %>% jsonlite::unserializeJSON()

# complete data

## verbose version ONLY (works)
outd_json <- outd %>% jsonlite::serializeJSON()
outd_re <- outd_json %>% jsonlite::unserializeJSON()

## test
outd_re_sum <- extract_curvep_data(outd_re, "summary")
