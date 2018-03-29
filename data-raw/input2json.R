data("zfishbeh")

data_comp <- zfishbeh %>%
  dplyr::group_by(endpoint, chemical) %>%
  dplyr::summarize(concs = list(concs), resps = list(resps)) %>% .[c(1,2),] %>%
  dplyr::ungroup()

inputl <-
  list(
    sample = 10,
    positive_threshold = list(seq(5, 10, by = 5)),
    negative_threshold = list(seq(5, 10, by = 5)),
    seed = 100,
    data = data_comp
  )

inputl %>%
  jsonlite::toJSON(pretty =  TRUE, auto_unbox =  TRUE)

