context("dataset simulation or summary")

# create initial dataset
data("zfishbeh")
dats1 <- zfishbeh %>%
  dplyr::group_by(endpoint, chemical, conc) %>%
  dplyr::ungroup()

data("zfishdev")
dats2 <- zfishdev %>%
  dplyr::group_by(endpoint, chemical, conc) %>%
  dplyr::ungroup()

test_that("summary, continuous",{
  d <- create_dataset(dats1)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 1))

  # add mask
  datsm <- dats1 %>% dplyr::mutate(mask = 0)
  d <- create_dataset(datsm)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 1))
})

test_that("summary, dichotomous",{
  d <- create_dataset(dats2)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 1))

  # add mask
  datsm <- dats2 %>% dplyr::mutate(mask = 0)
  d <- create_dataset(datsm)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 1))
})

test_that("simulation, dichotomous", {
  d <- create_dataset(dats2, n_samples = 3)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))

  # add mask
  datsm <- dats2 %>% dplyr::mutate(mask = 0)
  d <- create_dataset(datsm, n_samples = 3)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))
})

test_that("simulation, continous", {

  datsm <- dats1 %>% dplyr::mutate(mask = 0)

  # vdata = NULL
  d <- create_dataset(dats1, n_samples = 3)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))

  # vdata = NULL & mask
  d <- create_dataset(datsm, n_samples = 3)
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))


  # vdata != NULL
  d <- create_dataset(dats1, n_samples = 3, vdata = rnorm(100))
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))

  # vdata != NULL & mask
  d <- create_dataset(datsm, n_samples = 3,  vdata = rnorm(100))
  c <- d %>% dplyr::count(endpoint, chemical, conc) %>% dplyr::pull(n)
  expect_true(all(c == 3))

})
