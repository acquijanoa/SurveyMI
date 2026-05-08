suppressPackageStartupMessages({
  library(survey)
})

test_that("multinom_weights applies one nonnegative PSU factor per sampling weight", {
  set.seed(2)
  n_psu <- 8L
  n_per <- 5L
  d <- data.frame(
    strat = rep(c("S1", "S2"), each = n_psu * n_per / 2L),
    psu = rep(seq_len(n_psu), each = n_per),
    w = rep(runif(n_psu, 2, 12), each = n_per)
  )
  des <- svydesign(ids = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  W <- multinom_weights(des, B = 15L, normalize = FALSE)
  expect_true(all(W[["Wt_1"]] >= 0))
  ratio <- W[["Wt_1"]] / d$w
  for (pid in unique(d$psu)) {
    idx <- d$psu == pid
    expect_equal(length(unique(ratio[idx])), 1L)
  }
})

test_that("normalize=TRUE gives replicate columns with mean one", {
  set.seed(3)
  n <- 24L
  dd <- data.frame(
    strat = sample(letters[1:3], n, TRUE),
    psu = sample(10:30, n, TRUE),
    w = runif(n, 1, 9)
  )
  des <- svydesign(ids = ~psu, strata = ~strat, weights = ~w, data = dd, nest = TRUE)
  W <- multinom_weights(des, B = 20L, normalize = TRUE)
  for (k in seq_len(20L)) {
    expect_equal(mean(W[[paste0("Wt_", k)]]), 1)
  }
})

test_that("non-survey objects are rejected", {
  expect_error(
    multinom_weights(mtcars),
    "survey.design"
  )
})
