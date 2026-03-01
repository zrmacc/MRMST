test_that("Test two sample, no difference.", {
  withr::local_seed(101)
  arm <- rep(c(0, 1), each = 50)
  data <- GenData(
    covariates = data.frame(arm = arm),
    beta_event = 0,
    n_event = 2
  )
  two_sample <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")]
  )
  expect_gt(two_sample@Contrast$p[1], 0.05)
})

test_that("Test two sample, expected difference.", {
  withr::local_seed(101)
  arm <- rep(c(0, 1), each = 50)
  data <- GenData(
    covariates = data.frame(arm = arm),
    beta_event = 1.0,
    n_event = 2
  )
  two_sample <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")]
  )
  expect_lt(two_sample@Contrast$p[1], 0.05)
})





