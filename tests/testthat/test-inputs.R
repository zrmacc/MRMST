test_that("TwoSample errors when status and time column counts differ.", {
  arm <- rep(c(0, 1), each = 5)
  statuses <- matrix(1, nrow = 10, ncol = 2)
  times <- matrix(1, nrow = 10, ncol = 3)
  expect_error(
    TwoSample(arm = arm, statuses = statuses, times = times),
    "Number of status variables does not match"
  )
})

test_that("TwoSample errors when arm has more than two levels.", {
  arm <- rep(c(0, 1, 2), length.out = 30)
  data <- GenData(n_subj = 30, n_event = 2)
  expect_error(
    TwoSample(
      arm = arm,
      statuses = data[, c("status1", "status2")],
      times = data[, c("time1", "time2")]
    ),
    "Arm should have exactly two levels"
  )
})

test_that("OneSample errors when status and time column counts differ.", {
  statuses <- matrix(1, nrow = 5, ncol = 2)
  times <- matrix(1, nrow = 5, ncol = 3)
  expect_error(
    OneSample(statuses = statuses, times = times),
    "Number of status variables does not match"
  )
})

test_that("OneSample errors when there are zero observations.", {
  statuses <- matrix(integer(0), nrow = 0, ncol = 1)
  times <- matrix(numeric(0), nrow = 0, ncol = 1)
  expect_error(
    OneSample(statuses = statuses, times = times),
    "At least one observation is required"
  )
})

test_that("OneSample errors when tau is invalid.", {
  data <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  expect_error(
    OneSample(statuses = data$status, times = data$time, tau = 0),
    "tau must be a finite positive"
  )
  expect_error(
    OneSample(statuses = data$status, times = data$time, tau = -1),
    "tau must be a finite positive"
  )
})

test_that("OneSample errors when weights length or sign is invalid.", {
  withr::local_seed(1)
  data <- GenData(n_subj = 20, n_event = 2)
  expect_error(
    OneSample(
      statuses = data[, c("status1", "status2")],
      times = data[, c("time1", "time2")],
      weights = c(1, 2, 1)
    ),
    "Length of weights must equal"
  )
  expect_error(
    OneSample(
      statuses = data[, c("status1", "status2")],
      times = data[, c("time1", "time2")],
      weights = c(-1, 1)
    ),
    "Weights must be finite and non-negative"
  )
})

test_that("OneSample errors when statuses or times contain NA.", {
  data <- data.frame(time = c(1, 2, NA), status = c(1, 1, 0))
  expect_error(
    OneSample(statuses = data$status, times = data$time),
    "Missing values are not allowed"
  )
  data2 <- data.frame(time = c(1, 2, 3), status = c(1, NA, 0))
  expect_error(
    OneSample(statuses = data2$status, times = data2$time),
    "Missing values are not allowed"
  )
})
