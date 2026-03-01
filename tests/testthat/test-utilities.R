test_that("MaxEventTime single component.", {
  time <- c(0, 1, 2, 3)
  status <- c(1, 0, 1, 0)
  obs <- MaxEventTime(statuses = status, times = time)
  expect_equal(obs, 2)
})

test_that("MaxEventTime two components.", {
  time <- c(0, 1, 2, 3)
  status <- c(1, 0, 1, 0)
  times <- cbind(c(0, 1, 2, 3), c(3, 4, 5, 6))
  statuses <- cbind(status, status)
  obs <- MaxEventTime(statuses = statuses, times = times)
  expect_equal(obs, 5)
})

test_that("MaxEventTime when no events returns zero.", {
  time <- c(1, 2, 3)
  status <- c(0, 0, 0)
  obs <- MaxEventTime(statuses = status, times = time)
  expect_equal(obs, 0)
})

test_that("GenData returns expected columns.", {
  data <- GenData(n_subj = 10, n_event = 2)
  expect_true(all(c("time1", "time2", "status1", "status2") %in% names(data)))
  expect_equal(nrow(data), 10)
  expect_true(all(data$status1 %in% c(0, 1)))
  expect_true(all(data$status2 %in% c(0, 1)))
})
