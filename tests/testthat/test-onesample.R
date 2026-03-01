test_that("One sample single component gives expected AUC.", {
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 1)
  )
  obs <- OneSample(
    statuses = data$status,
    times = data$time
  )
  exp_auc <- 1.0 + 0.8 + 0.6 + 0.4 + 0.2
  expect_equal(obs@AUC, exp_auc)
  expect_equal(obs@K, 1L)
  expect_equal(obs@Tau, 5)
})

test_that("One sample with two components and custom tau.", {
  withr::local_seed(42)
  data <- GenData(n_subj = 50, n_event = 2)
  obs <- OneSample(
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    tau = 3.0
  )
  expect_equal(obs@K, 2L)
  expect_equal(obs@Tau, 3.0)
  expect_true(obs@AUC >= 0)
  expect_true(obs@SE >= 0)
})

test_that("OneSample is reproducible when seed is set.", {
  withr::local_seed(100)
  data <- GenData(n_subj = 40, n_event = 2)
  run1 <- OneSample(
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    seed = 12345
  )
  run2 <- OneSample(
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    seed = 12345
  )
  expect_equal(run1@AUC, run2@AUC)
  expect_equal(run1@SE, run2@SE)
  expect_equal(run1@Table, run2@Table)
})

test_that("TwoSample is reproducible when seed is set.", {
  withr::local_seed(200)
  arm <- rep(c(0, 1), each = 30)
  data <- GenData(
    covariates = data.frame(arm = arm),
    beta_event = 0.5,
    n_event = 2
  )
  run1 <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    seed = 99999
  )
  run2 <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    seed = 99999
  )
  expect_equal(run1@Arm0@AUC, run2@Arm0@AUC)
  expect_equal(run1@Arm0@SE, run2@Arm0@SE)
  expect_equal(run1@Arm1@AUC, run2@Arm1@AUC)
  expect_equal(run1@Arm1@SE, run2@Arm1@SE)
  expect_equal(run1@Contrast, run2@Contrast)
})

test_that("RMST with tau beyond last observed event time.", {
  # All events observed at 1,2,3,4,5; max event time is 5.
  data_all_events <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 1)
  )
  at_max <- OneSample(
    statuses = data_all_events$status,
    times = data_all_events$time,
    tau = 5,
    extend = FALSE
  )
  expect_equal(at_max@AUC, 3)

  # Tau beyond max, extend = FALSE: no extrapolation, same as tau = 5
  beyond_no_extend <- OneSample(
    statuses = data_all_events$status,
    times = data_all_events$time,
    tau = 10,
    extend = FALSE
  )
  expect_equal(beyond_no_extend@AUC, 3)
  expect_equal(beyond_no_extend@Tau, 10)

  # With extend = TRUE and all events, S(5) = 0 so extrapolation adds 0
  beyond_extend_all <- OneSample(
    statuses = data_all_events$status,
    times = data_all_events$time,
    tau = 10,
    extend = TRUE
  )
  expect_equal(beyond_extend_all@AUC, 3)
  expect_equal(beyond_extend_all@Tau, 10)

  # Last observation censored at 5: S(5) = 0.2, so extend adds (10 - 5) * 0.2 = 1
  data_censored <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 0)
  )
  beyond_extend <- OneSample(
    statuses = data_censored$status,
    times = data_censored$time,
    tau = 10,
    extend = TRUE
  )
  expect_equal(beyond_extend@AUC, 4, tolerance = 1e-6)
  expect_equal(beyond_extend@Tau, 10)
})

test_that("Survival curve drops to zero before tau: finite results.", {
  # All events at 1,2,3,4,5; S(5)=0. Tau=10 > 5. No extrapolation.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 1)
  )
  obs <- OneSample(
    statuses = data$status,
    times = data$time,
    tau = 10,
    extend = FALSE
  )
  expect_equal(obs@AUC, 3)
  expect_true(is.finite(obs@SE))
  expect_true(obs@SE >= 0)
  expect_true(all(is.finite(obs@Table$lower)))
  expect_true(all(is.finite(obs@Table$upper)))
})

test_that("All events at same time (survival to zero at first event): finite results.", {
  # Five subjects, all event at time 1. S(1)=0. At-risk can be zero after last time.
  data <- data.frame(
    time = c(1, 1, 1, 1, 1),
    status = c(1, 1, 1, 1, 1)
  )
  obs <- OneSample(
    statuses = data$status,
    times = data$time,
    tau = 2,
    seed = 42
  )
  expect_true(is.finite(obs@AUC))
  expect_true(obs@AUC >= 0)
  expect_true(is.finite(obs@SE))
  expect_true(obs@SE >= 0)
  expect_true(all(is.finite(obs@Table$lower)))
  expect_true(all(is.finite(obs@Table$upper)))
})

test_that("TwoSample contrast CIs use specified alpha.", {
  withr::local_seed(301)
  arm <- rep(c(0, 1), each = 40)
  data <- GenData(
    covariates = data.frame(arm = arm),
    beta_event = 0,
    n_event = 2
  )
  res_05 <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    alpha = 0.05,
    seed = 1
  )
  res_10 <- TwoSample(
    arm = data$arm,
    statuses = data[, c("status1", "status2")],
    times = data[, c("time1", "time2")],
    alpha = 0.10,
    seed = 1
  )
  width_05 <- res_05@Contrast$upper[1] - res_05@Contrast$lower[1]
  width_10 <- res_10@Contrast$upper[1] - res_10@Contrast$lower[1]
  expect_lt(width_10, width_05)
})


