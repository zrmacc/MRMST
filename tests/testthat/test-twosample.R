test_that("Test two sample, no difference.", {
  
  withr::local_seed(101)
  arm <- rep(c(0, 1), each = 50)
  data <- GenData(covariates = arm, beta_event = 0)
  two_sample <- TwoSample(
    arm = arm,
    statuses = data %>% dplyr::select(status1, status2),
    times = data %>% dplyr::select(time1, time2)
  )
  expect_gt(two_sample@Contrast$p[1], 0.05)
  
})


test_that("Test two sample, expected difference.", {
  
  withr::local_seed(101)
  arm <- rep(c(0, 1), each = 50)
  data <- GenData(covariates = arm, beta_event = 1)
  two_sample <- TwoSample(
    arm = arm,
    statuses = data %>% dplyr::select(status1, status2),
    times = data %>% dplyr::select(time1, time2)
  )
  expect_lt(two_sample@Contrast$p[1], 0.05)
  
})





