test_that("Test one sample.", {
  
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
  
})


