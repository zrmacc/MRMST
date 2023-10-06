test_that("Test max event-time function.", {
  
  # Case of 1 event-time.
  time <- c(0, 1, 2, 3)
  status <- c(1, 0, 1, 0)
  
  obs <- MaxEventTime(statuses = status, times = time)
  exp <- 2
  expect_equal(obs, exp)
  
  # Case of 2 event-time. 
  times <- cbind(c(0, 1, 2, 3), c(3, 4, 5, 6))
  statuses <- cbind(status, status)
  
  obs <- MaxEventTime(statuses = statuses, times = times)
  exp <- 5
  expect_equal(obs, exp)
  
})
