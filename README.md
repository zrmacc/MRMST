# Multiple Outcome Restricted Mean Survival Time

[![R-CMD-check](https://github.com/zrmacc/MRMST/workflows/R-CMD-check/badge.svg)](https://github.com/zrmacc/MRMST/actions)

Zachary R. McCaw <br>
Updated: 2023-10-06



### Overview

For each patient, suppose K time-to-event outcomes are available (e.g. T1, ..., TK), each potentially subject to censoring. An outcome is considered observed (status = 1) or censored (status = 0). For a specified time window (0, tau), the parameter of interest is the sum over outcomes of the restricted mean survival time (RMST) up to tau. That is, the cumulative area under the curve (AUC) for the K survival functions. This package performs inference on the difference and ratio of this parameter when comparing two independent treatment arms.


### Installation

``` r
devtools::install_github(repo = "zrmacc/MRMST")
```


### Data

Example data can be simulated using the `GenData` function:

``` r
library(MRMST)
data <- MRMST::GenData(n_subj = 100, n_event = 2)
head(data)
```

```
##   idx event_rate frailty censor_rate     time1 status1      time2 status2
## 1   1          1       1        0.25 0.4062701       1 0.01789208       0
## 2   2          1       1        0.25 2.8423212       0 0.14081056       1
## 3   3          1       1        0.25 0.9018589       1 0.13934203       1
## 4   4          1       1        0.25 0.8756545       1 0.03020559       1
## 5   5          1       1        0.25 0.3764735       0 0.38447476       1
## 6   6          1       1        0.25 0.7740696       1 0.27826969       0
```

### One-Sample Problem

To estimate the AUC for a single treatment arm in the case of $K=2$ time-to-event outcomes:


``` r
# Single arm estimation.
one_sample <- MRMST::OneSample(
  statuses = data[, c("status1", "status2")],
  times = data[, c("time1", "time2")]
)
show(one_sample)
```

```
## Multiple RMST with 2 components.
##     n   tau k   auc    se lower upper
## 1 100 3.286 2 2.061 0.147 1.773  2.35
```

### Two-Sample Problem

To perform inference on the AUC comparing two treatment arms:


``` r
# Simulate data with a treatment effect.
arm <- rep(c(0, 1), each = 50)
data <- MRMST::GenData(
  covariates = data.frame(arm = arm), 
  beta_event = 1.0,
  n_event = 2
)

# Comparison of two arms.
two_sample <- MRMST::TwoSample(
  arm = data$arm,
  statuses = data[, c("status1", "status2")],
  times = data[, c("time1", "time2")]
)
show(two_sample)
```

```
## Arm 0:
## Multiple RMST with 2 components.
##   arm  n   tau k   auc    se lower upper
## 1   0 50 1.329 2 1.388 0.096   1.2 1.576
## 
## Arm 1:
## Multiple RMST with 2 components.
##   arm  n   tau k   auc    se lower upper
## 1   1 50 1.329 2 0.616 0.063 0.492 0.741
## 
## Contrast:
##    stat    est    se  lower  upper p
## 1 A1-A0 -0.772 0.115 -0.997 -0.546 0
## 2 A1/A0  0.444 0.055  0.348  0.566 0
```
