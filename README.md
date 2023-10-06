# Multiple Outcome Restricted Mean Survival Time

Zachary R. McCaw <br>
Updated: 2023-10-06



### Overview

Suppose that, for each patient, $K$ time-to-event outcomes $(T_{1}, \dots, T_{K})$ outcomes are available, each potentially subject to censoring. Let $\delta_{k} = 1$ if the $k$th event is observed, and $\delta_{k} = 0$ if the $k$th event is censored. For a specified time window $(0, \tau)$, define the parameter:

$$\theta = \sum_{k=1}^{K}\int_{0}^{\tau}S_{k}(t)dt,$$

where $S_{k}(t) = \mathbb{P}(T_{k} > t)$. $\theta$ can be interpreted as the cumulative area under the curve (AUC) for the $K$ survival functions $(S_{1}, \dots, S_{K})$. This package performance inference on the difference and ratio of $\theta$, comparing two independent treatment arms.


### Installation

```r
devtools::install_github(repo = "zrmacc/MRMST")
```


### Data

Example data can be simulated using the `GenData` function:

```r
library(MRMST)
data <- MRMST::GenData(n_subj = 100, n_event = 2)
head(data)
```

```
##   idx event_rate frailty censor_rate      time1 status1       time2 status2
## 1   1          1       1        0.25 0.12967784       1 0.557029030       1
## 2   2          1       1        0.25 0.47786747       1 0.002097068       1
## 3   3          1       1        0.25 1.30706656       0 0.228389488       1
## 4   4          1       1        0.25 1.31941793       0 0.642464861       1
## 5   5          1       1        0.25 1.82593373       1 0.766493656       1
## 6   6          1       1        0.25 0.03761174       1 1.273435968       1
```

### One-Sample Problem

To estimate the AUC for a single treatment arm in the case of $K=2$ time-to-event outcomes:


```r
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
## 1 100 3.942 2 2.102 0.149  1.81 2.394
```

### Two-Sample Problem

To perform inference on the AUC comparing two treatment arms:


```r
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
## 1   0 50 1.815 2 1.657 0.128 1.406 1.908
## 
## Arm 1:
## Multiple RMST with 2 components.
##   arm  n   tau k   auc    se lower upper
## 1   1 50 1.815 2 0.747 0.066 0.618 0.875
## 
## Contrast:
##    stat    est    se  lower  upper p
## 1 A1-A0 -0.910 0.144 -1.192 -0.628 0
## 2 A1/A0  0.451 0.053  0.358  0.567 0
```






