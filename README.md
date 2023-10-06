# Multiple Outcome Restricted Mean Survival Time

Zachary R. McCaw <br>
Updated: 2023-10-06



### Overview

Suppose that, for each patient, $K$ time-to-event outcomes $(T_{1}, \dots, T_{K})$ outcomes are available, each potentially subject to censoring. Let $\delta_{k} = 1$ if the $k$th event is observed, and $\delta_{k} = 0$ if the $k$th event is censored. For a specified time window $(0, \tau)$, define the parameter:

$$
\theta = \sum_{k=1}^{K}\int_{0}^{\tau}S_{k}(t)dt,
$$
where $S_{k}(t) = \mathbb{P}(T_{k} > t)$. $\theta$ can be interpreted as the cumulative area under the curve (AUC) for the $K$ survival functions $(S_{1}, \dots, S_{K})$. This package performance inference on the difference and ratio of $\theta$, comparing two independent treatment arms.


### Data

Example data can be simulated using the `GenData` function:

```r
library(MRMST)
data <- MRMST::GenData(n_subj = 100, n_event = 2)
head(data)
```

```
##   idx event_rate frailty censor_rate      time1 status1      time2 status2
## 1   1          1       1        0.25 0.02403944       1 0.42087232       1
## 2   2          1       1        0.25 1.46719085       1 2.81605665       0
## 3   3          1       1        0.25 0.27992212       1 0.07674163       1
## 4   4          1       1        0.25 0.87976750       1 0.38247849       0
## 5   5          1       1        0.25 1.56244688       1 1.82365773       1
## 6   6          1       1        0.25 2.33398972       0 0.78085565       1
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
##     tau k  auc    se lower upper
## 1 3.972 2 1.79 0.134 1.528 2.053
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
##   arm   tau k   auc    se lower upper
## 1   0 2.496 2 1.893 0.139 1.622 2.165
## 
## Arm 1:
## Multiple RMST with 2 components.
##   arm   tau k   auc    se lower upper
## 1   1 2.496 2 0.794 0.078 0.641 0.946
## 
## Contrast:
##    stat    est    se  lower  upper p
## 1 A1-A0 -1.100 0.159 -1.411 -0.788 0
## 2 A1/A0  0.419 0.051  0.330  0.533 0
```






