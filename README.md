
# dcalasso: Fast divide-and-conquer Cox proportional hazards model with adaptive lasso

[![CRAN](https://www.r-pkg.org/badges/version/dcalasso)](https://CRAN.R-project.org/package=dcalasso)

## Questions that the package addresses

The package answers the following questions:

  - I want to fit a Cox proportional hazards model. But my dataset is
    too large to load to RAM. What should I do?
  - I want to fit a Cox proportional hazards model. But my dataset is
    too large to save as one piece and standard software requires to
    load the dataset as a whole piece. What should I do?
  - I want to do variable selection in Cox model. But my dataset is too
    large to make the computation feasible, too large to load to RAM,
    too large to save as one piece for Cox model variable selection
    method. What should I do? ….

Essentially - what should I do when my data for Cox model w/ or w/o
variable selection are too large?

## Methodology

The dcalasso package aims to fit Cox proportional hazards model to
extremely large, when both n and p are extremely large and n\>\>p. The
method and package have the following features:

  - The package tackles the Cox model fitting for extremely large data
    using the divide-and-conquer strategy, even when the data are too
    large to save as one file.  
  - This approach is able to achieve a fast computation. Meanwhile, it
    returns a set of results that are close to the precision as if the
    model was fitted to the dataset as a whole.  
  - The package could provide model fitting without variable selection
    as well as model fitting with variable selection. It returns results
    for both an unpenalized Cox model without variable selection and an
    adaptive LASSO-penalized variable selection for the Cox model.  
  - The adaptive LASSO variable selection has variable selection
    consistency and asymptotic normality.  
  - The method can be applied to both time-independent and
    time-dependent survival datasets.  
  - The package is flexible in terms of multi-core or single-core
    computation.

The method is detailed
[here](https://academic.oup.com/biostatistics/advance-article-abstract/doi/10.1093/biostatistics/kxz036/5572660?redirectedFrom=fulltext).
Briefly, the method first finds a divide-and-conquer Cox model estimate
without adaptive LASSO penalty by applying the divide-and-conquer
strategy with one-step estimation to the data that are divided into
subsets. Then it finds the divide-and-conquer adaptive LASSO estimate
based on the divide-and-conquer Cox estimate, using least square
approximation.

## Installation

Install development version from GitHub:

``` r
# install.packages("remotes")
install_github("celehs/dcalasso")
```

## Citation

Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac
Kohane, and Tianxi Cai. “[A Fast Divide-and-Conquer Sparse Cox
Regression.](https://arxiv.org/pdf/1804.00735.pdf)”. 2019 Sep 23. kxz036
