**Fast divide-and-conquer Cox proportional hazards model with adaptive lasso**

# Descriptions
  
The authors of the [methodological paper](https://arxiv.org/pdf/1804.00735.pdf) are [Yan Wang](https://www.researchgate.net/profile/Yan_Wang374) \<yaw719@mail.harvard.edu\>, [Chuan Hong](https://dbmi.hms.harvard.edu/person/postdoctoral-fellows/chuan-hong), [Nathan Palmer](https://dbmi.hms.harvard.edu/people/nathan-patrick-palmer), [Qian Di](https://scholar.google.com/citations?user=BpMY1OkAAAAJ&hl=en), [Joel Schwartz](https://www.hsph.harvard.edu/joel-schwartz/), [Isaac Kohane](https://dbmi.hms.harvard.edu/people/isaac-samuel-kohane), and [Tianxi Cai](https://www.hsph.harvard.edu/tianxi-cai/)
  
  The paper is archived here https://arxiv.org/pdf/1804.00735.pdf and has been accepted by Biostatistics.
  
  The package answers the following questions:
  
  I want to fit a Cox proportional hazards model. But my dataset is too large to load to RAM. What should I do?
  I want to fit a Cox proportional hazards model. But my dataset is too large to save as one piece. What should I do?
  I want to do variable selection in Cox model. But my dataset is too large to make the computation feasible, too large to load to RAM, too large to save as one piece for Cox model variable selection method. What should I do?
  ....
  
  Essentially - what should I do when my data for Cox model w/ or w/o variable selection are too large?

  The dcalasso package aims to fit extremely large Cox model with and without variable selection via adaptive LASSO, when both n and p are extremely large and n>>p, even if the data cannot be loaded or saved as a whole. The method first finds a divide-and-conquer Cox model estimate without adaptive LASSO penalty by applying the divide-and-conquer strategy with one-step estimation to the data that are divided into subsets. Then it finds the divide-and-conquer adaptive LASSO estimate based on the divide-and-conquer Cox estimate, using least square approximation. This method not only speeds up the computation when the dataset is extremely large, but also makes the computation feasible when the data is impossible to save or load as a whole. In the paper above, we show that the divide-and-conquer adaptive LASSO estimator has variable selection consistency and asymptotic normality property as the standard adaptive lasso estimator.
  
  The package can fit both adaptive lasso for time-independent Cox proportional hazards model and adaptive lasso for time-dependent Cox proportional hazards model. The package is flexible for data that can or cannot be saved or loaded as a whole. The package is also flexible in terms of using multi-core or single-core computation.
  
  The dcalasso is specialized in computing the divide-and-conquer adaptive lasso with fast computation, with the following novel features: (1) a divide-and-conquer strategy is applied for the computation of the initial unpenalized Cox proportional hazards model by splitting the observations into <tt>K</tt> chunks and processing each separately; (2) a fast linearization (least square approximation) is applied in the estimation of the shrinkage step further reducing the computational burden.
  
  When the interest is to estimate the Cox model without adaptive lasso shrinkage, this package is also useful and the user can simply take the unpenalized estimate result. Note that when n>>p, the computation for the adaptive lasso estimate takes a very small amount of time, i.e. the dcalasso is still computationally advantageous over the standard approach to fitting Cox model for extremely large dataset when the penalized estimate is not of interest.

# Installation

```r
# install.packages("devtools")
devtools::install_github("celehs/dcalasso")
library(dcalasso)
```

# Key function

```r
dcalasso(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subsets, na.action, 
  offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1)
```

where <tt>formula</tt> is a formula of a Cox model (see <tt>coxph</tt>), <tt>family</tt> specifies the family of the outcome (<tt>cox.ph</tt> currently), <tt>data</tt> specifies the dataset, if <tt>data</tt> is too large to load into the memory, <tt>data.rds</tt> specifies a vector of file locations where the data are split and saved, <tt>weights, subsets, na.action, offset</tt> are the same as the corresponding arguments in <tt>coxph</tt>, <tt>lambda</tt> is the penalization parameter for adaptive lasso (see <tt>glmnet</tt>), <tt>gamma</tt> is the exponent for the penalization of the adaptive penalty (see <tt>glmnet</tt>), <tt>K</tt> is the number of split which will be overwritten by <tt>length(data.rds)</tt> if <tt>data.rds</tt> is specified, <tt>iter.os</tt> is the number of iterations for the computation of the initial unpenalized estimator (default 2; a larger <tt>iter.os</tt> will result in longer computation but an unpenalized estimator closer to the results at convergence, and <tt>ncores</tt> is the number of cores used in the computation.

# Citation

Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "[A Fast Divide-and-Conquer Sparse Cox Regression.](https://arxiv.org/pdf/1804.00735.pdf)" arXiv preprint arXiv:1804.00735 (2018).
  
