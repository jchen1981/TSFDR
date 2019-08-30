# tsfdr
Two-stage false discovery rate control for confounder adjustment in genomics studies

## Overview
     The function implements the two-stage false discovery rate control
     for more powerful confounder adjustment in the analysis of genomic
     data. The method is based on the idea that the confounder(s)
     usually affect part of the genomic features, and thus adjusting
     the confounder(s) for ALL genomic features will be
     over-adjustment, leading to reduced statistical power.  The
     two-step procedure starts with performing the unadjusted analysis
     (first step - filtering) to narrow down the list of genomic
     features which are more likely to be affected by either the
     confounder or the variable of interest or both. In the second
     step, we conduct adjusted analysis on these 'top' candidates to
     reduce multiple testing burden. In other words, the unadjusted
     p-values tell us about the probability of the null hypotheses
     being false, and the multiple testing can be focused on those
     promising hypotheses. The procedure is theoretically guaranteed to
     control the false discovery rate while maximizing the power of
     discovery.

## Installation         

```
# install.packages(c("ggplot2", "reshape2", "doMC", "pbivnorm"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(c("limma", "qvalue"))
# install.packages("devtools")
devtools::install_github("jchen1981/TSFDR")
```



### An Example
We illustrate the usage of tsfdr package using simulated data.

```
# Load package
     require(qvalue)
     truth <- c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 850))
     x <- rnorm(50)
     z <- x + rnorm(50)
     z <- scale(z)

     y1 <- x %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y2 <- z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y3 <- x %*% t(rep(0.5, 50)) + z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y <- cbind(y1, y2, y3, matrix(rnorm(50 * 850), nrow = 50))

     # One stage procedure
     obj1 <- summary(lm(y ~ x + z))
     pv1 <- sapply(obj1, function(x) x$coefficient[2, 4])
     qv1 <- qvalue(pv1)$qvalue
     pos1 <- qv1 <= 0.05

     # Two stage procedure
     obj2 <- tsfdr(y, x, z, alpha = 0.05)
     pos2 <- obj2$pos

     sum(pos1 & !truth) / max(sum(pos1), 1)
     sum(pos2 & !truth) / max(sum(pos2), 1)
     sum(pos1 & truth)
     sum(pos2 & truth)
  
```
