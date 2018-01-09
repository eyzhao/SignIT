# SignIT
Mutation Signatures in Individual Tumours

## Installation

### Development Version

To install from GitHub, first make sure that you can install the rstan package and C++ toolchain by following [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). Once rstan is successfully installed, you can install SignIT from GitHub using the devtools package by executing the following in R:

```{r}
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("eyzhao/SignIT", recompile = TRUE, build_vignettes = FALSE)
```

You can switch build_vignettes to TRUE but it takes a lot longer to install and the vignettes are already available from CRAN. If installation fails, please let us know by filing an issue.
