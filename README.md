# SignIT
Mutation Signatures in Individual Tumours

## Tutorial

Installation instructions are below. For a tutorial on SignIT, please see the [Vignette](https://eyzhao.github.io/SignIT/vignette.html).

## Installation

### Development Version

To install from GitHub, please follow these steps:

First make sure that you can install the rstan package and C++ toolchain by following [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). Also, be sure to install BSgenome from Bioconductor. If it is not installed, you can use the following commands in R:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('BSgenome')
```

You may also need to install specific cancer genomes to fit your need. For example:

```{r}
BiocManager("BSgenome.Hsapiens.UCSC.hg19")
```

Once rstan and BSgenome are successfully installed, you can install SignIT from GitHub using the devtools package by executing the following in R:

```{r}
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("eyzhao/SignIT", recompile = TRUE, build_vignettes = FALSE)
```

A common error in `install_github` is 

```
Downloading GitHub repo eyzhao/SignIT@master
from URL https://api.github.com/repos/eyzhao/SignIT/zipball/master
Installation failed: error in running command
```

If this happens, try executing `options(unzip = 'internal')` and try again.

When installing, if you get a C++ error along the lines of `C++14 standard requested but CXX14 is not defined`, try editing `~/.R/Makevars` to read

```
CXX14 = g++ -std=c++1y
CXX14FLAGS = -O3 -Wno-unused-variable -Wno-unused-function
```

### Setting up R and Necessary Dependencies

A convenient way to set up a standalone installation of R with SignIT is to use an Anaconda environment. In a Bash shell, the following commands will check out such an environment and install it locally.

```{bash}
git clone https://github.com/eyzhao/bio-pipeline-dependencies.git
cd bio-pipeline-dependencies
git checkout tags/SignIT-paper-dependencies
make
```

After this is complete, run the following command to activate the environment and run R.

```{bash}
source miniconda3/bin/activate dependencies
R
```
