## Packages
library(tidyr)
library(gldrm)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp("src/Cpp_functions.cpp")
source("src/R_functions.R")
