## Packages
library(tidyr)
library(gldrm)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp("src/Cpp_functions.cpp")
Rcpp::sourceCpp("src/theta.cpp")
source("src/R_functions.R")

