# Simulation code corresponding to Figure 6 of the paper

library(devtools)
setwd("../npparam/")
load_all()
library(npparam)
source("../sim-npparam/synthetic.R")
source("../sim-npparam/estimators.R")
simulate.one <- function(m, snr=10){
  data <-  adaptive.iso(sample.size = 5000, m=m, sigma = 0.1)
  sigma <- 0.1
  cat("Fitting isotonic for sd", sigma, "\n")
  iso.est <- iso.lin.est(data$x, data$y)
  cat("Fitting sieve for sd", sigma, "\n")
  sieve.est <- .fit.sieve(data$x, data$y)
  cat("Fitting GBM for sd", sigma, "\n")
  gbm.est <- .fit.gbm(data$x, data$y)

  out <- list("iso.est"=iso.est,
              "sieve.est"=sieve.est,
              "gbm.est"=gbm.est)
  out
}

run.experiment <- function(seed, m, snr=10){
  set.seed(seed)
  test.x <- runif(1e4)
  data <-  adaptive.iso(sample.size = 10, m=m, sigma = 0)
  test.y <- data$fn(test.x)
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  x <- c("m","snr", "seed", "risk", "estimator", "L")
  colnames(df) <- x
  
  
  res <- replicate(30, simulate.one(m, snr), simplify = FALSE)
  
  est.risk <- unlist(lapply(res, function(m){
    mean((predict(m$iso.est, test.x)-test.y)^2)}))
  est.L <- unlist(lapply(res, function(m){m$iso.est$l.value}))
  out.df <- data.frame(m=m, snr=snr,  seed=seed, risk=est.risk, estimator="iso+lin",
                       L=est.L)
  df <- rbind(df, out.df)
  
  est.risk <- unlist(lapply(res, function(m){
    mean((m$sieve.est(test.x)-test.y)^2)}))
  out.df <- data.frame(m=m, snr=snr,  seed=seed, risk=est.risk, estimator="sieve",
                       L=0)
  df <- rbind(df, out.df)
  
  est.risk <- unlist(lapply(res, function(m){
    mean((m$gbm.est(test.x)-test.y)^2)}))
  out.df <- data.frame(m=m, snr=snr,  seed=seed, risk=est.risk, estimator="gbm",
                       L=0)
  df <- rbind(df, out.df)

  return(df)
}


# call Rscript cross_validation.R 100 500
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1]) # seed
cat("Base seed = ", seed, "\n")

m.piece <- as.numeric(args[2]) # m
cat("M = ", m.piece, "\n")
base.dir <- "../sim-npparam/results"
save.dir <- file.path(base.dir, "iso-exp1-m-piece")

dir.create(save.dir)
out.df <- run.experiment(seed, m.piece)

filename <- paste0("seed", seed, "N", m.piece, ".RData")
save(out.df, file=file.path(save.dir, filename))
