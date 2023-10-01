# Simulation code corresponding to Figure 2 of the paper (Scenario 2)

library(devtools)
setwd("../npparam/")
load_all()
library(npparam)
source("../sim-npparam/synthetic.R")
source("../sim-npparam/estimators.R")


simulate.one <- function(n){
  data <- adaptive.iso(sample.size = n, sigma = 0.1)
  cat("Fitting isotonic for sample size", n, "\n")
  iso.est <- iso.lin.est(data$x, data$y)
  cat("Fitting isotonic (CV) for sample size", n, "\n")
  iso.est.cv <- iso.lin.est(data$x, data$y, run.optim = FALSE)
  cat("Fitting sieve for sample size", n, "\n")
  sieve.est <- .fit.sieve(data$x, data$y)
  cat("Fitting GBM for sample size", n, "\n")
  gbm.est <- .fit.gbm(data$x, data$y)
  if(n > 2000){
    cat("Skipping KRR for sample size", n, "\n")
    krr.est <- NULL
    cat("Skipping random forest for sample size", n, "\n")
    rf.est <- NULL
  }else{
    cat("Fitting KRR for sample size", n, "\n")
    krr.est <- .fit.krr(data$x, data$y)
    cat("Fitting random forest for sample size", n, "\n")
    rf.est <- .fit.random.forest(data$x, data$y)
  }
  out <- list("iso.est"=iso.est,
              "iso.est.cv"=iso.est.cv,
              "sieve.est"=sieve.est,
              "krr.est"=krr.est,
              "rf.est"=rf.est,
              "gbm.est"=gbm.est)
  out
}

run.experiment <- function(seed, n){
  set.seed(seed)
  test.x <- runif(1e4)
  data <- adaptive.iso(sample.size = 10, sigma = 0.1)

  test.y <- data$fn(test.x)
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("n", "seed", "risk", "estimator", "L")
  colnames(df) <- x


  res <- replicate(10, simulate.one(n), simplify = FALSE)

  est.risk <- unlist(lapply(res, function(m){
    mean((predict(m$iso.est, test.x)-test.y)^2)}))
  est.L <- unlist(lapply(res, function(m){m$iso.est$l.value}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="iso+lin",
                       L=est.L)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((predict(m$iso.est.cv, test.x)-test.y)^2)}))
  est.L <- unlist(lapply(res, function(m){m$iso.est.cv$l.value}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="iso+lin cv",
                       L=est.L)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((m$sieve.est(test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="sieve",
                       L=0)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((m$gbm.est(test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="gbm",
                       L=0)
  df <- rbind(df, out.df)

  if (n <= 2000){
    est.risk <- unlist(lapply(res, function(m){
      mean((m$krr.est(test.x)-test.y)^2)}))
    out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="krr",
                         L=0)
    df <- rbind(df, out.df)

    est.risk <- unlist(lapply(res, function(m){
      mean((m$rf.est(test.x)-test.y)^2)}))
    out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="rf",
                         L=0)
    df <- rbind(df, out.df)
  }

  return(df)
}

# call Rscript cross_validation.R 100 500
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1]) # seed
cat("Base seed = ", seed, "\n")

sample.n <- as.numeric(args[2]) # seed
cat("Sample size = ", sample.n, "\n")
base.dir <- "../sim-npparam/results"
save.dir <- file.path(base.dir, "iso-exp2")

dir.create(save.dir)
out.df <- run.experiment(seed, sample.n)

filename <- paste0("seed", seed, "N", sample.n, ".RData")
save(out.df, file=file.path(save.dir, filename))
