# Simulation code corresponding to Figure 5 of the paper

library(devtools)
setwd("../npparam/")
load_all()
library(npparam)
source("../sim-npparam/synthetic.R")
source("../sim-npparam/estimators.R")


get.data <- function(data.type){
  if(data.type == "iso1"){
    lipschitz.fn
  }
  else if(data.type == "iso2"){
    adaptive.iso
  }
  else if(data.type == "cvx1"){
    generate.fn
  }
  else if(data.type == "cvx2"){
    adaptive.cvx
  }
}

get.estimator <- function(data.type){
  if((data.type == "iso1") || (data.type == "iso2")){
    iso.lin.est
  }
  else if((data.type == "cvx1") || (data.type == "cvx2")){
    cvx.quad.est
  }
}

simulate.one <- function(data, data.type="iso1"){
  estimator <- get.estimator(data.type)
  fn.est <- estimator(data$x, data$y)
  out <- list("fn.est"=fn.est)
  out
}

run.experiment <- function(seed, n, data.type){
  set.seed(seed)
  test.x <- runif(1e4)
  data.fn <- get.data(data.type)
  data <- data.fn(sample.size = 10, sigma = 0.1)
  
  test.y <- data$fn(test.x)
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("n", "seed", "risk", "estimator", "L")
  colnames(df) <- x
  
  # expose data
  data <- data.fn(sample.size = n, sigma = 0.1)
  res <- replicate(50, simulate.one(data, data.type), simplify = FALSE)
  
  est.risk <- unlist(lapply(res, function(m){
    mean((predict(m$fn.est, test.x)-test.y)^2)}))
  est.L <- unlist(lapply(res, function(m){m$fn.est$l.value}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="iso+lin",
                       L=est.L)
  df <- rbind(df, out.df)
  return(df)
}

# call Rscript cross_validation.R 100 500
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1]) # seed
cat("Base seed = ", seed, "\n")

sample.n <- as.numeric(args[2]) # seed
cat("Sample size = ", sample.n, "\n")

exp.type <- as.numeric(args[3]) # seed
cat("Experiment type = ", exp.type, "\n")
base.dir <- "../sim-npparam/results"
save.dir <- file.path(base.dir, paste0("robust-cv-", exp.type))

dir.create(save.dir)
out.df <- run.experiment(seed, sample.n, exp.type)

filename <- paste0("seed", seed, "N", sample.n, ".RData")
save(out.df, file=file.path(save.dir, filename))
