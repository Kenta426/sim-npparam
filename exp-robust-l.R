# Simulation code corresponding to Figure 4 of the paper

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

get.L <- function(data.type){
  if(data.type == "iso1"){
    seq(3*0.7, 3*1.3, length.out=500)
  }
  else if(data.type == "iso2"){
    seq(0.7, 1.3, length.out=500)
  }
  else if(data.type == "cvx1"){
    seq(8*0.7, 8*1.3, length.out=500)
  }
  else if(data.type == "cvx2"){
    seq(0.7, 1.3, length.out=500)
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

simulate.one <- function(data, data.type="iso1", L=1){
  estimator <- get.estimator(data.type)
  fn.est <- estimator(data$x, data$y, L=L)
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
  L.list <- get.L(data.type)
  res <- lapply(L.list, function(l){simulate.one(data, data.type, l)})
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
save.dir <- file.path(base.dir, paste0("robust-L-", exp.type))

dir.create(save.dir)
out.df <- run.experiment(seed, sample.n, exp.type)

filename <- paste0("seed", seed, "N", sample.n, ".RData")
save(out.df, file=file.path(save.dir, filename))
