# collection of estimators
library(Sieve)
library(randomForest)
library(gbm)
library(mgcv)

.fit.sieve <- function(x, y){
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]

  type <- 'cosine'; basisN <- 50
  sieve.model <- sieve_preprocess(X = x.train, basisN = basisN, type = type)
  sieve.fit<- sieve_solver(model = sieve.model, Y = y.train)
  sieve.prediction <- sieve_predict(model = sieve.fit,
                                    testX = x.val,
                                    testY = y.val)

  best_lambda_index <- which.min(sieve.prediction$MSE) #selected lambda

  predict.fn <- function(test.x){
    sieve.prediction <- sieve_predict(model = sieve.fit,
                                      testX = test.x)
    sieve.prediction$predictY[, best_lambda_index]
  }
  predict.fn
}

.fit.krr <- function(x, y){
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]

  type <- 'sobolev1' ##this kernel is K(x,y) = 1+min(x,y), you can find it in Wainwright's book
  krr_model <- Sieve:::KRR_preprocess(X = as.matrix(x.train, ncol = 1), type)

  eval.lambda <- function(lambda){
    estimate_beta <- Sieve:::KRR_cal_beta_C(U = krr_model$U, s = krr_model$s,
                                            Y = y.train, lambda = lambda)
    predicted_value_krr <- Sieve:::KRR_predict_C(trainX = as.matrix(x.train, ncol = 1),
                                                 testX = as.matrix(x.val, ncol = 1),
                                                 type = type, beta_hat = estimate_beta,
                                                 kernel_para = -1)
    mean((y.val-predicted_value_krr)^2)
  }
  val.risk <- sapply(krr_model$lambdas, eval.lambda)
  optimal.lambda <- krr_model$lambdas[which.min(val.risk)]

  function(test.x){
    estimate_beta <- Sieve:::KRR_cal_beta_C(U = krr_model$U,
                                            s = krr_model$s,
                                            Y = y.train,
                                            lambda = optimal.lambda)
    Sieve:::KRR_predict_C(trainX = as.matrix(x.train, ncol = 1),
                          testX = as.matrix(test.x, ncol = 1),
                          type = type, beta_hat = estimate_beta,
                          kernel_para = -1)
  }
}

.fit.random.forest <- function(x, y){
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]

  eval.tree <- function(tree){
    rf <- randomForest(Y ~ .,
                       data=data.frame(x=x.train, Y=y.train),
                       ntree = tree)
    predictY <-  predict(rf, newdata=data.frame(x=x.val))
    mean((y.val-predictY)^2)
  }
  ntrees <- c(50, 1e2, 5e2, 1e3, 5e3)
  val.risk <- sapply(ntrees, eval.tree)
  optimal.ntree <- ntrees[which.min(val.risk)]
  function(test.x){
    rf <- randomForest(Y ~ .,
                       data=data.frame(x=x.train, Y=y.train),
                       ntree = optimal.ntree)
    predict(rf, newdata=data.frame(x=test.x))
  }
}

.fit.gbm <- function(x, y){
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]

  best.risk <- Inf
  opt.shrinkage <- 0; opt.tree <- 0; opt.depth <- 0
  for(shrinkage in c(1e-2)){
    for(n.tree in c(1e2, 1e3, 2e3, 4e3, 8e3)){
      for(interaction.depth in c(2,5)){
        bm <- gbm(Y ~ ., data=data.frame(x=x.train, Y=y.train),
                  distribution = "gaussian",
                  interaction.depth = interaction.depth, n.trees = n.tree,
                  shrinkage = shrinkage, bag.fraction = 0.5)
        predictY <- predict(bm, newdata=data.frame(x=x.val))
        val.risk <- mean((predictY-y.val)^2)
        if (val.risk < best.risk) {
          best.risk <- val.risk;
          opt.shrinkage <- shrinkage; opt.tree <- n.tree;
          opt.depth <- interaction.depth
        }
      }
    }
  }

  function(test.x){
    bm <- gbm(Y ~ ., data=data.frame(x=x.train, Y=y.train),
              distribution = "gaussian",
              interaction.depth = opt.depth, n.trees = opt.tree,
              shrinkage = opt.shrinkage, bag.fraction = 0.5)
    predict(bm, newdata=data.frame(x=test.x))
  }
}


.fit.additive.sieve <- function(x, y){
  n <- dim(x)[1]; d <- dim(x)[2]; learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id,]; x.val <- x[-train.id,]
  y.train <- y[train.id]; y.val <- y[-train.id]

  type <- 'cosine'; basisN <- as.integer(2*d*sqrt(n))
  sieve.model <- sieve_preprocess(X = x.train, basisN = basisN, type = type,
                                  interaction_order = 1)
  sieve.fit<- sieve_solver(model = sieve.model, Y = y.train)
  sieve.prediction <- sieve_predict(model = sieve.fit,
                                    testX = x.val,
                                    testY = y.val)

  best_lambda_index <- which.min(sieve.prediction$MSE) #selected lambda

  predict.fn <- function(test.x){
    sieve.prediction <- sieve_predict(model = sieve.fit,
                                      testX = test.x)
    sieve.prediction$predictY[, best_lambda_index]
  }
  predict.fn
}

.fit.additive.krr <- function(x, y){
  n <- dim(x)[1]; d <- dim(x)[2]; learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id,]; x.val <- x[-train.id,]
  y.train <- y[train.id]; y.val <- y[-train.id]

  type <- 'sobolev1' ##this kernel is K(x,y) = 1+min(x,y), you can find it in Wainwright's book
  krr_model <- Sieve:::KRR_preprocess(X = x.train, type)

  eval.lambda <- function(lambda){
    estimate_beta <- Sieve:::KRR_cal_beta_C(U = krr_model$U, s = krr_model$s,
                                            Y = y.train, lambda = lambda)
    predicted_value_krr <- Sieve:::KRR_predict_C(trainX = x.train,
                                                 testX = x.val,
                                                 type = type, beta_hat = estimate_beta,
                                                 kernel_para = -1)
    mean((y.val-predicted_value_krr)^2)
  }
  val.risk <- sapply(krr_model$lambdas, eval.lambda)
  optimal.lambda <- krr_model$lambdas[which.min(val.risk)]

  function(test.x){
    estimate_beta <- Sieve:::KRR_cal_beta_C(U = krr_model$U,
                                            s = krr_model$s,
                                            Y = y.train,
                                            lambda = optimal.lambda)
    Sieve:::KRR_predict_C(trainX = x.train,
                          testX = test.x,
                          type = type, beta_hat = estimate_beta,
                          kernel_para = -1)
  }
}

.fit.additive.gbm <- function(x, y){
  n <- dim(x)[1]; d <- dim(x)[2]; learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id,]; x.val <- x[-train.id,]
  y.train <- y[train.id]; y.val <- y[-train.id]

  best.risk <- Inf
  opt.shrinkage <- 0; opt.tree <- 0; opt.depth <- 0
  for(shrinkage in c(1e-2)){
    for(n.tree in c(1e2, 1e3, 2e3, 4e3, 8e3)){
      for(interaction.depth in c(1)){
        bm <- gbm(Y ~ ., data=data.frame(x=x.train, Y=y.train),
                  distribution = "gaussian",
                  interaction.depth = interaction.depth, n.trees = n.tree,
                  shrinkage = shrinkage, bag.fraction = 0.5)
        predictY <- predict(bm, newdata=data.frame(x=x.val))
        val.risk <- mean((predictY-y.val)^2)
        if (val.risk < best.risk) {
          best.risk <- val.risk;
          opt.shrinkage <- shrinkage; opt.tree <- n.tree;
          opt.depth <- interaction.depth
        }
      }
    }
  }

  function(test.x){
    bm <- gbm(Y ~ ., data=data.frame(x=x.train, Y=y.train),
              distribution = "gaussian",
              interaction.depth = opt.depth, n.trees = opt.tree,
              shrinkage = opt.shrinkage, bag.fraction = 0.5)
    predict(bm, newdata=data.frame(x=test.x))
  }
}

.fit.gam <- function(x, y){
  n <- length(y)
  df <- data.frame(y=y, x)
  k <- 30
  gam.est <- gam(y~s(X1, k=k)+s(X2, k=k)+s(X3, k=k)+s(X4, k=k)+s(X5, k=k), data = df)
  function(test.x){
    predict.gam(gam.est, newdata = data.frame(test.x))
  }
}

.fit.gam.2d <- function(x, y){
  df <- data.frame(y=y, x)
  k <- 30
  gam.est <- gam(y~s(X1, k=k)+s(X2, k=k), data = df)
  function(test.x){
    predict.gam(gam.est, newdata = data.frame(test.x))
  }
}
