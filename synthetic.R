generate.fn <- function(sample.size = 100, sigma=0.1, freq=4){
  x <- runif(sample.size, 0, 1)
  y <- sin((x*2-1)*freq) + rnorm(sample.size, 0, sigma)
  idx <- order(x)
  list(x=x[idx], y=y[idx], fn=function(t){sin((t*2-1)*freq)})
}

lipschitz.fn <- function(sample.size = 100, sigma=0.1){
  x <- runif(sample.size, 0, 1)
  y <- rep(0, sample.size)
  fn <- function(x){
    sapply(x, function(t){
      if(t < 1/3){return(1 - 3*t)}
      else if(1/3 <= t && t < 2/3) return(-1+3*t)
      else{
        3 - 3*t
      }
    })
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  idx <- order(x)
  list(x=x[idx], y=y[idx], fn=fn)
}

adaptive.iso <- function(sample.size, L=1, m=3, sigma=0.1){
  x <- runif(sample.size, 0, 1)
  step <- ceiling((x*m)%%m)
  step <- (step > 0 & step <= m) * step + (step<=0) * 1 + (step > m) * m
  y <- step + L*x + rnorm(sample.size, 0, sigma)
  idx <- order(x)
  list(x=x[idx], y=y[idx],
       fn=function(t){
         step <- ceiling((t*m)%%m)
         step <- (step > 0 & step <= m) * step + (step<=0) * 1 + (step>m) * m
         step + L*t})
}

adaptive.cvx <- function(sample.size, L=1, m=3, sigma=0.1){
  x <- runif(sample.size, 0, 1)
  slopes <- (1:m)-2
  intercepts <- rep(0, m)
  step <- ceiling((x*m)%%m)
  step <- (step > 0 & step <= m) * step + (step<=0) * 1 + (step > m) * m

  if (m > 1){
    for (i in 2:m){
      intercepts[i] <- (slopes[i-1]-slopes[i])/m*(i-1) + intercepts[i-1]
    }
  }
  y <- slopes[step] * x + intercepts[step] + L*x^2 + rnorm(sample.size, 0, sigma)
  idx <- order(x)
  list(x=x[idx], y=y[idx],
       fn=function(t){
         step <- ceiling((t*m)%%m)
         step <- (step > 0 & step <= m) * step + (step<=0) * 1 + (step > m) * m
         slopes[step] * t + L * t^2 + intercepts[step]})
}

additive.sine <- function(sample.size, sigma=0.1){
  x <- matrix(runif(sample.size*5, 0, 1), ncol=5)
  freq.2 <- generate.fn(sample.size = 1, sigma=sigma, freq=2)
  freq.4 <- generate.fn(sample.size = 1, sigma=sigma, freq=4)
  freq.6 <- generate.fn(sample.size = 1, sigma=sigma, freq=6)
  freq.8 <- generate.fn(sample.size = 1, sigma=sigma, freq=8)
  freq.10 <- generate.fn(sample.size = 1, sigma=sigma, freq=10)

  fn <- function(t){
    freq.2$fn(t[, 1]) + freq.4$fn(t[, 2]) + freq.6$fn(t[, 3]) +
      freq.8$fn(t[, 4]) + freq.10$fn(t[, 5])
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  return(list(x=x, y=y,fn=fn))
}

additive.adaptive <- function(sample.size, sigma=0.1){
  x <- matrix(runif(sample.size*5, 0, 1), ncol=5)
  freq.2 <- adaptive.iso(sample.size = 1, sigma=sigma, L=0, m=1)
  freq.6 <- adaptive.iso(sample.size = 1, sigma=sigma,L=3, m=3)
  freq.8 <- adaptive.iso(sample.size = 1, sigma=sigma, L=1, m=3)

  fn <- function(t){
    freq.2$fn(t[, 1]) + freq.6$fn(1-t[, 2]) + freq.6$fn(t[, 3]) +
      freq.8$fn(1-t[, 4]) + freq.8$fn(t[, 5])
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  return(list(x=x, y=y,fn=fn))
}

additive.lipscitz <- function(sample.size, sigma=0.1){
  x <- matrix(runif(sample.size*5, 0, 1), ncol=5)
  d1 <- lipschitz.fn(1)
  fn <- function(t){
    d1$fn(t[,1]) - d1$fn(t[,2]) + abs(t[,3]) - abs(t[,4]) + 1
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  return(list(x=x, y=y,fn=fn))
}


additive.debug <- function(sample.size, sigma=0.1){
  x <- matrix(runif(sample.size*2, 0, 1), ncol=2)
  d1 <- lipschitz.fn(1)
  fn <- function(t){
    d1$fn(t[,1]) - d1$fn(t[,2])
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  return(list(x=x, y=y,fn=fn))
}

additive.debug.adaptive <- function(sample.size, sigma=0.1){
  x <- matrix(runif(sample.size*2, 0, 1), ncol=2)
  d1 <- adaptive.iso(1)
  fn <- function(t){
    d1$fn(t[,1]) - d1$fn(t[,2])
  }
  y <- fn(x) + rnorm(sample.size, 0, sigma)
  return(list(x=x, y=y,fn=fn))
}


# data <- adaptive.iso(1000)
# integrate(function(t){data$fn(t)^2},0,1)$value
# sqrt(integrate(function(t){data$fn(t)^2},0,1)$value/100)
# 
# data <- adaptive.iso(1000, sigma=0.3, L=1)
# plot(data$x, data$y)
# 
# data <- adaptive.cvx(1000, m=1, sigma=0.02, L=1)
# sqrt(integrate(function(t){data$fn(t)^2},0,1)$value/100)
# plot(data$x, data$y)
# lines(data$x, data$fn(data$x), col="red")
# data <- adaptive.cvx(1000, m=5, sigma=0.05, L=1)
# integrate(function(t){data$fn(t)^2},0,1)$value/50
# sqrt(integrate(function(t){data$fn(t)^2},0,1)$value/100)
# 
# plot(data$x, data$y)
# lines(data$x, data$fn(data$x), col="red")

