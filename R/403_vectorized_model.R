#' Vectorized Tensor Regression with Gaussian Graphical Model
#'
#' @param input A list with at least the elements *Y* and *x*. *Y* should be a list of *G* arrays of dimension $D + 2$ in which the last two dimension margins represent the time and the subject, respectively.
#' @param n.iter The number of iterations that the MCMC should run
#' @param n.burn The number of iterations at the beginning of the sampler that should be discarded as a "burn-in"
#' @param hyperparameters A data frame with one observation and eight variables: *a.tau*, *b.tau*, *a.lambda*, *b.lambda*, *a.zeta*, *b.zeta*, *a.sig*, and *b.sig*. Each should specify values to use for the hyperparameters in the model. If `NULL`, then the default values used are a.tau = D - 1 b.tau = 1^((1 / D) - 1), a.lambda = 3, b.lambda = 3^(1/(2*D)), a.zeta = 1, b.zeta = 0.01, a.sig = 1, b.sig = -log(0.95).
#'
#' @return A list with 9 elements, *B*, *W*, *lambda*, *tau*, *d*, *zeta*, *Sig*, *sig2y*, and *llik*. *B* is a list of length *n.iter* - *n.burn*, each with $G$ sublists. This will be improved in the future. *llik* is a vector of the log-likelihood values at each MCMC iteration.
#' @export
#'
#' @examples
vec_Y_x_with_GGM <- function(input, n.iter, n.burn,hyperparameters = NULL){

  if(n.burn >= n.iter) return("You should burn fewer than your total number of simulations!")
  require(tidyverse)

  n <- tail(dim(input$Y[[1]]),1) # Number of subjects
  if(is.vector(input$x)) input$x <- tcrossprod(input$x,rep(1,n))
  G <- length(input$Y) # Number of regions of interest
  TT <- dim(input$x)[1] # Number of time steps
  p <- sapply(input$Y,function(each_region) head(dim(each_region),-2),simplify=FALSE)
  D <- length(dim(input$Y[[1]])) - 2

  # > Vectorize the response ----
  input$Y <- sapply(input$Y,function(each_region) apply(each_region,(D+1):(D+2), identity))

  # Load libraries ----
  require(GIGrvg)
  require(mvtnorm)
  require(doParallel)
  require(abind)
  require(dlm)
  require(statmod)
  require(truncnorm)

  # Set hyperparameters ----
  if(is.null(hyperparameters)) {
    hyperparameters <- data.frame(
      a.tau = D - 1,
      b.tau = 1^((1 / D) - 1), # These values are from Guhaniyogi et al [2017]
      a.lambda = 3,
      b.lambda = 3^(1/(2*D)), # These values are from Guhaniyogi et al [2017]
      a.zeta = 1,
      b.zeta = 0.01, # These values are from Wang [2012]
      a.sig = 1,
      b.sig = -log(0.95) # These values are from Guhaniyogi et al [2017]
    )
  }

  a.tau <- hyperparameters$a.tau
  b.tau <- hyperparameters$b.tau
  a.lambda <- hyperparameters$a.lambda
  b.lambda <- hyperparameters$b.lambda
  a.zeta <- hyperparameters$a.zeta
  b.zeta <- hyperparameters$b.zeta
  a.sig <- hyperparameters$a.sig
  b.sig <- hyperparameters$b.sig

  # > Set up Parallelization ----
  num <- min(G, detectCores() - 1)

  cl <- makeCluster(num)
  registerDoParallel(cl)
  clusterEvalQ(cl = cl, library(mvtnorm))
  clusterEvalQ(cl = cl, library(GIGrvg))
  clusterEvalQ(cl = cl, library(doParallel))
  clusterEvalQ(cl = cl, library(abind))
  clusterEvalQ(cl = cl, library(truncnorm))

  # > Functions ----

  # > Set Storage ----
  results <-
    list(
      B = list(),
      W = list(),
      lambda = list(),
      tau = matrix(NA, G, n.iter),
      d = array(NA, dim = c(G, n, n.iter)),
      zeta = numeric(n.iter),
      Sig = array(NA, dim = c(G, G, n.iter)),
      sig2y = numeric(n.iter),
      llik = numeric(n.iter)
    )

  # > Set Initials ----
  d <- t(sapply(input$Y,function(each_region){
    apply(each_region,3,mean)
  }))

  betas <-  parallel::parSapplyLB(cl,seq(G),function(each_region,input){
    apply(input$Y[[each_region]], 1, function(each_voxel) {
      c(RcppArmadillo::fastLmPure(as.matrix(c(input$x)),c(each_voxel))$coefficients)
    })
  },input = input)

  Sig <- diag(1, G)  ## Suggested by Wang
  Upsilon <- matrix(1, G, G) ## Suggested by Wang
  diag(Upsilon) <- 0 ## Suggested by Wang
  zeta <- 3  ## Value used in Wang et al. [2012]
  tau <- rep(1, G) # Assuming unit value
  V <- lapply(p,prod)
  lambda <- rep(1,G)
  W <-
    sapply(betas, function(each_region) {
      rep(1,length(each_region))
    },simplify = FALSE) # List of length 5, each element a vector of length V_g of 1s
  sig2y <- 1 # Observational variance estimate

  ## Begin MCMC
  beginning_of_sampler <- proc.time()[3]
  # pb <- txtProgressBar(min=1,max=n.iter,style=3)
  for(s in 1:n.iter){
    # Do everything region-specific in parallel
    params <- foreach(each_region = seq(G)) %do% {
      # >> Draw tau ----
      xi <- a.tau - V[[each_region]]/2
      chi <- c(crossprod(betas[[each_region]],diag(1/W[[each_region]])) %*%
                 betas[[each_region]])
      tau.g <- rgig(1, xi, chi, 2 * b.tau)

      # >> Draw lambda ----
      lambda.g <- rgamma(1,a.lambda + length(betas[[each_region]]), b.lambda + sum(abs(betas[[each_region]]))/sqrt(tau.g))

      # >> Draw omega ----
      omega.g <- sapply(seq(length(betas[[each_region]])),function(each_element){
        GIGrvg::rgig(1,1/2,betas[[each_region]][each_element]^2 / tau.g,lambda.g^2)
      })

      # >> Draw beta ----
      beta_cov <- chol2inv(chol(diag(1/omega.g)/tau.g + diag(sum(input$x^2),length(omega.g))/sig2y))
      beta_mean <- c(beta_cov %*%
        apply(
          sapply(seq(dim(input$Y[[each_region]])[2]),function(each_time){
            sapply(seq(dim(input$Y[[each_region]])[3]),function(each_subject){
              (input$Y[[each_region]][,each_time,each_subject] -
                 d[each_region,each_subject])*
                input$x[each_time,each_subject]
            },simplify = "array")
          },simplify = "array"),
        1,sum) / sig2y)

      beta.g <- c(beta_mean + rnorm(length(beta_mean)) %*% chol(beta_cov))

      # Put all of the region-specific draws together
      list(
        tau = tau.g,
        omega = omega.g,
        B = beta.g,
        lambda = lambda.g
      )
    } # End dopar

    # Draw d ----

    y_hat <- mapply(function(y,parm){
      y - parm$B %o% input$x
    },y=input$Y,parm=params)

    inv_d_covar <- Sig + (TT*diag(sapply(p,prod)))^2 / sig2y
    #  #####################################################
    #  # Inverse Wishart prior on Sigma
    # (Or when the prior is applied to the covariance rather than the precision)
    # inv_d_covar <- chol2inv(chol(Sig)) + (TT*diag(sapply(p,prod)))^2 / sig2y
    #  #####################################################
    d_covar <- chol2inv(chol(inv_d_covar))
    # d_covar <- solve((inv_d_covar + t(inv_d_covar))/2)
    d_mean <- d_covar%*% t(sapply(y_hat, function(yh){apply(yh,3,sum)}))

    d <- apply(d_mean,2,function(dm){
      dm + rnorm(length(dm)) %*% chol((d_covar +
                                         t(d_covar))/2)
    })

    # >> Draw Sig ----
    matd <- t(d)
    S <- crossprod(matd,matd)
    # S <- chol2inv(chol(S)) # Will this work if I'm looking at a covariance rather than a precision?
    for(g in 1:G){
      delta <- rgamma(1,n/2 + 1,(S[g,g] + zeta)/2)
      S11inv <- chol2inv(chol(Sig[-g,-g]))
      # S11inv <- solve(Sig[-g,-g])
      varalpha <- chol2inv(chol((S[g,g] + zeta)*S11inv + diag(as.matrix(1/Upsilon[-g,g]))))
      # varalpha <- solve((S[g,g] + zeta)*S11inv + diag(as.matrix(1/Upsilon[-g,g])))
      meanalpha <- c(-varalpha%*%S[g,-g])
      s_alpha <- meanalpha + rnorm(G-1) %*% chol(varalpha)
      Sig[g,-g] <- s_alpha
      Sig[-g,g] <- s_alpha
      Sig[g,g] <- delta + (s_alpha)%*%S11inv%*%t(s_alpha)
      if(g < G){
        for(gg in (g+1):G){
          Upsilon[gg,g] <- Upsilon[g,gg] <- 1/rinvgauss(1,sqrt(zeta^2 / Sig[g,gg]^2),zeta^2)
        }
      }
    }
    zeta <- rgamma(1,a.zeta + G*(G+1)/2, b.zeta + sum(abs(Sig))/2)


    # Draw sig2y ----
    lik_sse <- mapply(function(parm,yy,dd){
       sum((yy - parm$B %o% input$x - array(1,dim = c(length(parm$B),TT)) %o% dd)^2)
    },yy = input$Y,parm=params,dd = split(d,row(d))) %>% unlist %>% sum

    sig2y <- 1/rgamma(1,
                      a.sig + (n * TT * sum(unlist(V))) / 2,
                      b.sig + lik_sse/2)

    # Grab the log-likelihood ----
    results$llik[s] <- -0.5*log(2*pi*sig2y)*n*TT*G*sum(unlist(V)) - 0.5/sig2y * sum(lik_sse)

    # Transform the way that the draws are stored
    betas <- lapply(params, function(z) {
      z$B
    })
    W <- lapply(params, function(z) {
      z$omega
    })
    lambda <- lapply(params, function(z) {
      z$lambda
    })
    tau <- sapply(params, function(z) {
      z$tau
    })

    # Store the results for each iteration
    results$B[[s]] <- betas
    results$W[[s]] <- W
    results$lambda[[s]] <- lambda
    results$tau[, s] <- tau
    results$d[, , s] <- d
    results$Sig[, , s] <- Sig
    results$zeta[s] <- zeta
    results$sig2y[s] <- sig2y

    # setTxtProgressBar(pb, s)
    if(s %% ceiling(n.iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Vectorised model -",
          " Iteration # ",s," of ",n.iter,
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - beginning_of_sampler, "  seconds #####\n",
          "##### Estimated time remaining: ", ((proc.time()[3] - beginning_of_sampler)/s)*(n.iter - s)," seconds #####\n"
        )
      )
    }

  } # End sampler

  if(n.burn > 0){
    results$B <- results$B[-(1:n.burn)]
    results$W <- results$W[-(1:n.burn)]
    results$lambda <- results$lambda[-(1:n.burn)]
    results$tau <- results$tau[, -(1:n.burn)]
    results$d <- results$d[, , -(1:n.burn)]
    results$Sig <- results$Sig[, , -(1:n.burn)]
    results$zeta <- results$zeta[-(1:n.burn)]
    results$sig2y <- results$sig2y[-(1:n.burn)]
  }

  results
}
