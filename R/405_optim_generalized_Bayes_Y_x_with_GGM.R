#' Generalized-Dimension Bayesian Tensor Response Regression with Gaussian Graphical Model
#'
#' @param input A list with at least the elements *Y* and *x*. *Y* should be a list of *G* arrays of dimension $D + 2$ in which the last two dimension margins represent the time and the subject, respectively.
#' @param n.iter The number of iterations that the MCMC should run
#' @param n.burn The number of iterations at the beginning of the sampler that should be discarded as a "burn-in"
#' @param ranks The number of ranks in the model. See Spencer et al. [2019] for details.
#' @param hyperparameters A data frame with one observation and eight variables: *a.tau*, *b.tau*, *a.lambda*, *b.lambda*, *a.zeta*, *b.zeta*, *a.sig*, and *b.sig*. Each should specify values to use for the hyperparameters in the model. If `NULL`, then the default values of a.tau = D - 1, b.tau = ranks^((1 / D) - 1), a.lambda = 3, b.lambda = 3^(1/(2*D)), a.zeta = 1, b.zeta = 0.01, a.sig = 1, b.sig = -log(0.95) are used.
#'
#' @return A list with 12 elements, *B*, *W*, *lambda*, *Phi*, *tau*, *d*, *zeta*, *Sig*, *sig2y*, *alpha*, *accept*, and *llik*. *B* is a list of length *n.iter* - *n.burn*, each with $G$ sublists. This will be improved in the future. *llik* is a vector of the log-likelihood values at each MCMC iteration.
#' @export
#'
#' @examples
#' BTRR_model_results <- BTR_Y_x_with_GGM(simulated_fmri_data,10,0,1)
#'
BTR_Y_x_with_GGM <- function(input, n.iter, n.burn, ranks,
                             hyperparameters = NULL, save_after = NULL,
                             save_llik = TRUE, results_list = NULL) {
  require(tidyverse)
  #  > Read in the data appropriately ----
  if(class(input$Y) == "array"){
    G <- 1
    lowest_dim_margin <- min(head(dim(input$Y),-2))
    n <- tail(dim(input$Y),1)
    p <- list(head(dim(input$Y),-2))
    D <- length(dim(input$Y)) - 2
  }else{
    lowest_dim_margin <- min(sapply(input$Y,function(each_region) head(dim(each_region),-2)))
    n <- tail(dim(input$Y[[1]]),1) # Number of subjects
    G <- length(input$Y) # Number of regions of interest
    p <- sapply(input$Y,function(each_region) head(dim(each_region),-2),simplify=FALSE)
    D <- length(dim(input$Y[[1]])) - 2
  }
  if(lowest_dim_margin < ranks) stop("The rank of your model cannot be larger than your smallest tensor dimension. Try a lower model rank.")
  if(is.vector(input$x)) input$x <- tcrossprod(input$x,rep(1,n))
  TT <- dim(input$x)[1] # Number of time steps

  # > Load necessary packages ----
  require(GIGrvg)
  # require(mvnfast)
  require(doParallel)
  require(abind)
  require(dlm)
  require(statmod)
  require(truncnorm)
  require(Rcpp)
  require(RcppArmadillo)

  # > Functions ----
  # source("R/900_misc.R")

  # > Set up Parallelization ----
  if(G > 1){
    num <- min(G, detectCores() - 1)
    if(num > 1){
      cl <- makeCluster(num)
      registerDoParallel(cl)
      clusterEvalQ(cl = cl, library(GIGrvg))
      clusterEvalQ(cl = cl, library(doParallel))
      clusterEvalQ(cl = cl, library(abind))
      clusterEvalQ(cl = cl, library(truncnorm))
      clusterEvalQ(cl = cl, library(Rcpp))
      clusterEvalQ(cl = cl, library(RcppArmadillo))
      clusterEvalQ(cl = cl, library(tidyverse))
      # clusterEvalQ(cl = cl, source("R/900_misc.R"))
    }
  }

  # > Hyperparameters ----
  if(is.null(hyperparameters)) {
    hyperparameters <- data.frame(
      a.tau = D - 1,
      b.tau = ranks^((1 / D) - 1), # These values are from Guhaniyogi et al [2017]
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

  ################################################
  # These are the hyperparameters for an inverse Wishart prior for covariance of D
  # nu <- G
  # VV <- diag(n,G)
  ################################################

  # > Set Storage ----
  results <-
    list(
      B = list(),
      W = list(),
      lambda = list(),
      Phi = list(),
      tau = matrix(NA, G, n.iter),
      d = array(NA, dim = c(G, n, n.iter)),
      zeta = numeric(n.iter),
      Sig = array(NA, dim = c(G, G, n.iter)),
      llik = vector(mode = "numeric",length = n.iter),
      sig2y = numeric(n.iter),
      alpha = numeric(n.iter),
      accept = vector(mode = "numeric",length = G)
    )

  # > Set Initials ----
  # >> Continuing from other result files ----
  if(!is.null(results_list)){
    S <- length(results_list$B)
    d <- results_list$d[,,S]
    betas <- results_list$B[[S]]
    if(G > 1){
      Sig <- results_list$Sig[,,S]
      zeta <- results_list$zeta[S]
      Upsilon <- matrix(1, G, G) ## Suggested by Wang
      diag(Upsilon) <- 0 ## Suggested by Wang
      for(g in seq(G)){
        if(g < G){
          for(gg in (g+1):G){
            Upsilon[gg,g] <- Upsilon[g,gg] <- 1/rinvgauss(1,sqrt(zeta^2 / Sig[g,gg]^2),zeta^2)
          }
        }
      }
    }
    tau <- results_list$tau[,S]
    alpha.g <-
      seq(ranks ^ (-D), ranks ^ (-.1), length.out = 10) # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
    lambda <- results_list$lambda[[S]]
    W <- results_list$W[[S]]
    phi <- results_list$Phi[[S]]
    sig2y <- results_list$sig2y[S]
    if(ranks > 1){
      Xi <- sapply(phi,function(each_rank){
        out <- numeric(ranks - 1)
        out[1] <- each_rank[1]
        for(r in seq(2,ranks - 1)){
          out[r] <- each_rank[r] / (prod(head((1 - each_rank),r-1)))
        }
        return(out)
      # },simplify = FALSE)
      })
    }else{
      Xi <- matrix(1,1,G)
    }
    cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(ranks - 1),simplify = FALSE)
  }else{
    # >> Using best-guess initial values ----
    d <- matrix(0,G,n)
    if(D > 2){
      B <-  parallel::parSapplyLB(cl,seq(G),function(each_region,input,TT,n,D){
        apply(array(input$Y[[each_region]],dim=c(dim(input$Y[[each_region]])[seq(D)], prod(dim(input$Y[[each_region]])[(D+1):(D+2)]))), seq(D), function(each_voxel) {
          # lm(each_voxel ~ c(x))$coefficients[2]
          ifelse(sd(each_voxel) != 0,
                 c(RcppArmadillo::fastLmPure(as.matrix(c(input$x)),each_voxel)$coefficients), # This is much faster
                 0)
        })
      },input = input,TT=TT,n=n,D=D)
    }else{
      if(G > 1){
        B <-  sapply(seq(G),function(each_region,input){
          each_region_y <- input$Y[[each_region]]
          yyg <-
            array(c(each_region_y), dim = c(dim(each_region_y)[seq(D)], prod(dim(each_region_y)[(D+1):(D+2)])))
          apply(yyg, seq(D), function(each_voxel) {
            ifelse(sd(each_voxel) != 0,
                   c(RcppArmadillo::fastLmPure(as.matrix(c(input$x)),each_voxel)$coefficients),
                   0)

          })
        },input = input)
      }else{
        B <- apply(input$Y,seq(D),function(each_voxel){
          c(RcppArmadillo::fastLmPure(as.matrix(c(input$x)),c(each_voxel))$coefficients)
        })
      }
    }

    ## The initial values set using the MLE
    betas <- sapply(seq(G),function(each_region){
      sapply(seq(D),function(each_dim){
        if(G > 1) {
          svd_B <- svd(mode_k_matriz(B[[each_region]],each_dim))
        }else{
          svd_B <- svd(mode_k_matriz(B,each_dim))
        }
        sapply(seq(ranks),function(each_rank){
          (svd_B$d[each_rank])^(1/D) * svd_B$u[,each_rank]
        })
      },simplify = F)
    },simplify = FALSE)

    if(G > 1){
      y_prime <- mapply(function(yy,bb){
        activation_residuals <- yy - bb%o%input$x
        aggregated_activation_residuals <- apply(activation_residuals,length(dim(yy)),mean)
        return(aggregated_activation_residuals)
      },yy = input$Y,bb = B)
      Sig <- solve(crossprod(y_prime)/n)
      Sig_L <- diag(1,G*(G-1)/2)
      Upsilon <- matrix(1, G, G) ## Suggested by Wang
      diag(Upsilon) <- 0 ## Suggested by Wang
      zeta <- 3  ## Value used in Wang et al. [2012]
    }
    tau <- rep(1, G) # Assuming unit value
    alpha.g <-
      seq(ranks ^ (-D), ranks ^ (-.1), length.out = 10) # This is the grid of alpha values suggested by Guhaniyogi et al. [2015]
    lambda <-
      sapply(1:G, function(each_region) {
        matrix(1, ranks, 2)
      }, simplify = F) # List of length G, each element a list of length k = 2, each element a matrix with dim = ranks x 2
    W <-
      lapply(betas, function(each_region) {
        sapply(each_region, function(each_dim) {
          array(1, dim = dim(each_dim))
        }, simplify = FALSE)
      }) # List of length 5, each element a list of length 2, each element a matrix dim p_j, ranks
    phi <-
      sapply(1:G, function(each_region) {
        rep(1 / ranks, ranks)
      }, simplify = F) # List of length G, each element a matrix of ranks rows by 2 columns
    sig2y <- 1 # Observational variance estimate
    if(ranks > 1){
      Xi <- sapply(seq(G),function(x) rep(.6,ranks - 1))
      if(!is.matrix(Xi)) Xi <- matrix(Xi,nrow = 1)
      cov_Metro <- sapply(seq(G),function(x) 0.01 * diag(ranks - 1),simplify = FALSE)
    }else{
      Xi <- matrix(1,1,G)
    }
  }

accept <- rep(0,G)

  # > Run MCMC ----
beginning_of_sampler <- proc.time()[3]
  for (s in 1:n.iter) {
    # >> Griddy-Gibbs ----
    params <- foreach(each_region = seq(G)) %dopar% {
      # Set up to sample alpha IF RANKS > 1
      if(ranks > 1){
        M = 4 # Number of reference sets per grid value of alpha
        l.weights <- sapply(alpha.g, function(proposed) {
          bw <- mapply(function(b, w) {
            abind(b, w, along = 3)
          },
          b = betas[[each_region]],
          w = W[[each_region]],
          SIMPLIFY = F)
          chi <- sapply(bw, function(each_dim) {
            apply(each_dim, 2, function(each_position) {
              t(each_position[, 1]) %*%
                diag(1 / each_position[, 2]) %*%
                each_position[, 1]
            })
          })
          ### INCLUDE RANK 1 change:
          if(ranks == 1){
            chi <- sum(chi)
          }else{
            chi <- apply(chi, 1, sum)
          }
          ## Draw Phi proposals
          ##### Phi under a stick-breaking prior
          old_Xi_g <- Xi[,each_region]
          phi.l <- sapply(seq(M),function(m){
            new_Xi_g <- c(old_Xi_g + cov_Metro[[each_region]] %*%
                            rnorm(ranks - 1))
            while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
              new_Xi_g <- c(old_Xi_g + cov_Metro[[each_region]] %*%
                              rnorm(ranks - 1))
            }
            new_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
              stick_break_log_posterior(new_Xi_g,cr, betas[[each_region]],W[[each_region]],tau[each_region],proposed)
            }))
            old_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
              stick_break_log_posterior(old_Xi_g,cr, betas[[each_region]],W[[each_region]],tau[each_region],proposed)
            }))
            if(exp(new_post_dens - old_post_dens) > runif(1)) old_Xi_g <- new_Xi_g
            stick_values(old_Xi_g)
          })

          ## Draw tau proposals
          ### ANOTHER RANK 1 CHANGE
          if(ranks == 1){
            chi2 <- chi / phi.l
          }else{
            chi2 <- apply(phi.l, 2, function(each_proposal) {
              chi / each_proposal
            })
            chi2 <- colSums(chi2)
          }
          tau.l <- rgig(M, a.tau - ranks * sum(p[[each_region]])/2, chi2, 2 * b.tau)
          refs <- list(phi = phi.l, tau = tau.l)
          ## Evaluate the densities
          lik.mean.tensor <-
            composeParafac(betas[[each_region]]) %o% input$x + array(1, dim = p[[each_region]]) %o% rep(1,TT) %o% d[each_region,]
          l.lik <-
            sum(dnorm(c(input$Y[[each_region]]), c(lik.mean.tensor), sqrt(sig2y), log = T))
          if(ranks == 1){
            l.bdens <- apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
              # Log prior density for all betas
              sapply(each_proposal[-1], function(each_rank_phi) {
                sum(unlist(sapply(bw, function(each_dim) {
                  apply(each_dim, 2, function(each_rank_bw) {
                    dnorm(
                      each_rank_bw[, 1],
                      0,
                      each_proposal[1] * each_rank_phi * each_rank_bw[, 2],
                      log = T
                    )
                  })
                })))
              })
            })
          }else{
            l.bdens <-
              colSums(apply(rbind(refs$tau, refs$phi), 2, function(each_proposal) {
                # Log prior density for all betas
                sapply(each_proposal[-1], function(each_rank_phi) {
                  sum(unlist(sapply(bw, function(each_dim) {
                    apply(each_dim, 2, function(each_rank_bw) {
                      dnorm(
                        each_rank_bw[, 1],
                        0,
                        each_proposal[1] * each_rank_phi * each_rank_bw[, 2],
                        log = T
                      )
                    })
                  })))
                })
              }))
          }

          l.tau <-
            dgamma(refs$tau, a.tau, b.tau, log = T) # Log prior density for tau

          ### RANK 1 CHANGE:
          if(ranks == 1){
            l.phi <- sapply(refs$phi,function(each_proposal){
              lgamma(ranks * proposed) - ranks * lgamma(proposed) + sum((rep(proposed, ranks) - 1) * log(each_proposal))
            })
          }else{
            l.phi <-
              apply(refs$phi, 2, function(each_proposal) {
                lgamma(ranks * proposed) - ranks * lgamma(proposed) + sum((rep(proposed, ranks) - 1) * log(each_proposal))
              })
          }
          # Log prior density for phi
          apply(cbind(l.phi, l.tau, l.bdens), 1, sum) + l.lik
        })
        mean.lweights <- apply(l.weights, 2, mean)
        weights <- exp(mean.lweights - max(mean.lweights))
        alpha <- sample(alpha.g, 1, prob = weights)
      }else{
        alpha <- 0
      }

      # >> Draw phi ----
      bg <- betas[[each_region]]
      wg <- W[[each_region]]
      ch <- mapply(function(b, w) {
        apply(abind(b, w, along = 3), 2, function(each_rank) {
          crossprod(each_rank[, 1], diag(1 / each_rank[, 2])) %*% each_rank[, 1]
        })
      }, b = bg, w = wg)

      ##### Phi under a stick-breaking prior
      old_Xi_g <- Xi[,each_region]
      if(ranks == 1){
        phi.g <- 1
      }else{
      accept <- results$accept[each_region]
      new_Xi_g <- c(old_Xi_g + cov_Metro[[each_region]]%*%rnorm(ranks - 1))
      while(length(new_Xi_g[new_Xi_g <= 0]) > 0){
        new_Xi_g <- c(old_Xi_g + cov_Metro[[each_region]]%*%rnorm(ranks - 1))
      }
      new_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
        stick_break_log_posterior(new_Xi_g,cr, betas[[each_region]],W[[each_region]],tau[each_region],alpha)
      }))
      old_post_dens <- sum(sapply(seq(ranks - 1),function(cr){
        stick_break_log_posterior(old_Xi_g,cr, betas[[each_region]],W[[each_region]],tau[each_region],alpha)
      }))
      if(exp(new_post_dens - old_post_dens) > runif(1)){
        old_Xi_g <- new_Xi_g
        accept <- accept + 1
      }
      phi.g <- stick_values(old_Xi_g)
      }

      # >> Draw tau ----
      xi <- a.tau - ranks * sum(p[[each_region]])/2
      ### RANK 1 ADJUSTMENT
      if(ranks == 1){
        chi <- sum(ch)
      }else{
        chi <- sum(apply(ch, 1, sum) / phi.g)
      }
      tau.g <- rgig(1, xi, chi, 2 * b.tau)

      # >> Draw lambda ----
      sumabsb <- sapply(bg, function(each_dim) {
        apply(each_dim, 2, function(each_rank) {
          sum(abs(each_rank))
        })
      })
      if(ranks == 1){
        lambda.g <- sapply(sumabsb,function(each_dim){
          rgamma(ranks,
                 a.lambda + p[[each_region]],
                 b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
        })
        lambda.g <- t(lambda.g)
      }else{
        lambda.g <-
          apply(sumabsb, 2, function(each_dim) {
            rgamma(ranks,
                   a.lambda + p[[each_region]],
                   b.lambda + (phi.g * tau.g) ^ (-.5) * each_dim)
          })
      }

      # >> Draw omega ----
      omega.g <- sapply(seq(D), function(each_dim) {
        sapply(seq(ranks), function(each_rank) {
          ch <-
            sapply(bg[[each_dim]][, each_rank], function(each_value) {
              each_value ^ 2 / (tau.g * phi.g[each_rank])
            })
          rgig(p[[each_region]][each_dim], 0.5, ch, lambda.g[each_rank, each_dim])
        }, simplify = T)
      }, simplify = F)

      # >> Draw beta ----
      # This next part has to be done sequentially, so the for loops are unavoidable
      for (each_dim in seq(D)) {
        for (each_rank in 1:ranks) {
          if(ranks == 1){
            expected <- array(1,dim=c(p[[each_region]],TT)) %o% d[each_region,]
          }else{
            expected <- composeParafac(lapply(bg, function(each_dim_I) {each_dim_I[, -each_rank,drop = FALSE]})) %o% input$x + array(1,dim=c(p[[each_region]],TT)) %o% d[each_region,]
          }
          if(G > 1){
            y.til <- (input$Y[[each_region]] - expected) %>%
            apply((D+1):(D+2),mode_k_matriz,each_dim) %>%
              array(dim=c(p[[each_region]][each_dim],prod(p[[each_region]][-each_dim]),TT,n))
          }else{
            y.til <- (input$Y - expected) %>%
            apply((D+1):(D+2),mode_k_matriz,each_dim) %>%
              array(dim=c(p[[each_region]][each_dim],prod(p[[each_region]][-each_dim]),TT,n))
          }
          betas_by_rank <-
            sapply(seq(ranks), function(rr) {
              sapply(seq(D), function(dd) {
                bg[[dd]][, rr]
              },simplify = FALSE)
            },simplify = FALSE) # This restructuring allows for the calculation
          vec_outer_other_betas <- c(Reduce(outer,betas_by_rank[[each_rank]][-each_dim]))
          var <-  1 /  (sum(input$x^2 %o% vec_outer_other_betas^2) / sig2y + (1 / (tau.g * phi.g[each_rank]) / diag(omega.g[[each_dim]][,each_rank]) ))
          mean_beta <- var %*% apply(y.til,1,function(dim_margin){
            sum(sapply(seq(TT),function(each_time){
              sapply(seq(n),function(each_subject){
                input$x[each_time,each_subject] * vec_outer_other_betas * dim_margin[,each_time,each_subject]
              })
            })) / sig2y
          })
          bg[[each_dim]][, each_rank] <-
            rnorm(p[[each_region]][each_dim], mean_beta, sqrt(diag(var)))
          if (each_dim > 1 && each_rank > 1) {
            bg[[each_dim]][1, each_rank] <-
              rtruncnorm(
                1,
                b = bg[[each_dim]][1, (each_rank - 1)],
                mean = (mean_beta)[1],
                sd = sqrt(var)[1]
              )
          }
        }
      }

      list(
        al = alpha,
        phi = phi.g,
        tau = tau.g,
        omega = omega.g,
        betas = bg,
        lambda = lambda.g,
        Xi = old_Xi_g[drop = FALSE],
        accept = accept
      )
    } # End dopar

    # >> Draw d ----
     y_hat <- mapply(function(y,parm){
       y - composeParafac(parm$betas) %o% input$x
     },y=input$Y,parm=params)

    inv_d_covar <- Sig + (TT*diag(sapply(p,prod)))^2 / sig2y
   #  #####################################################
   #  # Inverse Wishart prior on Sigma
   # (Or when the prior is applied to the covariance rather than the precision)
   # inv_d_covar <- chol2inv(chol(Sig)) + (TT*diag(sapply(p,prod)))^2 / sig2y
   #  #####################################################
   d_covar <- chol2inv(chol(inv_d_covar))
   d_mean <- d_covar%*% t(sapply(y_hat, function(yh){apply(yh,D+2,sum)}))

   d <- apply(d_mean,2,function(dm){
     dm + rnorm(length(dm)) %*% chol(d_covar)
   })

    # >> Draw Sig ----
    matd <- t(d)
    S <- crossprod(matd,matd)
    for(g in 1:G){
      delta <- rgamma(1,n/2 + 1,(S[g,g] + zeta)/2)
      S11inv <- chol2inv(chol(Sig[-g,-g]))
      varalpha <- chol2inv(chol((S[g,g] + zeta)*S11inv + diag(as.matrix(1/Upsilon[-g,g]))))
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

    ################################################
    # Inverse Wishart Prior on Sigma
   # inv_Sig <- MCMCpack::rwish(n+nu,VV + tcrossprod(d))
   # Sig <- chol2inv(chol((inv_Sig + t(inv_Sig))/2 + 1e-8))
    ################################################

    betas <- lapply(params, function(z) {
      z$betas
    })

    # >> Draw sig2y ----
    if(G > 1){
      sum_sq_diff <- mapply(function(y,parms,dd){
        (y - composeParafac(parms$betas) %o% input$x - array(1,dim=head(dim(y),-1))%o%dd)^2
      },y = input$Y,parms = params,dd = split(d,row(d))) %>% unlist %>% sum
    }else{
      sum_sq_diff <-
        (input$Y - composeParafac(params[[1]]$betas) %o% input$x)^2 %>%
        c %>%
        sum
    }

    sig2y <- 1/rgamma(1,
                      a.sig + n*TT*sum(sapply(p,prod))/2,
                      b.sig + (1/2)*sum(sum_sq_diff))

    W <- lapply(params, function(z) {
      z$omega
    })
    lambda <- lapply(params, function(z) {
      z$lambda
    })
    phi <- lapply(params, function(z) {
      z$phi
    })
    tau <- sapply(params, function(z) {
      z$tau
    })

    if(ranks > 1){
      Xi <- sapply(params, function(z) z$Xi[drop = FALSE])
      if(!is.matrix(Xi)) Xi <- t(as.matrix(Xi))
    }
    accept_all <- sapply(params,function(z) z$accept)

    # >> Get the log-likelihood ----
    if(save_llik == TRUE){
      llik <- -0.5*log(2*pi*sig2y)*n*TT*G*sum(sapply(p,prod)) - 0.5/sig2y * sum(sum_sq_diff)
      results$llik[s] <- llik
    }

    results$B[[s]] <- betas
    results$W[[s]] <- W
    results$lambda[[s]] <- lambda
    results$Phi[[s]] <- phi
    results$tau[, s] <- tau
    results$d[, , s] <- d
    results$Sig[, , s] <- Sig
    results$zeta[s] <- zeta
    results$sig2y[s] <- sig2y
    results$alpha[s] <- params[[1]]$al
    results$accept <- accept_all

    if(!is.null(save_after)){
      if(s %% save_after == 0){
        save(results,file = paste0("../Tests/405_",D,"D_rank_",ranks,"_first_",s,"_samples_",format(Sys.Date(),"%Y%m%d"),".RData"))
      }
    }

    if(s %% ceiling(n.iter/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time()," - Rank = ", ranks,
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
    results$Phi <- results$Phi[-(1:n.burn)]
    results$tau <- results$tau[, -(1:n.burn)]
    results$d <- results$d[, , -(1:n.burn)]
    results$Sig <- results$Sig[, , -(1:n.burn)]
    results$zeta <- results$zeta[-(1:n.burn)]
    results$sig2y <- results$sig2y[-(1:n.burn)]
  }

  results$accept <- results$accept / n.iter

  results
} # End MCMC function
