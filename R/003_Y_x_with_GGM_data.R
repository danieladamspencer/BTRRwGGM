#' Simulate data for the Bayesian Tensor Response Regression with Gaussian Graphical Modeling
#'
#' @param subjects A scalar. The number of subjects in the sample.
#' @param regions A scalar. The number of distinct tensor responses for each subject at each time step.
#' @param max_time A scalar. The number of time steps to be simulated from.
#' @param avg_margin_size A vector of length *D* that gives the average margin lengths for the tensor reponse *regions* for a given time and subject.
#' @param SNR A scalar. The signal-to-noise ratio as defined by Welvaert and Rosseel, 2013
#' @param CNR A scalar. The contrast-to-noise ratio as defined by Welvaert and Rosseel, 2013
#' @param conn_regions A scalar. The number of pairs of regions that have significant connectivity.
#' @param conn_level A scalar. The strength of the correlation between the connected regions, between 0 and 1
#'
#' @return A list with five elements, the response *Y*, the covariate *x*, and the true values for $d_{i,g}$, $\mathbf{B}_g$ and $\boldsymbol{\Sigma}$.
#' @export
#'
#' @examples
#' simulated_data <- Y_x_with_GGM_data()
Y_x_with_GGM_data <- function(subjects = 50, regions = 5, max_time = 100,
                              avg_margin_size = c(10,10,10),SNR = 1, CNR = 1,
                              conn_regions = 1, conn_level = 0.95){
  require(neuRosim)
  d_covar <- diag(regions)
  conns <- matrix(NA,conn_regions,2)
  for(g in seq(conn_regions)){
    conns[g,] <- sample(seq(regions),size=2,replace = FALSE)
    if(g >= 2){
      for(gg in seq(g - 1)){
        while(identical(conns[g,], conns[gg,]) |
              identical(conns[g,], rev(conns[gg,]))) {
          conns[g,] <- sample(seq(regions),size=2,replace = FALSE)}
      }
    }
    d_covar[conns[g,1],conns[g,2]] <-
      d_covar[conns[g,2],conns[g,1]] <- conn_level
  }
  d <- matrix(rnorm(subjects*regions),subjects,regions) %*% d_covar * SNR

  p <- sapply(seq(regions),function(r) sapply(avg_margin_size,rpois,n=1),
              simplify = FALSE)
  B <- sapply(p,function(region_idx){
    neuRosim::specifyregion(dim=region_idx,
                            coord = runif(length(region_idx))*region_idx,
                            radius = min(round(min(region_idx*0.1)),
                                         sample(3:5,1)),
                            form = "sphere")*CNR
  },simplify = FALSE)
  x <- canonicalHRF(seq(max_time),param=list(
    a1 = round(max_time*0.12), # Delay of response relative to onset
    a2 = round(max_time*0.5), # Delay of undershoot relative to onset
    b1 = 2, # Dispersion of response
    b2 = 1, # Dispersion of undershoot
    c = 0.5 # Scale of undershoot
  ))

  Y <- mapply(function(b,dd){
    outer(outer(b,x,FUN = `*`),dd,FUN = `+`) +
      array(rnorm(prod(dim(b))*max_time*subjects),
            dim = c(dim(b),max_time,subjects))
  },b=B,dd=split(d,col(d)),SIMPLIFY = FALSE)

  output <-
    list(
      Y = Y,
      x = x,
      true_d = d,
      true_B = B,
      true_d_covar = d_covar
    )
  return(output)
}


#' Simulate fMRI Data for the Bayesian Tensor Response Regression with Gaussian Graphical Modeling
#'
#' @param subjects A scalar. The number of subjects in the sample.
#' @param regions A scalar. The number of distinct tensor responses for each subject at each time step.
#' @param max_time A scalar. The number of time steps to be simulated from.
#' @param avg_margin_size A vector of length *D* that gives the average margin lengths for the tensor reponse *regions* for a given time and subject.
#' @param SNR A scalar. The signal-to-noise ratio as defined by Welvaert and Rosseel, 2013
#' @param CNR A scalar. The contrast-to-noise ratio as defined by Welvaert and Rosseel, 2013
#' @param conn_regions A scalar. The number of pairs of regions that have significant connectivity.
#' @param conn_level A scalar. The strength of the correlation between the connected regions, between 0 and 1
#' @param block_period A scalar. The length of the block period for the fMRI block trials
#'
#' @return A list with five elements, the response *Y*, the covariate *x*, and the true values for $d_{i,g}$, $\mathbf{B}_g$ and $\boldsymbol{\Sigma}$.
#' @export
#'
#' @examples
#' simulated_fmri_data <- Y_x_with_GGM_fmri_data()
Y_x_with_GGM_fmri_data <- function(subjects = 50, regions = 5, max_time = 100,
                              avg_margin_size = c(10,10,10),SNR = 1, CNR = 1,
                              conn_regions = 1, conn_level = 0.95,
                              block_period = max_time / 10){
  require(neuRosim)
  require(fmri)
  d_covar <- diag(regions)
  conns <- matrix(NA,conn_regions,2)
  for(g in seq(conn_regions)){
    conns[g,] <- sample(seq(regions),size=2,replace = FALSE)
    if(g >= 2){
      for(gg in seq(g - 1)){
        while(identical(conns[g,], conns[gg,]) |
              identical(conns[g,], rev(conns[gg,]))) {
          conns[g,] <- sample(seq(regions),size=2,replace = FALSE)}
      }
    }
    d_covar[conns[g,1],conns[g,2]] <-
      d_covar[conns[g,2],conns[g,1]] <- conn_level
  }
  # if(!is.null(scenario)){
  #   regions <- 10
  #   if(scenario == 1){
  #     d_covar <- matrix(
  #       c(
  #         1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  #         1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  #         0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  #         0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,
  #         0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,
  #         0.0,0.0,-1.0,-1.0,-1.0,1.0,0.0,0.0,0.0,0.0,
  #         0.0,0.0,-1.0,-1.0,-1.0,1.0,1.0,0.0,0.0,0.0,
  #         0.0,0.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,0.0,0.0,
  #         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,
  #         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0
  #       ), nrow = 10, ncol = 10, byrow = TRUE
  #     )
  #   }
  # }

  d <- matrix(rnorm(subjects*regions),subjects,regions) %*% d_covar * SNR

  p <- sapply(seq(regions),function(r) sapply(avg_margin_size,rpois,n=1),
              simplify = FALSE)
  B <- sapply(p,function(region_idx){
    neuRosim::specifyregion(dim=region_idx,
                            coord = runif(length(region_idx))*region_idx,
                            radius = min(round(min(region_idx*0.1)),
                                         sample(3:5,1)),
                            form = "sphere")*CNR
  },simplify = FALSE)

  k <- round(max_time / block_period)
  obs_x <- vector("numeric",length = max_time)
  for(i in seq(0,(k-1))){
    obs_x[1+((i*block_period):(i*block_period + round(block_period/2)))] <- 1
  }
  hrf <- canonicalHRF(seq(max_time))
  x <- convolve(obs_x,rev(hrf),type = "o")[seq(max_time)]

  Y <- mapply(function(b,dd){
    outer(outer(b,x,FUN = `*`),dd,FUN = `+`) +
      array(rnorm(prod(dim(b))*max_time*subjects),
            dim = c(dim(b),max_time,subjects))
  },b=B,dd=split(d,col(d)),SIMPLIFY = FALSE)

  output <-
    list(
      Y = Y,
      x = x,
      true_d = d,
      true_B = B,
      true_d_covar = d_covar
    )
  return(output)
}
