#' Compose a tensor from its CANDECOMP/PARAFAC (CP) decomposition
#'
#' This function takes a list of length D containing all of
#' the components of the CP decomposition and returning a D-dimensiona
#' tensor.
#'
#' @param bb A list of length D in which each element is a p_d by R matrix
#'
#' @return A single array-class tensor. In two-dimensions, this will be returned as a matrix.
#' @export
#'
composeParafac <- function(bb){
  DD <- length(bb)
  pp <- lapply(bb,nrow)
  RR <- lapply(bb,ncol)
  stopifnot(all(unlist(RR) == unlist(RR)[1]))
  RRR <- unique(unlist(RR))
  cp_summands <- sapply(1:RRR,function(r){
    rank_betas <- lapply(bb,function(b){b[,r]})
    Reduce("%o%",rank_betas)
  },simplify = "array")
  apply(cp_summands,1:DD,sum)
}


tensor_image_2D <- function(mat_object){
  require(tidyverse)
  require(magrittr)
  n_rows <- dim(mat_object)[1]
  n_cols <- dim(mat_object)[2]
  mat_object %>%
    t %>%
    as.data.frame %>%
    magrittr::set_colnames(as.character(1:n_rows)) %>%
    cbind(row_num=1:n_cols) %>%
    tidyr::gather(key=col_num,"value",-row_num) %>%
    dplyr::mutate(col_num = as.numeric(col_num)) %>%
    ggplot() + geom_tile(aes(x=col_num,y=row_num,fill=value,color=value)) +
    scale_color_distiller(palette="Greys") +
    scale_fill_distiller(palette="Greys") +
    theme_bw() + # theme(legend.position = "none") +
    labs(y="",x="")
}

mode_k_matriz <- function(tens,k){
  if(k > length(dim(tens))) stop("You cannot have a mode greater than the dimension of the tensor!")
  t(apply(tens,k,c))
}

# Function to calculate Phi_g given Xi_g
stick_values <- function(breaks){
  if(!is.numeric(breaks)) stop("You need to input a numeric variable.")
  # if(length(breaks) < 2) stop("Your input should be of length at least 2.")
  out <- numeric(length(breaks) + 1)
  for(i in seq(length(breaks))){
    out[i] <- breaks[i]*prod(1 - breaks[seq(length(breaks)) < i])
  }
  out[length(breaks) + 1] <- 1 - sum(out[seq(length(breaks))])
  return(out)
}

# Log posterior function for the stick-breaking construction of Phi_g
stick_break_log_posterior <- function(Xi_g, current_rank, betas_g, omega_g, tau_g,alpha_g){
  require(abind)
  model_rank <- ncol(betas_g[[1]])
  if(current_rank == model_rank) stop("This only needs to be done for r = 1,...,R-1!")
  dimension_lengths <- sapply(betas_g,function(x) dim(x)[1])
  part_a <- log(Xi_g[current_rank])*(-sum(dimension_lengths)/2)
  part_b <- (alpha_g - (model_rank - current_rank)*sum(dimension_lengths)/2 - 1)  * log(1 - Xi_g[current_rank])
  BWB <- mapply(function(b,w){
    bw_bind <- abind(b,w,along = 3)
    apply(bw_bind,2,function(bwb){
      t(bwb[,1])%*%diag(1/bwb[,2])%*%bwb[,1]
    })
  },b = betas_g,w = omega_g)
  dim_sum_BWB <- apply(BWB,1,sum)
  part_c <- (1/Xi_g[current_rank]) * dim_sum_BWB[current_rank]
  greater_ranks <- if(current_rank < model_rank - 1) seq(current_rank + 1, model_rank - 1)
  if(is.null(greater_ranks)){part_d = 0}else{
    part_d <- sum(unlist(sapply(greater_ranks,function(each_rank){
      (1 / (Xi_g[each_rank] * prod(1 - Xi_g[seq(each_rank - 1)]))) * dim_sum_BWB[each_rank]
    })))
  }
  part_e <- (1 / stick_values(Xi_g)[model_rank])*dim_sum_BWB[model_rank]
  out <- part_a + part_b - (1/tau_g) * (part_c + part_d + part_e)
  return(out)
}

# A c++ implementation of lapply
# require(Rcpp)
# require(inline)
# src <- '
#    Rcpp::List input(data);
# Rcpp::Function f(fun);
# Rcpp::List output(input.size());
# std::transform(input.begin(), input.end(), output.begin(), f);
# output.names() = input.names();
# return output;
# '
# cpp_lapply <- cxxfunction(signature(data = "list", fun = "function"),
#                           src, plugin = "Rcpp")


make_region_array <- function(y){
  lapply(sapply(1:G, function(each_region) {
    lapply(y, function(each_subject) {
      each_subject[[each_region]]
    })
  }, simplify = F), function(each_subject_and_region) {
    array(unlist(each_subject_and_region), dim = c(
      dim(each_subject_and_region[[1]]),
      length(each_subject_and_region)
    ))
  })
}

single_region_array <- function(y, region) {
  sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array")
}

require(tidyverse)
std_single_region_array <- function(y,region,Dim,TT,p,n){
  a <- sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array")
  b <- apply(a,Dim+2,function(each_subject){
      y_mean <- apply(each_subject,seq(Dim),mean)
      y_sd <- apply(each_subject,seq(Dim),sd)
      centered_y <- apply(each_subject,Dim+1,function(each_time){
        (each_time - y_mean) / y_sd
      })
      centered_y[is.nan(centered_y)] <- 0
      array(centered_y,dim = dim(each_subject)[seq(Dim+1)])
    })
  cc <- array(b,dim = c(p[[region]],TT,n))
  return(cc)
}

vec_std_single_region_array <- function(y,region,p){
  sapply(y, function(each_subject) {
    each_subject[[region]]
  },simplify = "array") %>%
    apply(D+2,function(each_subject){
      y_mean <- apply(each_subject,seq(D),mean)
      y_sd <- apply(each_subject,seq(D),sd)
      centered_y <- apply(each_subject,D+1,function(each_time){
        (each_time - y_mean) / y_sd
      })
      centered_y[is.nan(centered_y)] <- 0
      array(centered_y,dim = c(p[[region]],TT))
    }) %>%
    array(dim = c(prod(p[[region]]),TT,n))
}
