
#' Calculate the FDR and Power of the test given the true signal.
#'
#' @param Theta_true True value of signal matrix.
#' @param Theta_hat Estimated value of signal matrix.
#' @param pvalue_mat Estimated p-value of Theta_true.
#' @param sig_level The significance level.
#'
#' @return A vector including estimated FDR, Power and mean bias over nonzero locations.
#' @export
#'
#' @examples
cal_eval_metric <- function(Theta_true,Theta_hat,pvalue_mat,sig_level=0.05){

  loc_active <- which(abs(Theta_true) > 1e-5,arr.ind = T)
  loc_zero <- which(abs(Theta_true) <= 1e-5, arr.ind = T)

  bias <- Reduce('+',
                 apply(loc_active,1,function(x){
                   abs(Theta_hat[x[1],x[2]] - Theta_true[x[1],x[2]])
                 }))

  num_sig <- sum(pvalue_mat <= sig_level)
  num_tp <- Reduce('+',
                   apply(loc_active,1,function(x){
                     as.numeric(pvalue_mat[x[1],x[2]] <= sig_level)
                   }))
  num_fp <- num_sig - num_tp

  power <-  num_tp / nrow(loc_active)

  if(sum(pvalue_mat <= sig_level) ==0){
    fdr <- 1
  }else{
    fdr <- num_fp/num_sig
  }

  bias_mean <- bias/length(loc_active)
  eval_result <- c(fdr,power,bias_mean)
  names(eval_result) <- c('fdr','power','bias_mean')
  return(eval_result)
}


#' Get the p-value threshold for the mediation effect test
#'
#' @param pvalue_Theta_mat P value matrix of exposures to mediators.
#' @param pvalue_beta P value matrix of mediators to outcome.
#' @param screen_N Number of pairs used for multiple testing.
#' @param alpha The nominal False Discovery Rate level.
#'
#' @return
#' \item{pvalue_screen_JST_mat}{P value matrix of NIE using JST.}
#' \item{pvalue_cutoff}{P value threshold that controls the FDR at `\alpha` level.}
#'
#' @export get_nie_pvalue_screen_JST
#'
#' @examples
get_nie_pvalue_screen_JST <- function(pvalue_Theta_mat,pvalue_beta,screen_N=NULL,
                                      alpha=0.1){

  loc_screen0 <- which(apply(pvalue_Theta_mat,2,function(x)sum(x==1)<length(x)))

  p <- nrow(pvalue_Theta_mat); q <- ncol(pvalue_Theta_mat)
  pvalue_nie_mat_min <- matrix(1,p,q)
  for(i in 1:p){
    for(j in loc_screen0){
      pvalue_nie_mat_min[i,j] <- min(pvalue_Theta_mat[i,j],pvalue_beta[j])
    }
  }
  p_vec <- as.vector(pvalue_nie_mat_min)
  c <- max(p_vec[rank(p_vec,ties.method = 'min')<=screen_N])
  screen_loc <- which(pvalue_nie_mat_min <= c,arr.ind = T)


  pvalue_screen_JST_mat <- matrix(1,p,q)
  for(it in 1:nrow(screen_loc)){
    row_i <- screen_loc[it,1]; col_j <- screen_loc[it,2]
    pvalue_screen_JST_mat[row_i,col_j] <- max(pvalue_Theta_mat[row_i,col_j],pvalue_beta[col_j])
  }

  fdr_vec <- c();
  p_vec <- setdiff(sort(unique(as.vector(pvalue_screen_JST_mat)),decreasing = F),1)
  for(i in 1:min(1000,length(p_vec))){
    t <- p_vec[i]
    denom <- sum(pvalue_screen_JST_mat <= t)
    numer <- nrow(screen_loc)*t
    fdr_tmp <- numer / denom
    fdr_vec <- c(fdr_vec,fdr_tmp)
    if(t > 0.1) break
  }

  if(min(fdr_vec)>alpha){
    warning(paste0('minimal fdr is higher than alpha = ',alpha,'. The minimal fdr is ',min(fdr_vec)))
    sig_level <- p_vec[which.min(fdr_vec)]
  }else{
    sig_level <- max((p_vec[1:length(fdr_vec)])[fdr_vec <= alpha])
  }

  return(list(
    pvalue_screen_JST_mat = pvalue_screen_JST_mat,
    pvalue_cutoff = sig_level
  ))

}

