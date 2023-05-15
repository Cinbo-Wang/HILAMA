
#'  A function to identify the nonzero active paths from high dimensional exposures, high dimensional mediators and a continuous outcome
#'  while guaranteeing false discovery rate (FDR) control in finite sample settings.
#'
#' This function provides false discovery rate (FDR) control and support recovery in high-dimensional mediation models with
#' continuous outcome, where the exposure matrix, mediator matrix are potentially high-dimensional.
#'
#' @param X A matrix of input exposure data.
#' @param M_multi A matrix of input mediator data.
#' @param Y A vector of response data.
#' @param is.parallel A logical value indicating whether to use parallel computation, default is TRUE.
#' @param core_num The number of CPU cores to use in parallel computation, default is depends on your system.
#' @param verbose A logical value indicating whether to display progress info during computation, default is TRUE.
#' @param is.adap A logical value. If FALSE (default), the initial estimate in \code{infer_ddlasso()} uses LASSO. Otherwise, adaptive LASSO is used.
#'
#' @return
#' A list containing \code{nie_mat_hat}, \code{pvalue_Theta}, \code{pvalue_beta}, \code{nde_vec_hat},
#' and \code{pvalue_nde}, as well as \code{coef_ls}, which contains the intermediate calculation information.
#' @export
#'
#' @examples
hilama <- function(X, M_multi, Y, is.parallel=T, core_num=4, verbose=T, is.adap=F){

  n <- nrow(X); p <- ncol(X); q <- ncol(M_multi)

  ### Step1: X->M: column-wise ddlasso
  if(verbose) cat('Start step 1:  X -> M \n')
  score_all_ls <- esti_score_all(X,rho=0.5,is.parallel = is.parallel,core_num = core_num)
  result_Theta_ls <- list()
  for(j in 1:q){
    Mj <- M_multi[,j]
    coef_hat_df <- infer_ddlasso(X,Mj,rho=0.5,is.parallel = is.parallel,core_num = core_num,
                                 score_ls = score_all_ls,is.adap = is.adap)
    result_Theta_ls[[j]] <- coef_hat_df
    # if(j %% 10 ==0 ) cat('j:',j,'.\n')
  }

  pvalue_Theta<- Theta_hat <- matrix(0,nrow=p,ncol=q,dimnames = list(colnames(X),colnames(M_multi)))
  for(j in 1:q){
    tmp <- result_Theta_ls[[j]]
    u0 <- abs(tmp[,'dd_coef'] / tmp[,'std'])

    pvalue <- 2 * stats::pnorm(u0,lower.tail = F)
    pvalue_Theta[,j] <- pvalue
    Theta_hat[,j] <- tmp[,'dd_coef']
  }

  #### Step2: (X,M)->Y: XY, ddlasso
  # 针对每个Xj, 计算nodewise lasso

  if(verbose) cat('Start step 2: (X,M) -> Y \n')
  X_M <- cbind(X,M_multi)
  # plot(svd(X_M)$d)
  beta_gamma_hat_df <- infer_ddlasso(X=X_M,Y,rho=0.5,is.parallel = is.parallel,
                                     core_num = core_num, is.adap = is.adap)
  rownames(beta_gamma_hat_df) <- colnames(X_M)

  ### 存储中间结果
  coef_ls <- list(X2M_ls = result_Theta_ls,
                  beta_gamma_hat_df = beta_gamma_hat_df)

  #### Step3: 整理结果，计算NIE and pvalue adjust by BH and JST
  # coef_ls <- readRDS(tmp_file)
  if(verbose) cat('Start step 3: integrate results \n')
  gamma_hat <- beta_gamma_hat_df$dd_coef[1:p]
  beta_hat <- beta_gamma_hat_df$dd_coef[(p+1):(p+q)]
  pvalue_gamma <- beta_gamma_hat_df$pvalue[1:p]
  pvalue_beta <- beta_gamma_hat_df$pvalue[(p+1):(p+q)]

  nie_mat_hat <- Theta_hat %*% diag(beta_hat)
  rownames(nie_mat_hat)  <- colnames(X)
  colnames(nie_mat_hat)  <- colnames(M_multi)

  nde_vec_hat <- gamma_hat
  pvalue_nde <- pvalue_gamma
  names(nde_vec_hat) <- names(pvalue_nde) <- colnames(X)

  result_ls <- list(
    nie_mat_hat = nie_mat_hat, pvalue_Theta = pvalue_Theta,
    pvalue_beta = pvalue_beta,
    nde_vec_hat = nde_vec_hat, pvalue_nde = pvalue_nde,
    coef_ls = coef_ls
  )

  return(result_ls)
}




