

#' Estimate the score matrix for a given data matrix.
#'
#' Estimate the score matrix of a data matrix using node-wise Lasso on the transformed data matrix.
#' @param X A matrix data.
#' @param nlambda Number of hyper-parameter to search.
#' @param ratio_singular_min Used to check the existence of hidden confounder. If \eqn{\frac{\lambda_k}{\lambda_{k+1}}} is always smaller than this value,
#' then no hidden confounder exists and only single debiased procedure is applied.
#' @param Khat The number of possible hidden confounders.
#' @param parallel A logical value indicating whether to use parallel computation, default is TRUE.
#' @param core_num The number of CPU cores to use in parallel computation, default is depends on your system.
#'
#' @return
#' @export
#'
#' @examples
#'

esti_score_all <- function(X,nlambda=100,ratio_singular_min=1.5,Khat=NULL,parallel=F,core_num = 4){
  require(glmnet);require(doParallel);require(foreach);require(dplyr);require(stats)

  X <- scale(X,center=TRUE,scale=FALSE)

  n <- nrow(X); p <- ncol(X)
  svd.X <- svd(X)
  U <- svd.X$u; D <- svd.X$d; V <- svd.X$v

  D_lead <- lead(D)
  K_bar <- floor(min(n,p)/2)
  ratio_singular <- D[1:K_bar] / D_lead[1:K_bar]
  if(is.null(Khat)){
    # plot(ratio_singular);plot(D)
    if(max(ratio_singular)>ratio_singular_min){
      K_hat <- which.max(ratio_singular)
      D_tilde <- rep(1,ncol(U)); D_tilde[1:K_hat] <- 0
      F_dc <- diag(n) - U %*% diag(1 - D_tilde) %*% t(U)

    }else{
      K_hat <- 0
      F_dc <- diag(n)
    }
  }else{
    K_hat <- Khat
    D_tilde <- rep(1,ncol(U)); D_tilde[1:K_hat] <- 0
    F_dc <- diag(n) - U %*% diag(1 - D_tilde) %*% t(U)
  }

  X_tilde <- F_dc %*% X

  # nodewise lasso procedure
  Gamma = NULL
  denominator = NULL
  sequence = (nlambda:1) / nlambda * 3 * sqrt(log(p) / n)

  if(parallel){
    # library(doParallel);library(foreach)
    type <- ifelse(.Platform$OS.type=='windows','PSOCK','FORK')
    core_num <- ifelse(is.null(core_num),detectCores(logical = F)-1,core_num)

    cl <- makeCluster(core_num,type)
    registerDoParallel(cl)
    result_ls <- foreach(j = 1:p,.packages = 'glmnet')%dopar%{
      lambda_j = sqrt(sum(X_tilde[, j]^2) / n) * sequence
      fit_j = glmnet(x = X_tilde[, -j], y = X_tilde[, j], family = "gaussian",lambda = lambda_j, intercept = FALSE)
      temp = as.vector(crossprod(X_tilde[, j] - predict(fit_j, newx = X_tilde[, -j]), X_tilde[, j]))

      list(gamma = fit_j$beta,
           denom = temp)
    }
    stopImplicitCluster();
    stopCluster(cl)

    for(j in 1:p){
      tmp_ls <- result_ls[[j]]
      Gamma = rbind(Gamma, tmp_ls$gamma)
      denominator = rbind(denominator, tmp_ls$denom)
    }
  }else{
    for (j in 1:p){
      lambda_j = sqrt(sum(X_tilde[, j]^2) / n) * sequence
      fit_j = glmnet(x = X_tilde[, -j], y = X_tilde[, j], family = "gaussian",lambda = lambda_j, intercept = FALSE)

      Gamma = rbind(Gamma, fit_j$beta)
      temp = as.vector(crossprod(X_tilde[, j] - predict(fit_j, newx = X_tilde[, -j]), X_tilde[, j]))
      denominator = rbind(denominator, temp)
    }
  }


  index = 1
  IC = Inf
  for (i in 1:nlambda){
    Omega = matrix(0, ncol = p, nrow = p)
    gamma = matrix(Gamma[, i], nrow = (p - 1), ncol = p)

    Omega[upper.tri(Omega, diag = FALSE)] = - gamma[upper.tri(gamma, diag = FALSE)]  # 上三角
    Omega[lower.tri(Omega, diag = FALSE)] = - gamma[lower.tri(gamma, diag = TRUE)]
    diag(Omega) = 1
    rm(gamma)

    Omega = t(t(Omega) / (denominator[, i] / n))

    # symmetrization
    Omega = Omega * (abs(Omega) <= abs(t(Omega))) + t(Omega) * (abs(Omega) > abs(t(Omega)))
    #Omega = Matrix::Matrix(Omega, sparse = TRUE)
    # IC
    eig_min = RSpectra::eigs_sym(as(Omega, "matrix"), k = 1, which = "SA", opts = list(retvec = FALSE))$values
    # 	The k smallest (algebraic) eigenvalues, considering any negative sign.
    if (length(eig_min) == 0) {
      #eig_min = -Inf # algorithm may not be converged # 原本这里是—Inf
      eig.ls <- eigen(Omega)
      eig_min <- min(eig.ls$values)
    }

    if (eig_min <= 0){
      next
    }else {
      IC_new = - determinant(Omega)$modulus[1] + log(n) / n * (sum(Omega != 0) - p) / 2 +
        sum(Omega * crossprod(X_tilde) / n)
      if (IC_new <= IC){
        IC = IC_new
        index = i
      }
    }
  }

  gamma = Matrix::Matrix(Gamma[, index], nrow = (p - 1), ncol = p, sparse = TRUE)
  denominator = denominator[, index]
  rm(Gamma, Omega)

  score_ls <- list()
  for(j in 1:p){
    z_j <- as.vector(X_tilde[, j] - X_tilde[, -j] %*% gamma[, j])
    score_ls[[j]] <- z_j
  }
  return(list(score_ls=score_ls,
              Khat = K_hat))

}



#' Infer the coefficient of high dimensional linear model with latent confounder using Decorrelating and Debiasing approach.
#'
#' This function estimates and infers the coefficients of a high-dimensional linear model with a latent confounder using the ecorrelating and Debiasing approach proposed by Sun, Ma, and Xia (2023).
#'
#' @param X A data matrix.
#' @param Y Continuous response vector.
#' @param ratio_singular_min Used to check the existence of hidden confounder. If \eqn{\frac{\lambda_k}{\lambda_{k+1}}} is always smaller than this value,
#' then no hidden confounder exists and only single debiased procedure is applied.
#' @param Khat The number of possible hidden confounders.
#' @param score_ls Object returned from esti_score_all() function.
#'
#' @references Sun, Y., Ma, L., & Xia, Y. (2023). A decorrelating and debiasing approach to simultaneous inference for high-dimensional confounded models. Journal of the American Statistical Association, 1-12.
#' @return
#' @export
#'
#' @examples
infer_ddlasso <- function(X,Y,ratio_singular_min=1.5,Khat=NULL,
                          score_ls = NULL){
  require(glmnet);require(stats)

  if(is.null(score_ls)) {
    stop('Please supply the score_ls object returned by esti_score_all() function.')
  }
  n <- nrow(X); p <- ncol(X)
  X <- scale(X,center=TRUE,scale=FALSE) # score_ls is obtained after centering X.
  Y <- Y - mean(Y)

  svd.X <- svd(X)
  U <- svd.X$u; D <- svd.X$d; V <- svd.X$v
  D_lead <- lead(D)
  K_bar <- min(floor(min(n,p)/2),20)

  ## decide whether transform the original data (whether hidden confounder exits.)
  ratio_singular <- D[1:K_bar] / D_lead[1:K_bar]
  if(is.null(Khat)){
    # plot(ratio_singular);plot(D)
    if(max(ratio_singular)>ratio_singular_min){
      K_hat <- which.max(ratio_singular)
      D_tilde <- rep(1,ncol(U)); D_tilde[1:K_hat] <- 0
      F_dc <- diag(n) - U %*% diag(1 - D_tilde) %*% t(U)
    }else{
      K_hat <- 0
      F_dc <- diag(n)
    }
  }else{
    K_hat <- Khat
    D_tilde <- rep(1,ncol(U)); D_tilde[1:K_hat] <- 0
    F_dc <- diag(n) - U %*% diag(1 - D_tilde) %*% t(U)
  }

  # Decorrelation
  Xtilde <- F_dc %*% X; Ytilde <- F_dc %*% Y
  cv.fit <- cv.glmnet(x=Xtilde, y=Ytilde,lambda.min.ratio=0.01,intercept=T)
  beta_hat_init <- as.matrix(coef(cv.fit,s=cv.fit$lambda.min)[-1])
  error_vec <- Ytilde - predict(cv.fit,Xtilde,s='lambda.min')
  residual_init <- Y - predict(cv.fit,X,s='lambda.min')

  signal_hat = which(beta_hat_init != 0)
  sigma_hat = norm(error_vec,type='2')**2
  if (length(signal_hat) == 0){
    sigma_hat = sigma_hat / sum(F_dc * F_dc)
    sigma_xi_hat = sqrt(sigma_hat)
  }else {
    UDV_list = svd(Xtilde[, signal_hat], nv = 0)
    rank = sum(UDV_list$d != 0)
    projection = tcrossprod(UDV_list$u[, 1:rank])
    projection = (diag(1, n) - projection) %*% F_dc
    sigma_hat = sigma_hat / sum(projection * projection)
    sigma_xi_hat = sqrt(sigma_hat)
  }

  ## Double debiased
  coef_mat <- matrix(0,nrow=p,ncol=3)
  colnames(coef_mat) <- c('init_coef','dd_coef','sd')
  for(id in 1:p){
    z_score <- score_ls[[id]]

    bias <- t(z_score) %*% F_dc  %*% residual_init / (t(z_score) %*% F_dc %*% X[,id])
    coef_debiased <- beta_hat_init[id] + bias
    sd_coef_id <- sigma_xi_hat / sqrt(sum(z_score**2))
    coef_mat[id,] <- c(beta_hat_init[id],coef_debiased,sd_coef_id)

  }

  z_stat <- coef_mat[,2] / coef_mat[,3]
  pvalue <- 2*(1- pnorm(abs(z_stat)))

  coef_mat <- cbind(coef_mat,z_stat,pvalue)
  coef_df <- as.data.frame(coef_mat)
  coef_df$init_coef <- coef_mat[,1]
  coef_df$dd_coef <- coef_mat[,2]
  coef_df$sd <- coef_mat[,3]
  return(coef_df)
}



#' A function to identify the nonzero active paths from high dimensional exposures, high dimensional mediators and a continuous outcome
#' while guaranteeing false discovery rate (FDR) control in finite sample settings.
#'
#' This function provides false discovery rate (FDR) control and support recovery in high-dimensional mediation models with
#' continuous outcome, where the exposure matrix, mediator matrix are potentially high-dimensional.
#'
#' @param X A matrix of input exposure data.
#' @param M_multi A matrix of input mediator data.
#' @param Y A vector of continuous response data.
#' @param Z A matrix of baseline covariates.
#' @param parallel A logical value indicating whether to use parallel computation, default is TRUE.
#' @param core_num The number of CPU cores to use in parallel computation, default is depends on your system.
#' @param alpha0 P value threshold  to screen unnecessary mediators for the following mediator model.
#' @param cut_off_ratio_singular_value Used to check the existence of hidden confounder. If \eqn{\frac{\lambda_k}{\lambda_{k+1}}} is always smaller than this value,
#' then no hidden confounder exists and only single debiased procedure is applied.
#' @param K_med Number of hidden confounder in the mediator model Y ~ M_multi.
#' @param K_or Number of hidden confounder in the outcome model Y ~ M_multi + X.
#' @param verbose A logical value indicating whether to display progress info during computation, default is TRUE.
#'
#' @return
#' \item{nie_mat_hat}{}
#' \item{Theta_hat}{}
#' \item{pvalue_Theta}{}
#' \item{beta_hat}{}
#' \item{pvalue_beta}{}
#' \item{nde_vec_hat}{}
#' \item{pvalue_nde}{}
#' \item{K_med}{}
#' \item{K_or}{}

#' @author Xinbo Wang
#'
#' @export hilama_fast
#'
#' @examples
#' \dontrun{
#' p <- 300
#' q <- 100
#' eta <- 1
#' rho <- 0.5
#' s <- 2
#' p_z <- 3
#' r_p <- r_pq <- 0.1
#' n <- 400

#' Sigma_E <- matrix(0,nrow=p,ncol=p)
#' for(i in 1:(p-1)){
#'   for(j in (i+1):p){
#'     Sigma_E[i,j] <- (rho)**abs(i-j)
#'   }
#' }
#' Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)

#' Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
#'               nrow=s,byrow = T)
#' Psi_z <- matrix(rnorm(n=p_z*p,mean=0,sd=1),
#'                 nrow=p_z,byrow = T)

#' Theta <- matrix(0,nrow=p,ncol=q)
#' row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
#' col_active_ls <- vector('list',length(row_active))
#' for(i in 1:length(row_active)){
#'   col_active_ls[[i]] <- sample(1:q,size=sample(5:50,size=1))
#' }
#' for(i in 1:length(row_active)){
#'   row_loc <- row_active[i]
#'   for(col_loc in col_active_ls[[i]]){
#'     Theta[row_loc,col_loc] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
#'   }
#' }

#' Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
#'               nrow=s,byrow=T)
#' Phi_z <- matrix(rnorm(p_z*q,mean=eta,sd=1)*sample(c(1,-1),size=p_z*q,replace = T),
#'                 nrow=p_z,byrow=T)

#' col_active <- sort(unique(unlist(col_active_ls)))
#' col_silence <- setdiff(1:q,col_active)
#' s10 <- min(length(col_silence),ceiling(q*r_pq*0.2))
#' if(s10>0){
#'   med_loc_false <- sample(col_silence,size=s10)
#' }else{
#'   med_loc_false <- NULL
#' }
#' col_active_prob <- table(unlist(col_active_ls)) / length(unlist(col_active_ls))
#' med_loc_true <- sample(as.numeric(names(col_active_prob)),size=q * r_pq - s10,prob = col_active_prob) # 被target越多，越有可能成为true mediator

#' beta <- rep(0,q)
#' for(i in c(med_loc_true,med_loc_false)){
#'   beta[i] <-  rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
#' }

#' exp_loc_true <- sample(1:p,ceiling(p*r_pq))
#' gamma  <- rep(0,p);
#' for(i in exp_loc_true){
#'   gamma[i] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
#' }

#' phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
#' phi_z <- rnorm(p_z,eta,1)*sample(c(1,-1),size=p_z,replace = T)



#' H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
#' Z <- matrix(rnorm(n*p_z,0,1),nrow=n,byrow=T)
#' E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigma_E)
#' X <- Z %*% Psi_z + H %*% Psi + E

#' E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
#' M_multi <- X %*% Theta + Z %*% Phi_z + H %*% Phi + E_tilde

#' e <- matrix(rnorm(n,0,1),nrow=n)
#' Y <- X %*% gamma + M_multi %*% beta + Z %*% phi_z + H %*% phi + e

#' result_ls <- hilama_fast(X, M_multi, Y,Z,parallel=T,core_num=4)


#' nie_mat_hat <- result_ls$nie_mat_hat
#' pvalue_Theta_mat <- result_ls$pvalue_Theta
#' pvalue_beta <- result_ls$pvalue_beta
#' nie_mat_true <- Theta %*% diag(beta)

#' loc_screen <- which(apply(pvalue_Theta_mat,2,function(x)sum(x==1)<length(x)))
#' screen_N <- ceiling(0.1*p*length(loc_screen))
#' result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat,pvalue_beta,screen_N=screen_N,
#'                                                   alpha=0.1)
#' pvalue_nie_screen_JST_mat <- result_screen_jst_ls$pvalue_screen_JST_mat
#' result_eval_screen_JST <- cal_eval_metric(nie_mat_true,pvalue_nie_screen_JST_mat,sig_level=result_screen_jst_ls$pvalue_cutoff)
#' result_eval_screen_JST
#'
#' }



#' @import glmnet
#' @import foreach
#' @import doParallel
#' @import stats
#'
hilama_fast <- function(X, M_multi, Y, Z = NULL,
                        parallel=T, core_num=4, alpha0 = 0.1,
                        cut_off_ratio_singular_value = 1.5, K_med=NULL,K_or = NULL,verbose=T){

  n <- nrow(X); p <- ncol(X); q <- ncol(M_multi)

  if(verbose) message('Start step 1: (X,M) -> Y \n')
  X_M <- cbind(X,M_multi)

  if(!is.null(Z)){
    P_Z <- Z %*% solve(t(Z)%*%Z) %*% t(Z)
    X_M_tilde <- (diag(n)- P_Z) %*% X_M
    Y_tilde <- (diag(n) - P_Z) %*% Y
  }else{
    X_M_tilde <- X_M
    Y_tilde <- Y
  }

  score_X_M_ls <- esti_score_all(X = X_M_tilde,ratio_singular_min = cut_off_ratio_singular_value,
                                 Khat = K_or,parallel = parallel,core_num = core_num)
  K_or <- score_X_M_ls$Khat
  beta_gamma_hat_df <- infer_ddlasso(X=X_M_tilde,Y=Y_tilde,score_ls = score_X_M_ls$score_ls)

  rownames(beta_gamma_hat_df) <- colnames(X_M)
  pvalue_M <- beta_gamma_hat_df$pvalue[-(1:p)]
  loc_screen_M <- which(pvalue_M <= alpha0)

  if(length(loc_screen_M)>0){
    if(verbose) message('Start step 2:  X -> M \n')
    if(!is.null(Z)){
      P_Z <- Z %*% solve(t(Z)%*%Z) %*% t(Z)
      M_multi_tilde <- (diag(n)- P_Z) %*% M_multi
      X_tilde <- (diag(n) - P_Z) %*% X
    }else{
      M_multi_tilde <- M_multi
      X_tilde <- X
    }

    score_X_ls <- esti_score_all(X=X_tilde,ratio_singular_min = cut_off_ratio_singular_value,
                                 Khat = K_med,parallel = parallel,core_num = core_num)
    K_med <- score_X_ls$Khat
    result_Theta_ls <- vector('list',length=ncol(M_multi))
    for(j in loc_screen_M){
      Mj <- M_multi_tilde[,j]
      coef_hat_df <- infer_ddlasso(X=X_tilde,Y=Mj, score_ls = score_X_ls$score_ls)
      result_Theta_ls[[j]] <- coef_hat_df
      # if(j %% 10 ==0 ) cat('j:',j,'.\n')
    }
    pvalue_Theta <- matrix(1,nrow=p,ncol=q,dimnames = list(colnames(X),colnames(M_multi)))
    Theta_hat <- matrix(0,nrow=p,ncol=q,dimnames = list(colnames(X),colnames(M_multi)))
    Z_Theta_hat <- matrix(0,nrow=p,ncol=q,dimnames = list(colnames(X),colnames(M_multi)))
    for(j in loc_screen_M){
      tmp <- result_Theta_ls[[j]]
      u0 <- abs(tmp[,'dd_coef'] / tmp[,'sd'])

      pvalue <- 2* pnorm(u0,lower.tail = F)
      pvalue_Theta[,j] <- pvalue
      Theta_hat[,j] <- tmp[,'dd_coef']
      Z_Theta_hat[,j] <- u0
    }

    coef_ls <- list(X2M_ls = result_Theta_ls,
                    beta_gamma_hat_df = beta_gamma_hat_df)

    if(verbose) message('Start step 3: integrate results \n')
    gamma_hat <- beta_gamma_hat_df$dd_coef[1:p]
    beta_hat <- beta_gamma_hat_df$dd_coef[(p+1):(p+q)]
    pvalue_gamma <- beta_gamma_hat_df$pvalue[1:p]
    pvalue_beta <- beta_gamma_hat_df$pvalue[(p+1):(p+q)]

    Z_gamma_beta_hat <- abs(beta_gamma_hat_df$dd_coef / beta_gamma_hat_df$sd)

    nie_mat_hat <- Theta_hat %*% diag(beta_hat)
    rownames(nie_mat_hat)  <- colnames(X)
    colnames(nie_mat_hat)  <- colnames(M_multi)

    nde_vec_hat <- gamma_hat
    pvalue_nde <- pvalue_gamma
    names(nde_vec_hat) <- names(pvalue_nde) <- colnames(X)

    result_ls <- list(
      nie_mat_hat = nie_mat_hat,
      Theta_hat = Theta_hat, pvalue_Theta = pvalue_Theta,
      beta_hat = beta_hat, pvalue_beta = pvalue_beta,
      nde_vec_hat = nde_vec_hat, pvalue_nde = pvalue_nde,
      K_med = K_med, K_or = K_or
    )
    return(result_ls)
  }else{
    warning('No significant mediators!')
    return(NULL)
  }
}



