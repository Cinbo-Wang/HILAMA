# HILAMA

`HILAMA` is a R package to identify the nonzero acitive paths from high dimensional exposures, high dimensional mediators and continous outcome based on the Decorrelate \& Debiase method [1], which can provide a guarantee of  
false discovery rate control in finite sample settings.

## Installation

You can install the development version of `HILAMA` from Github via the `devtools` package.
```R
devtools::install_github('Cinbo-Wang/HILAMA')
```

## Examples

We first generate the simulation data.

```R
generate_param <- function(p,q,eta,kappa,s,r_p,r_pq,r_h){
  # p: dimension of exposure
  # q: dimension of mediator
  # s: dimension of hidden confounder
  # eta: effect size of hidden confounder
  # kappa: correlation size among exposure
  # r_h ：proportion of exposures and mediators affected by hidden confounder
  # r_p,r_qp: proportion of variables affecting the mediator, response.
  
  Sigma_E <- matrix(0,nrow=p,ncol=p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Sigma_E[i,j] <- (kappa)**abs(i-j)
    }
  }
  Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)
  
  Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
                nrow=s,byrow = T)
  
  if((p-ceiling(p*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:p,size=ceiling(p*r_h),replace = F)
      loc_zero <- setdiff(1:p,loc_active)
      Psi[i,loc_zero] <- 0
    }
  }
  
  Theta <- matrix(0,nrow=p,ncol=q)
  row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
  col_active_ls <- vector('list',length(row_active))
  for(i in 1:length(row_active)){
    col_active_ls[[i]] <- sample(1:q,size=sample(5:20,size=1))
  }
  col_active <- sort(unique(unlist(col_active_ls)))
  for(i in 1:length(row_active)){
    row_loc <- row_active[i]
    for(col_loc in col_active_ls[[i]]){
      Theta[row_loc,col_loc] <- runif(1,0.5,1.5)*sample(c(1,-1),1)  # Theta_hat[i,j]
    }
  }
  
  Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
                nrow=s,byrow=T)
  if((q-ceiling(q*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:q,size=ceiling(q*r_h),replace = F)
      loc_zero <- setdiff(1:q,loc_active)
      Phi[i,loc_zero] <- 0
    }
  }
  
  col_silence <- setdiff(1:q,col_active)
  s10 <- min(length(col_silence),ceiling(q*r_pq*0.2))
  if(s10>0){
    med_loc_false <- sample(col_silence,size=s10)
  }else{
    med_loc_false <- NULL
  }
  col_active_prob <- table(unlist(col_active_ls)) / length(unlist(col_active_ls))
  med_loc_true <- sample(as.numeric(names(col_active_prob)),size=q * r_pq - s10,prob = col_active_prob) # 被target越多，越有可能成为true mediator
  
  beta <- rep(0,q)
  for(i in c(med_loc_true,med_loc_false)){
    beta[i] <- runif(1,0.5,1.5)*sample(c(1,-1),1) 
  }
  
  exp_loc_true <- sample(1:p,ceiling(p*r_pq))
  gamma  <- rep(0,p);
  for(i in exp_loc_true){
    gamma[i] <- runif(1,0.5,1.5)*sample(c(1,-1),1)
  }
  
  phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
  
  return(list(
    Psi = Psi, Sigma_E=Sigma_E,Theta=Theta,Phi=Phi,
    gamma = gamma,beta=beta,phi = phi
  ))
  
}

generate_dataset <- function(n,param_ls){
  s <- nrow(param_ls[['Psi']]); p <- ncol(param_ls[['Psi']]); q <- ncol(param_ls[['Theta']])
  
  H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
  E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = param_ls[['Sigma_E']])
  X <- H %*% param_ls[['Psi']] + E
  
  E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
  M <- X %*% param_ls[['Theta']] + H %*% param_ls[['Phi']] + E_tilde
  
  e <- matrix(rnorm(n,0,1),nrow=n)
  Y <- X %*% param_ls[['gamma']] + M %*% param_ls[['beta']] + H %*% param_ls[['phi']]+e
  
  return(list(
    X = X, M = M, Y = Y,
    Theta=param_ls[['Theta']],beta=param_ls[['beta']],
    gamma =param_ls[['gamma']]
  ))
}

s <- 3; r_p <- r_pq <- 0.1; 
kappa <- 0.6; p <- 200; q <- 20
eta <- 1
r_h <- 0.5
n <- 300
param_ls <- generate_param(p,q,eta,kappa,s,r_p,r_pq,r_h)
data_ls <- generate_dataset(n,param_ls)
X <- data_ls$X
M_multi <- data_ls$M
Y <- data_ls$Y

Theta_true <- data_ls$Theta
beta_true <- data_ls$beta
gamma_true <- data_ls$gamma
nie_mat_true <- Theta_true %*% diag(beta_true)

```

start estimation using `HILAMA`

```R
require(HILAMA)
result_ls <- hilama(X, M_multi, Y,is.parallel=T,core_num=5,verbose=T,is.adap = F)

nie_mat_hat <- result_ls$nie_mat_hat
pvalue_Theta_mat <- result_ls$pvalue_Theta
pvalue_beta <- result_ls$pvalue_beta

### screen + JST test + BH 
screen_N <- ceiling(0.1*p*q)
result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat,pvalue_beta,screen_N=screen_N,
                                                  fdr_level=0.1)
pvalue_nie_screen_JST_mat <- result_screen_jst_ls$pvalue_screen_JST_mat
sig_level <- result_screen_jst_ls$pvalue_cutoff

result_eval_screen_JST <- cal_eval_metric_pvalue(nie_mat_true,nie_mat_hat,
                                                pvalue_nie_screen_JST_mat,sig_level = sig_level)
result_eval_screen_JST
```

## References

[1] Y. Sun, L. Ma, and Y. Xia. A decorrelating and debiasing approach to simultaneous
inference for high-dimensional confounded models, 2022.


