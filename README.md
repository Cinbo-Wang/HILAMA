# HILAMA

`HILAMA` is a R package to identify the nonzero acitive paths from high dimensional exposures, high dimensional mediators and continous outcome based on the Decorrelating \& Debiasing Approach [1], which can provide a guarantee of false discovery rate control in finite sample settings.

## Installation

You can install the development version of `HILAMA` from Github via the `devtools` package.
```R
devtools::install_github('Cinbo-Wang/HILAMA')
```

## Examples

We first generate the simulation data.

```R
p <- 300
q <- 100
eta <- 1
rho <- 0.5
s <- 2
p_z <- 3
r_p <- r_pq <- 0.1
n <- 400

Sigma_E <- matrix(0,nrow=p,ncol=p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    Sigma_E[i,j] <- (rho)**abs(i-j)
  }
}
Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)

Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
              nrow=s,byrow = T)
Psi_z <- matrix(rnorm(n=p_z*p,mean=0,sd=1),
                nrow=p_z,byrow = T)

Theta <- matrix(0,nrow=p,ncol=q)
row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
col_active_ls <- vector('list',length(row_active))
for(i in 1:length(row_active)){
  col_active_ls[[i]] <- sample(1:q,size=sample(5:50,size=1))
}
for(i in 1:length(row_active)){
  row_loc <- row_active[i]
  for(col_loc in col_active_ls[[i]]){
    Theta[row_loc,col_loc] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
  }
}

Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
              nrow=s,byrow=T)
Phi_z <- matrix(rnorm(p_z*q,mean=eta,sd=1)*sample(c(1,-1),size=p_z*q,replace = T),
                nrow=p_z,byrow=T)

col_active <- sort(unique(unlist(col_active_ls)))
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
  beta[i] <-  rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
}

exp_loc_true <- sample(1:p,ceiling(p*r_pq))
gamma  <- rep(0,p);
for(i in exp_loc_true){
  gamma[i] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
}

phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
phi_z <- rnorm(p_z,eta,1)*sample(c(1,-1),size=p_z,replace = T)



H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
Z <- matrix(rnorm(n*p_z,0,1),nrow=n,byrow=T)
E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigma_E)
X <- Z %*% Psi_z + H %*% Psi + E

E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
M_multi <- X %*% Theta + Z %*% Phi_z + H %*% Phi + E_tilde

e <- matrix(rnorm(n,0,1),nrow=n)
Y <- X %*% gamma + M_multi %*% beta + Z %*% phi_z + H %*% phi + e
```

start estimation using `HILAMA`

```R
require(HILAMA)
result_ls <- hilama_fast(X, M_multi, Y,Z,parallel=T,core_num=4)


## screen + JST test + BH 
nie_mat_hat <- result_ls$nie_mat_hat
pvalue_Theta_mat <- result_ls$pvalue_Theta
pvalue_beta <- result_ls$pvalue_beta
nie_mat_true <- Theta %*% diag(beta)

screen_N <- ceiling(0.1*p*q)
result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat,pvalue_beta,screen_N=screen_N,
                                                  alpha=0.1)
pvalue_nie_screen_JST_mat <- result_screen_jst_ls$pvalue_screen_JST_mat
result_eval_screen_JST <- cal_eval_metric(nie_mat_true,nie_mat_hat,pvalue_nie_screen_JST_mat,sig_level=result_screen_jst_ls$pvalue_cutoff)
result_eval_screen_JST

# Get significant NIE and NDE
result_sig_ls <- get_sig_nie_nde(result_ls)
result_sig_ls$result_nie
result_sig_ls$result_nde
```

## References

[1] Sun, Y., Ma, L., & Xia, Y. (2023). A decorrelating and debiasing approach to simultaneous inference for high-dimensional confounded models. Journal of the American Statistical Association, 1-12.

