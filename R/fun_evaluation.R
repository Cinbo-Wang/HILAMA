cal_eval_metric_pvalue <- function(Theta_true,Theta_hat,pvalue_mat,sig_level=0.05){

  loc_active <- which(abs(Theta_true) > 1e-5,arr.ind = T)
  loc_zero <- which(abs(Theta_true) <= 1e-5, arr.ind = T)

  # 仅计算active bias
  bias <- Reduce('+',
                 apply(loc_active,1,function(x){
                   abs(Theta_hat[x[1],x[2]] - Theta_true[x[1],x[2]])
                 }))
  power <- Reduce('+',
                  apply(loc_active,1,function(x){
                    as.numeric(pvalue_mat[x[1],x[2]] <= sig_level)
                  })) / nrow(loc_active)  # =也算显著的

  if(sum(pvalue_mat <= sig_level) ==0){
    fdr <- 1
  }else{
    fdr <- Reduce('+',apply(loc_zero, 1, function(x){
      as.numeric(pvalue_mat[x[1],x[2]]<=sig_level)
    }))/sum(pvalue_mat <= sig_level)
  }

  bias_mean <- bias/length(loc_active)
  eval_result <- c(fdr,power,bias_mean)
  names(eval_result) <- c('fdr','power','bias_mean')
  return(eval_result)
}

cal_eval_metric_size <- function(nie_mat_true,nie_mat_hat,nie_scaled_mat_hat,fdr_level=0.1){

  p <- nrow(nie_mat_true); q <- ncol(nie_mat_hat)

  nie_data <- data.frame(index = 1:(p*q),
                         nie_true = as.vector(nie_mat_true),
                         nie_hat = as.vector(nie_mat_hat),
                         nie_hat_scaled = as.vector(nie_scaled_mat_hat))
  loc_active <- which(nie_data$nie_true != 0)
  loc_zero <- which(nie_data$nie_true == 0)


  bias <- abs(nie_data$nie_true[loc_active] - nie_data$nie_hat[loc_active])
  bias_mean <- mean(bias)


  ## 评价方式2：控制fdr为0.1后，计算tpr
  fdr_vec <- c()
  nie_vec <- sort(unique(abs(nie_data$nie_hat_scaled)),decreasing = T)

  if(length(nie_vec)==1){
    # 所有元素都为0
    eval_result <- c(0,0)
  }else{
    for(i in 2:min(1000,length(nie_vec))){
      pre_P <-  sum(abs(nie_data$nie_hat_scaled) >= nie_vec[i])
      TP <- Reduce('+',sapply(loc_active,function(x){
        return(abs(nie_data$nie_hat_scaled[x])>=nie_vec[i])
      }))
      FP <- pre_P - TP
      fdr <- FP /  pre_P
      fdr_vec <- c(fdr_vec,fdr)
    }

    if(sum(fdr_vec<fdr_level)==0){
      # 当最小fdr都大于fdr_sig level, 输出最小fdr下，最小的threshold
      nie_thre_fdr <- min(nie_vec[which(fdr_vec == min(fdr_vec))+1])
    }else{
      nie_thre_fdr <- min(nie_vec[which(fdr_vec <= fdr_level) + 1])
    }

    pre_P <- sum(abs(nie_data$nie_hat_scaled) >= nie_thre_fdr)
    TP <- Reduce('+',sapply(loc_active,function(x){
      return(abs(nie_data$nie_hat_scaled[x])>= nie_thre_fdr)
    }))
    FP <- pre_P - FP

    fdr <- FP / pre_P
    tpr <- TP / length(loc_active)

    eval_result <- c(fdr,tpr)
  }

  names(eval_result) <- c('fdr','power','bias_mean')

  return(eval_result)
}


get_nie_pvalue_screen_JST <- function(pvalue_Theta_mat,pvalue_beta,screen_N=NULL,
                                      fdr_level=0.1){

  p <- nrow(pvalue_Theta_mat); q <- ncol(pvalue_Theta_mat)

  pvalue_nie_mat_min <- pvalue_nie_mat_JST <- matrix(0,p,q)
  for(i in 1:p){
    for(j in 1:q){
      pvalue_nie_mat_min[i,j] <- min(pvalue_Theta_mat[i,j],pvalue_beta[j])
      pvalue_nie_mat_JST[i,j] <- max(pvalue_Theta_mat[i,j],pvalue_beta[j])
    }
  }
  p_vec <- as.vector(pvalue_nie_mat_min)
  c <- max(p_vec[rank(p_vec,ties.method = 'min')<=screen_N]) # for ties, we take the minimal value

  screen_loc <- which(pvalue_nie_mat_min <= c,arr.ind = T)
  pvalue_screen_JST_mat <- matrix(1,p,q)
  for(it in 1:nrow(screen_loc)){
    row_i <- screen_loc[it,1]; col_j <- screen_loc[it,2]
    pvalue_screen_JST_mat[row_i,col_j] <- max(pvalue_Theta_mat[row_i,col_j],pvalue_beta[col_j])
  }

  # 仅在满足min(pa,pb)的pair上进行BH adjust
  fdr_vec <- c();
  p_vec <- sort(unique(as.vector(pvalue_screen_JST_mat)),decreasing = F)
  for(i in 1:min(1000,length(p_vec))){
    t <- p_vec[i]
    denom <- sum(pvalue_screen_JST_mat <= t)  # 改为包含等号！
    numer <- nrow(screen_loc)*t  # 仅考虑min pvalue < c的
    fdr_tmp <- numer / denom
    fdr_vec <- c(fdr_vec,fdr_tmp)
    if(t > 0.1) break
  }

  if(min(fdr_vec)>fdr_level){
    warning(paste0('Minimal fdr is higher than the nominal fdr level:',fdr_level,'. The minimal fdr is',min(fdr_vec)))
    sig_level <- p_vec[which.min(fdr_vec)]
  }else{
    sig_level <- max((p_vec[1:length(fdr_vec)])[fdr_vec <= fdr_level])
  }

  return(list(
    pvalue_screen_JST_mat = pvalue_screen_JST_mat,
    pvalue_cutoff = sig_level
  ))

}


