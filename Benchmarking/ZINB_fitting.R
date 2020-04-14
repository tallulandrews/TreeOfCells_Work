simulate_counts_from_ZINB <- function(nc,ng){
  ## set seed for mu, r , and d per gene ##
  mu_vec <- 10^(runif(ng,0,4))
  r_vec <- rgamma(n = ng,shape = 1,rate = 0.5)
  d_vec <- runif(ng,0,0.9)
  
  ## build simulation row by row using gene specific values of mu r and d ##
  counts <- t(sapply(1:ng,function(g){
    
    # Generate counts from negative binomial #
    gNB <- rnbinom(nc, size = r_vec[g], mu = mu_vec[g])
    
    # Generate technical dropouts by simulating from binomial; 1=dropout #
    gDO <- rbinom(nc, size = 1, prob = d_vec[g])
    
    # Generate zero-inflated counts #
    gcount <- gNB * (1-gDO)
    return(gcount)
  }))
  
  row.names(counts) <- paste('gene',1:ng,sep='')
  colnames(counts) <- paste('cell',1:nc,sep='')
  param_df <- data.frame(mu=mu_vec, r=r_vec, d=d_vec, N=rep(nc, length(mu_vec)))
  rownames(param_df) <- rownames(counts);
  return(list(mat=counts, params= param_df))
}

simulate_counts_from_NB <- function(nc,ng){
  ## set seed for mu, r , and d per gene ##
  mu_vec <- 10^(runif(ng,0,4))
  r_vec <- rgamma(n = ng,shape = 1,rate = 0.5)
  
  ## build simulation row by row using gene specific values of mu r and d ##
  counts <- t(sapply(1:ng,function(g){
    
    # Generate counts from negative binomial #
    gNB <- rnbinom(nc, size = r_vec[g], mu = mu_vec[g])
    
    gcount <- gNB
    return(gcount)
  }))
  
  row.names(counts) <- paste('gene',1:ng,sep='')
  colnames(counts) <- paste('cell',1:nc,sep='')
  return(list(mat=counts, params=data.frame(mu=mu_vec, r=r_vec, N=rep(nc, length(mu_vec)))))
}

library(TreeOfCells)

set.seed(1192)
sim <- simulate_counts_from_ZINB(3000, 50000)
fit <- TreeOfCells::fit_ZINB_to_matrix(sim$mat)


sim2 <- simulate_counts_from_NB(3000, 50000)
clean <- which(rowSums(sim2$mat) > 0)
sim2$mat <- sim2$mat[clean,]
sim2$params <- sim2$params[clean,]
fit2 <- TreeOfCells::fit_NB_to_matrix(sim2$mat)

err <- (fit2[,1:2] - sim2$params[,1:2])/sim2$params[,1:2]
err <- abs(err)*100;

sim_fit_out <- list(ZINB_sim=sim$params, ZINB_fit=fit, NB_sim=sim2$params, NB_fit=fit2)
saveRDS(sim_fit_out, "Fitting_Errors_nc3000.rds")


set.seed(1192)
sim <- simulate_counts_from_ZINB(300, 50000)
fit <- TreeOfCells::fit_ZINB_to_matrix(sim$mat)


sim2 <- simulate_counts_from_NB(300, 50000)
clean <- which(rowSums(sim2$mat) > 0)
sim2$mat <- sim2$mat[clean,]
sim2$params <- sim2$params[clean,]
fit2 <- TreeOfCells::fit_NB_to_matrix(sim2$mat)

err <- (fit2[,1:2] - sim2$params[,1:2])/sim2$params[,1:2]
err <- abs(err)*100;

sim_fit_out <- list(ZINB_sim=sim$params, ZINB_fit=fit, NB_sim=sim2$params, NB_fit=fit2)
saveRDS(sim_fit_out, "Fitting_Errors_nc300.rds")

set.seed(1192)
sim <- simulate_counts_from_ZINB(30, 50000)
fit <- TreeOfCells::fit_ZINB_to_matrix(sim$mat)


sim2 <- simulate_counts_from_NB(30, 50000)
clean <- which(rowSums(sim2$mat) > 0)
sim2$mat <- sim2$mat[clean,]
sim2$params <- sim2$params[clean,]
fit2 <- TreeOfCells::fit_NB_to_matrix(sim2$mat)

err <- (fit2[,1:2] - sim2$params[,1:2])/sim2$params[,1:2]
err <- abs(err)*100;

sim_fit_out <- list(ZINB_sim=sim$params, ZINB_fit=fit, NB_sim=sim2$params, NB_fit=fit2)
saveRDS(sim_fit_out, "Fitting_Errors_nc30.rds")


