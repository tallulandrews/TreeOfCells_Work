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
  return(list(mat=counts, params=data.frame(mu=mu_vec, r=r_vec, d=d_vec, N=rep(nc, length(mu_vec)))))
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


sim <- simulate_counts_from_ZINB(3000, 30000)
fit <- TreeOfCells::fit_ZINB_to_matrix(sim$mat)


sim2 <- simulate_counts_from_NB(3000, 30000)
fit2 <- TreeOfCells::fit_NB_to_matrix(sim$mat)
