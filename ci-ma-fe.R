packages <- c("MASS", "stats", "base", "mvtnorm", "foreach", "doParallel", "Matrix", "matrixcalc", "expm", "rstudioapi")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) { install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org") }
invisible(lapply(packages, library, character.only = TRUE))

loc.current <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else if(Sys.getenv("RSTUDIO")=="1") {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}
code.dir <- loc.current()
setwd(code.dir)


ci_ma_fe_simul <- function(Nsim, dim, true_bet, Num_ext, num_study, study_size, A, mu_y, sigma2_y, df_y, mu_x, sigma2_x, x_corr, num_boot, seed_num){
  set.seed(seed_num)
  d <- list(0)
  S <- list(0); S_exact<- list(0); S_boot <- list(0)
  x_study <- list(0); x_study_red <- list(0)
  y_study <- list(0); y_study_ext <- array(0,c(Num_ext,num_study))
  x_ext <- array(0, c(Num_ext, dim)); x_red_ext  <- list(0)
  E <- array(0, c(dim+1, dim+1)); 
  U <- list(0); U_boot <- list(0)
  beta_k <- list(0);
  U_stack <-array(0, c(Nsim, sum(unlist(d)), dim+1))
  
  #target
  hat_bet <- array(0, c(Nsim, dim+1)); hat_bet_exact <- hat_bet; hat_bet_rootn<- hat_bet 
  hat_Phi_boot <- array(0, c(Nsim, dim+1, dim+1)); hat_Phi_exact_boot <- hat_Phi_boot; hat_Phi_rootn_boot <- hat_Phi_boot
  bet_cover_boot <- rep(0,dim); bet_cover_exact_boot <- rep(0,dim); bet_cover_rootn_boot <- rep(0,dim);
  
  r1 =  x_corr; r2 =r1 ; r3 = r2; r4 = r2; r5 = r2; r6 = r2; r7 = r2;
  Sigma2_x = sigma2_x*matrix(
    c(1, r1, r2, r3, r4, r5, r6, r7,
      r1, 1, r1, r2, r3, r4, r5, r6,
      r2, r1, 1, r1, r2, r3, r4, r5,
      r3, r2, r1, 1, r1, r2, r3, r4,
      r4, r3, r2, r1, 1, r1, r2, r3,
      r5, r4, r3, r2, r1, 1, r1, r2,
      r6, r5, r4, r3, r2, r1, 1, r1,
      r7, r6, r5, r4, r3, r2, r1, 1),
    nrow=dim,  ncol=dim) 
  
  Sigma2_x_intercept = sigma2_x*matrix(
    c(0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  1, r1, r2, r3, r4, r5, r6, r7,
      0, r1, 1, r1, r2, r3, r4, r5, r6,
      0, r2, r1, 1, r1, r2, r3, r4, r5,
      0, r3, r2, r1, 1, r1, r2, r3, r4,
      0, r4, r3, r2, r1, 1, r1, r2, r3,
      0, r5, r4, r3, r2, r1, 1, r1, r2,
      0, r6, r5, r4, r3, r2, r1, 1, r1,
      0, r7, r6, r5, r4, r3, r2, r1, 1),
    nrow=dim+1,  ncol=dim+1) 
  
  y_mc <- list(0)
  x_mc <- list(0)  
  x_mc_red <- list(0)
  
  U_exact <- list(0)
  E_exact <- Sigma2_x_intercept + c(1,mu_x)%*%t(c(1,mu_x));  #This is the exact matrix of E(X'X)
  num_mc_int <- 50000
  
  for(kk in 1:num_study){
    d[[kk]] <- length(A[[kk]])
    beta_k[[kk]] <- rep(0, d[[kk]])
    S[[kk]] <- matrix(0, d[[kk]], d[[kk]])
    x_red_ext[[kk]] <- array(0,c(Num_ext, d[[kk]])) 
    x_study[[kk]] <- array(0, c(study_size[kk],dim+1)); 
    x_study_red[[kk]] <- array(0, c(study_size[kk],d[[kk]])); 
    y_study[[kk]] <- rep(0, study_size[kk])
    
    x_mc[[kk]] <- MASS::mvrnorm(n=num_mc_int, mu=mu_x, Sigma2_x)
    y_mc[[kk]] <- x_mc[[kk]]%*%true_bet + sqrt(sigma2_y)*rt(n=num_mc_int, df= df_y) 
    x_mc_red[[kk]] <- x_mc[[kk]][,A[[kk]]]
    reg_mc <- lm(y_mc[[kk]]~x_mc_red[[kk]])
    S_exact[[kk]] <- vcov(reg_mc)
    U_exact[[kk]] <- solve(E_exact[c(1,A[[kk]]+1),c(1,A[[kk]]+1)])%*% (E_exact[c(1,A[[kk]]+1),])
    
    x_study[[kk]]  <- MASS::mvrnorm(n=study_size[kk], mu=mu_x, Sigma= Sigma2_x)
    y_study[[kk]]  <- x_study[[kk]]%*%true_bet + sqrt(sigma2_y)*rt(n=study_size[kk], df= df_y)  
    
    x_study_red[[kk]] <- x_study[[kk]][,A[[kk]]]
    reg <- lm(y_study[[kk]]~x_study_red[[kk]])
    beta_k[[kk]] <- as.vector(reg$coefficients)
    names(beta_k[[kk]]) <- c("(Intercept)",paste("X", as.character(A[[kk]]), sep=""))
    S[[kk]] <- vcov(reg)
  }
  
  U_exact_stack <- do.call(rbind, U_exact)
  V_exact <- as.matrix(Matrix::bdiag(S_exact))
  U_prime_exact  <- as.matrix(expm::sqrtm(solve(V_exact))%*%U_exact_stack)
  
  for(i in 1:Nsim){
    set.seed(seed_num + i)
    x_ext <- MASS::mvrnorm(n=Num_ext, mu=mu_x, Sigma=Sigma2_x)
    E <- (1/Num_ext)*t(cbind(1,x_ext))%*%cbind(1,x_ext)
    
    for(kk in 1:num_study){
      #use external datapoints to construct U_k's
      x_red_ext[[kk]] <- x_ext[, A[[kk]]]
      U[[kk]] <- solve(E[c(1,A[[kk]]+1),c(1,A[[kk]]+1) ] )%*% (E[c(1,A[[kk]]+1), ] )
    }
    
    U_stack <- do.call(rbind, U)
    Z <- unlist(beta_k)
    V <- as.matrix(Matrix::bdiag(S))
    
    hat_bet[i,] <- as.vector(solve(t(U_stack)%*%solve(V)%*%U_stack)%*%t(U_stack)%*%solve(V)%*%Z)
    hat_bet_exact[i,] <-  as.vector(solve(t(U_exact_stack)%*%solve(V_exact)%*%U_exact_stack)%*%t(U_exact_stack)%*%solve(V_exact)%*%Z)
    hat_bet_rootn[i,] <-  as.vector(solve(t(U_stack)%*%solve(V_exact)%*%U_stack)%*%t(U_stack)%*%solve(V_exact)%*%Z)
    
    # Bootstrap
    bet_boot_k <- rep(0,dim+1)
    bet_boot_rootn_k <- rep(0,dim+1)
    bet_boot_exact_k <- rep(0,dim+1)
    
    Var_bet_boot <- array(0, c(dim+1,dim+1))
    Var_bet_rootn_boot <- array(0, c(dim+1,dim+1))
    Var_bet_exact_boot <- array(0, c(dim+1,dim+1))
    
    #this is a design matrix for our new linear model V^{-1/2}Z=U'b+e, e~(0,I)
    Z_prime  <- expm::sqrtm(solve(V))%*%Z
    Z_prime_rootn  <- expm::sqrtm(solve(V_exact))%*%Z
    Z_prime_exact  <- expm::sqrtm(solve(V_exact))%*%Z
    
    U_prime  <- as.matrix(expm::sqrtm(solve(V))%*%U_stack)
    U_prime_rootn  <- as.matrix(expm::sqrtm(solve(V_exact))%*%U_stack)
    
    # Initialize the mean and covariance matrix
    mean_bet_boot <- rep(0, dim+1)
    mean_bet_exact_boot <- rep(0, dim+1)
    mean_bet_rootn_boot <- rep(0, dim+1)
    
    Z_boot <- unlist(beta_k)
    
    for(bb in 1:num_boot){
      chunklength <- sum(unlist(d)+1)
      d_boot_index <- sample(1:chunklength, chunklength, replace=TRUE)
      
      U_prime_boot <- U_prime[d_boot_index,] 
      Z_prime_boot <- Z_prime[d_boot_index,]
      U_prime_exact_boot <- U_prime_exact[d_boot_index,] 
      Z_prime_exact_boot <- Z_prime_exact[d_boot_index,]
      
      bet_boot_k <- solve(t(U_prime[d_boot_index,] )%*%U_prime[d_boot_index,] )%*%t(U_prime[d_boot_index,] )%*%Z_prime[d_boot_index,]
      bet_boot_rootn_k <- solve(t(U_prime_rootn[d_boot_index,])%*%U_prime_rootn[d_boot_index,])%*%t(U_prime_rootn[d_boot_index,])%*%Z_prime_rootn[d_boot_index,]
      bet_boot_exact_k <- solve(t(U_prime_exact[d_boot_index,])%*%U_prime_exact[d_boot_index,])%*%t(U_prime_exact[d_boot_index,]) %*% Z_prime_exact[d_boot_index,]
      
      delta_boot <- bet_boot_k - mean_bet_boot
      mean_bet_boot <- mean_bet_boot + delta_boot / bb
      Var_bet_boot <- Var_bet_boot + (delta_boot %*% t(bet_boot_k - mean_bet_boot))
      
      delta_rootn_boot <- bet_boot_rootn_k - mean_bet_rootn_boot
      mean_bet_rootn_boot <- mean_bet_rootn_boot + delta_rootn_boot / bb
      Var_bet_rootn_boot <- Var_bet_rootn_boot + (delta_rootn_boot %*% t(bet_boot_rootn_k - mean_bet_rootn_boot))
      
      delta_exact_boot <- bet_boot_exact_k - mean_bet_exact_boot
      mean_bet_exact_boot <- mean_bet_exact_boot + delta_exact_boot / bb
      Var_bet_exact_boot <- Var_bet_exact_boot + (delta_exact_boot %*% t(bet_boot_exact_k - mean_bet_exact_boot))
    }
    Var_bet_boot <- Var_bet_boot/(num_boot-1)
    Var_bet_rootn_boot <- Var_bet_rootn_boot/(num_boot-1)
    Var_bet_exact_boot <- Var_bet_exact_boot/(num_boot-1)
    
    hat_Phi_boot[i,,] <- Var_bet_boot
    hat_Phi_rootn_boot[i,,] <- Var_bet_rootn_boot
    hat_Phi_exact_boot[i,,] <- Var_bet_exact_boot
    
    for(ii in 1:dim){
      if( (true_bet[ii] > hat_bet[i,ii+1] - 1.96*sqrt(diag(hat_Phi_boot[i,,]))[ii+1]) &&  (true_bet[ii] < hat_bet[i,ii+1] + 1.96*sqrt(diag(hat_Phi_boot[i,,]))[ii+1]) ){
        bet_cover_boot[ii] <- bet_cover_boot[ii]+1 }
      if( (true_bet[ii] > hat_bet_exact[i,ii+1] - 1.96*sqrt(diag(hat_Phi_exact_boot[i,,]))[ii+1]) &&  (true_bet[ii] < hat_bet_exact[i,ii+1] + 1.96*sqrt(diag(hat_Phi_exact_boot[i,,]))[ii+1]) ){
        bet_cover_exact_boot[ii]<- bet_cover_exact_boot[ii]+1 }
      if( (true_bet[ii] > hat_bet_rootn[i,ii+1] - 1.96*sqrt(diag(hat_Phi_rootn_boot[i,,]))[ii+1]) &&  (true_bet[ii] < hat_bet_rootn[i,ii+1] + 1.96*sqrt(diag(hat_Phi_rootn_boot[i,,]))[ii+1]) ){
        bet_cover_rootn_boot[ii]<- bet_cover_rootn_boot[ii]+1 }
    }
  }
  
  output_table <- rbind(true_bet, #true
                        apply(hat_bet, 2,mean)[-1], #estim
                        true_bet-apply(hat_bet, 2, mean)[-1], #bias
                        (apply(hat_bet, 2,sd)^2)[-1],#var
                        (true_bet - apply(hat_bet, 2,mean)[-1])^2+(apply(hat_bet, 2,sd)[-1])^2,#mse
                        bet_cover_boot/Nsim) #coverage
  
  rownames(output_table) <- c("True Value", "Estimate", "Bias", "Variance", "MSE", "Coverage")
  
  output_table_exact <- rbind(true_bet, #true
                              apply(hat_bet_exact, 2,mean)[-1], #estim
                              true_bet-apply(hat_bet_exact, 2, mean)[-1], #bias
                              (apply(hat_bet_exact, 2,sd)^2)[-1],#var
                              (true_bet-apply(hat_bet_exact, 2,mean)[-1])^2+(apply(hat_bet_exact, 2,sd)[-1])^2,#mse
                              bet_cover_exact_boot/Nsim) #coverage
  
  output_table_rootn <- rbind(true_bet, #true
                              apply(hat_bet_rootn, 2,mean)[-1], #estim
                              true_bet-apply(hat_bet_rootn, 2, mean)[-1], #bias
                              (apply(hat_bet_rootn, 2,sd)^2)[-1],#var
                              (true_bet-apply(hat_bet_rootn, 2,mean)[-1])^2+(apply(hat_bet_rootn, 2,sd)[-1])^2,#mse
                              bet_cover_rootn_boot/Nsim) #coverage
  
  rownames(output_table_exact)= rownames(output_table); rownames(output_table_rootn)= rownames(output_table);
  total_output<- list(output_table,output_table_exact,output_table_rootn,hat_bet)
  names(total_output) <- c("Real", "Oracle", "Rootn", "Estimates")
  
  return(total_output)
}

cl <- detectCores()
registerDoParallel(cl)

dim <- 8
true_bet <- c(-4, 2, 0, -2, 4, -2, 0, 2)
mu_x <- c(10, -20, 30, -10, 20, -30, 10, -20)
sigma2_x <- 10
sigma2_y <- 3
df_y <- 4.1
num_study <- 20
sample_size <- c(rep(300, num_study/4), rep(400, num_study/4), rep(500, num_study/4), rep(600, num_study/4))

A <- vector("list", num_study)
base_patterns <- list(
  c(1,2),       
  c(3,4),      
  c(5,6),       
  c(7,8),      
  c(1,2,7,8),   
  c(3,4,5,6)   
)

for (i in 1:num_study) {
  A[[i]] <- base_patterns[[((i - 1) %% length(base_patterns)) + 1]]
}

results_50000 <- foreach(rho = seq(from=0, to=0.9, by=0.1), .combine = list, .packages = c("MASS", "Matrix", "expm")) %dopar% {
  ci_ma_fe_simul(
    Nsim = 5000,
    dim = dim,
    true_bet = true_bet,
    Num_ext = 50000,
    num_study = num_study,
    study_size = sample_size,
    A = A,
    mu_y = 5,
    sigma2_y = sigma2_y,
    df_y = df_y,
    mu_x = mu_x,
    sigma2_x = sigma2_x,
    x_corr = rho,
    num_boot = 1000,
    seed_num = 12345678
  )
}

results_10000 <- foreach(rho = seq(from=0, to=0.9, by=0.1), .combine = list, .packages = c("MASS", "Matrix", "expm")) %dopar% {
  ci_ma_fe_simul(
    Nsim = 5000,
    dim = dim,
    true_bet = true_bet,
    Num_ext = 10000,
    num_study = num_study,
    study_size = sample_size,
    A = A,
    mu_y = 5,
    sigma2_y = sigma2_y,
    df_y = df_y,
    mu_x = mu_x,
    sigma2_x = sigma2_x,
    x_corr = rho,
    num_boot = 1000,
    seed_num = 12345678
  )
}

save.image("gendata-ci-ma-fe.Rdata")