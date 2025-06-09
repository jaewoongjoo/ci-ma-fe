packages <- c("MASS", "stats", "base", "mvtnorm", "foreach", "doParallel", "Matrix", "matrixcalc", "expm", "rstudioapi", "latex2exp")

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
  beta_k <- list(0);  beta_k_boot <- list(0)
  U_stack <-array(0, c(Nsim, sum(unlist(d)), dim+1))

  #target
  hat_bet <- array(0, c(Nsim, dim+1)); hat_bet_exact <- hat_bet; 
  hat_Phi_boot <- array(0, c(Nsim, dim+1, dim+1)); hat_Phi_exact_boot <- hat_Phi_boot; 
  bet_cover_boot <- rep(0,dim); bet_cover_exact_boot <- rep(0,dim); 

  for(kk in 1:num_study){
   d[[kk]] <- length(A[[kk]])
   beta_k[[kk]] <- rep(0, d[[kk]])
   S[[kk]] <- matrix(0, d[[kk]], d[[kk]])
   x_red_ext[[kk]] <- array(0,c(Num_ext, d[[kk]])) 
   x_study[[kk]] <- array(0, c(study_size[kk],dim+1)); 
   x_study_red[[kk]] <- array(0, c(study_size[kk],d[[kk]])); 
   y_study[[kk]] <- rep(0, study_size[kk])
  }

  for(i in 1:Nsim){
   Sigma2_x = sigma2_x * (matrix(x_corr, nrow=dim, ncol=dim) + diag(dim) * (1 - x_corr))
   Sigma2_x_intercept = matrix(0, nrow=dim+1, ncol=dim+1)
   Sigma2_x_intercept[2:(dim+1), 2:(dim+1)] = Sigma2_x
    
   x_ext <- MASS::mvrnorm(n=Num_ext, mu=mu_x, Sigma=Sigma2_x)
   E <- (1/Num_ext)*t(cbind(1,x_ext))%*%cbind(1,x_ext)
   E_exact <- Sigma2_x_intercept + c(1,mu_x)%*%t(c(1,mu_x));  
   U_exact <- list(0)
    
   y_mc <- list(0)
   x_mc <- list(0)  
   x_mc_red <- list(0)
    
   num_mc_int <- 50000
   for(kk in 1:num_study){
    x_mc[[kk]] <- MASS::mvrnorm(n=num_mc_int, mu=mu_x, Sigma2_x)
    y_mc[[kk]] <- x_mc[[kk]]%*%true_bet + sqrt(sigma2_y)*rt(n=num_mc_int, df= df_y) 
    x_mc_red[[kk]] <- x_mc[[kk]][,A[[kk]]]
    reg_mc <- lm(y_mc[[kk]]~x_mc_red[[kk]])
    S_exact[[kk]] <- vcov(reg_mc)
   }
 
   for(kk in 1:num_study){
    #construct study statistics: \beta^{(k)}'s and S_k's.
    x_study[[kk]]  <- MASS::mvrnorm(n=study_size[kk], mu=mu_x, Sigma= Sigma2_x)
    y_study[[kk]]  <- x_study[[kk]]%*%true_bet + sqrt(sigma2_y)*rt(n=study_size[kk], df= df_y)  
    
    x_study_red[[kk]] <- x_study[[kk]][,A[[kk]]]
    reg <- lm(y_study[[kk]]~x_study_red[[kk]])
    beta_k[[kk]] <- as.vector(reg$coefficients)
    names(beta_k[[kk]]) <- c("(Intercept)",paste("X", as.character(A[[kk]]), sep=""))
    S[[kk]] <- vcov(reg)
    
    #use external datapoints to construct U_k's
    x_red_ext[[kk]] <- x_ext[, A[[kk]]]
    U[[kk]] <- solve(E[c(1,A[[kk]]+1),c(1,A[[kk]]+1) ] )%*% (E[c(1,A[[kk]]+1), ] )
    U_exact[[kk]] <- solve(E_exact[c(1,A[[kk]]+1),c(1,A[[kk]]+1)])%*% (E_exact[c(1,A[[kk]]+1), ] )
   }

   U_stack <- do.call(rbind, U)
   U_exact_stack <- do.call(rbind, U_exact)
  
   Z <- unlist(beta_k)
   V <- as.matrix(Matrix::bdiag(S))
   V_exact <- as.matrix(Matrix::bdiag(S_exact))

   #final: estimate of beta (including intercept)
   hat_bet[i,] <- as.vector(solve(t(U_stack)%*%solve(V)%*%U_stack)%*%t(U_stack)%*%solve(V)%*%Z)
   hat_bet_exact[i,] <-  as.vector(solve(t(U_exact_stack)%*%solve(V_exact)%*%U_exact_stack)%*%t(U_exact_stack)%*%solve(V_exact)%*%Z)

  #Bootstrap
  chunklength <- sum(unlist(d)+1)
  bet_boot_k <- rep(0,dim+1)
  bet_boot_exact_k <- rep(0,dim+1)
  Var_bet_boot <- array(0, c(dim+1,dim+1))
  Var_bet_exact_boot <- array(0, c(dim+1,dim+1))
  Z_boot <- unlist(beta_k)
  
  #this is a design matrix for our new linear model V^{-1/2}Z=U'b+e, e~(0,I)
  Z_prime  <- sqrtm(solve(V))%*%Z_boot
  U_prime  <- as.matrix(sqrtm(solve(V))%*%U_stack)
  Z_prime_exact  <- sqrtm(solve(V_exact))%*%Z_boot
  U_prime_exact  <- as.matrix(sqrtm(solve(V_exact))%*%U_exact_stack)

  # Initialize the mean and covariance matrix
  mean_bet_boot <- rep(0, dim+1)
  mean_bet_exact_boot <- rep(0, dim+1)
  
  count <- 0  
  
  for(bb in 1:num_boot){
    d_boot_index <- sample(1:chunklength, chunklength, prob= rep(1/chunklength, chunklength), replace=TRUE)
    U_prime_boot <- U_prime[d_boot_index,] 
    U_prime_exact_boot <- U_prime_exact[d_boot_index,] 
    Z_prime_boot <- Z_prime[d_boot_index,]
    Z_prime_exact_boot <- Z_prime_exact[d_boot_index,]
    
    singular_flag <- FALSE
    
    bet_boot_k <- tryCatch({
      solve(t(U_prime_boot) %*% U_prime_boot) %*% t(U_prime_boot) %*% Z_prime_boot
    }, error = function(e) {
      singular_flag <<- TRUE
      return(NULL)
    })
    
    bet_boot_exact_k <- tryCatch({
      solve(t(U_prime_exact_boot) %*% U_prime_exact_boot) %*% t(U_prime_exact_boot) %*% Z_prime_exact_boot
    }, error = function(e) {
      singular_flag <<- TRUE
      return(NULL)
    })

    if (singular_flag) {
      next
    }
    
    count <- count + 1  

    delta_boot <- bet_boot_k - mean_bet_boot
    mean_bet_boot <- mean_bet_boot + delta_boot / count
    Var_bet_boot <- Var_bet_boot + (delta_boot %*% t(bet_boot_k - mean_bet_boot))
    
    delta_exact_boot <- bet_boot_exact_k - mean_bet_exact_boot
    mean_bet_exact_boot <- mean_bet_exact_boot + delta_exact_boot / count
    Var_bet_exact_boot <- Var_bet_exact_boot + (delta_exact_boot %*% t(bet_boot_exact_k - mean_bet_exact_boot))
  }
  
  Var_bet_boot <- Var_bet_boot/(count-1)
  Var_bet_exact_boot <- Var_bet_exact_boot/(count-1)
  
  hat_Phi_boot[i,,] <- Var_bet_boot
  hat_Phi_exact_boot[i,,] <- Var_bet_exact_boot
  
  for(ii in 1:dim){
    if( (true_bet[ii] > hat_bet[i,ii+1] - 1.96*sqrt(diag(hat_Phi_boot[i,,]))[ii+1]) &&  (true_bet[ii] < hat_bet[i,ii+1] + 1.96*sqrt(diag(hat_Phi_boot[i,,]))[ii+1]) ){
      bet_cover_boot[ii] <- bet_cover_boot[ii]+1 }
    if( (true_bet[ii] > hat_bet_exact[i,ii+1] - 1.96*sqrt(diag(hat_Phi_exact_boot[i,,]))[ii+1]) &&  (true_bet[ii] < hat_bet_exact[i,ii+1] + 1.96*sqrt(diag(hat_Phi_exact_boot[i,,]))[ii+1]) ){
      bet_cover_exact_boot[ii]<- bet_cover_exact_boot[ii]+1 }
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
 
 rownames(output_table_exact)= rownames(output_table); 
 total_output<- list(output_table,output_table_exact,hat_bet)
 names(total_output) <- c("Real", "Oracle", "Estimates")

 return(total_output)
}

cl <- 4
registerDoParallel(cl)


###1. create figure S-6
###1-1. create the left subfigure for figure S-6 (p=6)

dim <- 6
true_bet <- c(-4, 2, 0, -2, 4, -2)
mu_x <- c(10, -20, 30, -10, 20, -30)
sigma2_x <- 10
sigma2_y <- 3
df_y <- 4.1
num_study <- 10
sample_size <- c(300,300,300,400,400,400,500,500,600,600)

A <- list(0)
A[[1]] = c(1,2)
A[[2]] = c(3,4)
A[[3]] = c(5,6)
A[[4]] = c(1,6)
A[[5]] = c(2,3,4,5)
A[[6]] = c(1,2)
A[[7]] = c(3,4)
A[[8]] = c(5,6)
A[[9]] = c(1,6)
A[[10]] = c(2,3,4,5)

results_50000 <- foreach(rho=c(0.3, 0.6, 0.8, 0.9), .combine = list, .errorhandling = "remove")%dopar%{
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

rho_0.3 <-  results_50000[[1]][[1]][[1]][[2]];
rho_0.6 <-  results_50000[[1]][[1]][[2]];
rho_0.8 <-  results_50000[[1]][[2]];
rho_0.9 <-  results_50000[[2]];

bet_estim <- rbind(rho_0.3$Real[1,], rho_0.3$Real[2,], rho_0.6$Real[2,], rho_0.8$Real[2,], rho_0.9$Real[2,])

bet_var <- rbind(rho_0.3$Real[4,],  rho_0.6$Real[4,],  rho_0.8$Real[4,],  rho_0.9$Real[4,],
                 rho_0.3$Oracle[4,],rho_0.6$Oracle[4,],rho_0.8$Oracle[4,],rho_0.9$Oracle[4,])

bet_mse <- rbind(rho_0.3$Real[5,],  rho_0.6$Real[5,],  rho_0.8$Real[5,],  rho_0.9$Real[5,],
                 rho_0.3$Oracle[5,],rho_0.6$Oracle[5,],rho_0.8$Oracle[5,],rho_0.9$Oracle[5,])

bet_cov <- rbind(rho_0.3$Real[6,],  rho_0.6$Real[6,],  rho_0.8$Real[6,],  rho_0.9$Real[6,],
                 rho_0.3$Oracle[6,],rho_0.6$Oracle[6,],rho_0.8$Oracle[6,],rho_0.9$Oracle[6,])

true_bet <- rho_0.3$Real[1,]
num_beta <- length(true_bet)
beta_expressions <- do.call(c, lapply(1:num_beta, function(i) as.expression(substitute(beta[i], list(i = i)))))

rownames(bet_estim) <- c("True Value", "rho0.3", "rho0.6", "rho0.8", "rho0.9")
colnames(bet_var) <- colnames(bet_estim); colnames(bet_mse)<- colnames(bet_estim);
rownames(bet_var) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_cov) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_mse) <- rownames(bet_var); 
colnames(bet_cov) <- paste0("beta_", 1:num_beta)

colors <- c("red","black")
pch_values <- c(0,1,2,5)
pch_values_filled <- c(15, 16, 17, 18) 

pdf_filename_panel3 <- paste0("panel3-k10-p6.pdf")
pdf(pdf_filename_panel3 , width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )

coverage_real <- bet_cov[1:4,]
coverage_oracle <- bet_cov[5:8,]

plot(NULL, xlim=c(1, num_beta), ylim=c(0.88, 1.00), xlab="", ylab="", main="", xaxt='n', cex=1.4,  cex.lab=1.4)
axis(1, at = 1:num_beta, cex.axis=1.4, labels = beta_expressions)
abline(h=0.95, lty=2, col="red")

for (i in 1:4){
  points(1:num_beta, as.numeric(coverage_real[i, ]), col=colors[1],   pch=pch_values_filled[i], cex=1)
  points(1:num_beta, as.numeric(coverage_oracle[i, ]), col=colors[2], pch=pch_values[i], cex=1.6)
}

legend("topleft", 
       legend=c(TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')),
       pch=pch_values)
dev.off()


###1-2. create the middle subfigure for figure S-6 (p=5)
dim <- 5
true_bet <- c(-4, 2, 0, -2, 4)
mu_x <- c(10, -20, 30, -10, 20)

A <- list(0)
A[[1]] = c(1,2)
A[[2]] = c(3)
A[[3]] = c(4,5)
A[[4]] = c(1)
A[[5]] = c(2,3)
A[[6]] = c(4,5)
A[[7]] = c(1,2,3)
A[[8]] = c(4,5)
A[[9]] = c(1,2)
A[[10]] = c(3,4,5)

results_50000 <- foreach(rho=c(0.3, 0.6, 0.8, 0.9), .combine = list, .errorhandling = "remove")%dopar%{
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

rho_0.3 <-  results_50000[[1]][[1]][[1]][[2]];
rho_0.6 <-  results_50000[[1]][[1]][[2]];
rho_0.8 <-  results_50000[[1]][[2]];
rho_0.9 <-  results_50000[[2]];

bet_estim <- rbind(rho_0.3$Real[1,], rho_0.3$Real[2,], rho_0.6$Real[2,], rho_0.8$Real[2,], rho_0.9$Real[2,])

bet_var <- rbind(rho_0.3$Real[4,],  rho_0.6$Real[4,],  rho_0.8$Real[4,],  rho_0.9$Real[4,],
                 rho_0.3$Oracle[4,],rho_0.6$Oracle[4,],rho_0.8$Oracle[4,],rho_0.9$Oracle[4,])

bet_mse <- rbind(rho_0.3$Real[5,],  rho_0.6$Real[5,],  rho_0.8$Real[5,],  rho_0.9$Real[5,],
                 rho_0.3$Oracle[5,],rho_0.6$Oracle[5,],rho_0.8$Oracle[5,],rho_0.9$Oracle[5,])

bet_cov <- rbind(rho_0.3$Real[6,],  rho_0.6$Real[6,],  rho_0.8$Real[6,],  rho_0.9$Real[6,],
                 rho_0.3$Oracle[6,],rho_0.6$Oracle[6,],rho_0.8$Oracle[6,],rho_0.9$Oracle[6,])

true_bet <- rho_0.3$Real[1,]
num_beta <- length(true_bet)
beta_expressions <- do.call(c, lapply(1:num_beta, function(i) as.expression(substitute(beta[i], list(i = i)))))

rownames(bet_estim) <- c("True Value", "rho0.3", "rho0.6", "rho0.8", "rho0.9")
colnames(bet_var) <- colnames(bet_estim); colnames(bet_mse)<- colnames(bet_estim);
rownames(bet_var) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_cov) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_mse) <- rownames(bet_var); 
colnames(bet_cov) <- paste0("beta_", 1:num_beta)

colors <- c("red","black")
pch_values <- c(0,1,2,5)
pch_values_filled <- c(15, 16, 17, 18) 

pdf_filename_panel3 <- paste0("panel3-k10-p5.pdf")
pdf(pdf_filename_panel3 , width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )

coverage_real <- bet_cov[1:4,]
coverage_oracle <- bet_cov[5:8,]

plot(NULL, xlim=c(1, num_beta), ylim=c(0.88, 1.00), xlab="", ylab="", main="", xaxt='n', cex=1.4,  cex.lab=1.4)
axis(1, at = 1:num_beta, cex.axis=1.4, labels = beta_expressions)
abline(h=0.95, lty=2, col="red")

for (i in 1:4){
  points(1:num_beta, as.numeric(coverage_real[i, ]), col=colors[1],   pch=pch_values_filled[i], cex=1)
  points(1:num_beta, as.numeric(coverage_oracle[i, ]), col=colors[2], pch=pch_values[i], cex=1.6)
}

legend("topleft", 
       legend=c(TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')),
       pch=pch_values)
dev.off()

###1-3. create the right subfigure for figure S-6 (p=4)
dim <- 4
true_bet <- c(-4, 2, 0, -2)
mu_x <- c(10, -20, 30, -10)

A <- list(0)
A[[1]] = c(1,2)
A[[2]] = c(3,4)
A[[3]] = c(1,4)
A[[4]] = c(2,3)
A[[5]] = c(1,2)
A[[6]] = c(3,4)
A[[7]] = c(1,4)
A[[8]] = c(2,3)
A[[9]] = c(1,2)
A[[10]] = c(3,4)

results_50000 <- foreach(rho=c(0.3, 0.6, 0.8, 0.9), .combine = list, .errorhandling = "remove")%dopar%{
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

rho_0.3 <-  results_50000[[1]][[1]][[1]][[2]];
rho_0.6 <-  results_50000[[1]][[1]][[2]];
rho_0.8 <-  results_50000[[1]][[2]];
rho_0.9 <-  results_50000[[2]];

bet_estim <- rbind(rho_0.3$Real[1,], rho_0.3$Real[2,], rho_0.6$Real[2,], rho_0.8$Real[2,], rho_0.9$Real[2,])

bet_var <- rbind(rho_0.3$Real[4,],  rho_0.6$Real[4,],  rho_0.8$Real[4,],  rho_0.9$Real[4,],
                 rho_0.3$Oracle[4,],rho_0.6$Oracle[4,],rho_0.8$Oracle[4,],rho_0.9$Oracle[4,])

bet_mse <- rbind(rho_0.3$Real[5,],  rho_0.6$Real[5,],  rho_0.8$Real[5,],  rho_0.9$Real[5,],
                 rho_0.3$Oracle[5,],rho_0.6$Oracle[5,],rho_0.8$Oracle[5,],rho_0.9$Oracle[5,])

bet_cov <- rbind(rho_0.3$Real[6,],  rho_0.6$Real[6,],  rho_0.8$Real[6,],  rho_0.9$Real[6,],
                 rho_0.3$Oracle[6,],rho_0.6$Oracle[6,],rho_0.8$Oracle[6,],rho_0.9$Oracle[6,])

true_bet <- rho_0.3$Real[1,]
num_beta <- length(true_bet)
beta_expressions <- do.call(c, lapply(1:num_beta, function(i) as.expression(substitute(beta[i], list(i = i)))))

rownames(bet_estim) <- c("True Value", "rho0.3", "rho0.6", "rho0.8", "rho0.9")
colnames(bet_var) <- colnames(bet_estim); colnames(bet_mse)<- colnames(bet_estim);
rownames(bet_var) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_cov) <- c("rho0.3_Real", "rho0.6_Real", "rho0.8_Real", "rho0.9_Real","rho0.3_exact", "rho0.6_exact", "rho0.8_exact", "rho0.9_exact") ; 
rownames(bet_mse) <- rownames(bet_var); 
colnames(bet_cov) <- paste0("beta_", 1:num_beta)

colors <- c("red","black")
pch_values <- c(0,1,2,5)
pch_values_filled <- c(15, 16, 17, 18) 

pdf_filename_panel3 <- paste0("panel3-k10-p4.pdf")
pdf(pdf_filename_panel3 , width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )

coverage_real <- bet_cov[1:4,]
coverage_oracle <- bet_cov[5:8,]

plot(NULL, xlim=c(1, num_beta), ylim=c(0.88, 1.00), xlab="", ylab="", main="", xaxt='n', cex=1.4,  cex.lab=1.4)
axis(1, at = 1:num_beta, cex.axis=1.4, labels = beta_expressions)
abline(h=0.95, lty=2, col="red")

for (i in 1:4){
  points(1:num_beta, as.numeric(coverage_real[i, ]), col=colors[1],   pch=pch_values_filled[i], cex=1)
  points(1:num_beta, as.numeric(coverage_oracle[i, ]), col=colors[2], pch=pch_values[i], cex=1.6)
}

legend("topleft", 
       legend=c(TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')),
       pch=pch_values)
dev.off()