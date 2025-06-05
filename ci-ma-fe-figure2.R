packages <- c("expm")

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

url <- "https://osf.io/nb3y2/download"
download.file(url, destfile = "2013.csv", mode = "wb")
data <- read.csv("2013.csv")
rawdata <- data[-which(is.na(data$SMK_STAT)),]
attach(rawdata)

set.seed(12345678)
#variables
ggtp <- as.numeric(GAMMA_GTP)
map <- BP_LWST + 1/3*(BP_HIGH - BP_LWST)
bmi <- as.numeric(WEIGHT/((HEIGHT/100)^2))
drinking <- as.numeric(DRK_YN=="Y")
smoking <- as.numeric(SMK_STAT)
blds <- as.numeric(BLDS)
tc <- as.numeric(TOT_CHOLE)

#construct data frame
data <- data.frame(ggtp, map, bmi, drinking, smoking, blds, tc)
mydata <- data[sample(nrow(data), replace=FALSE),]
#form external data
N <- length(mydata$bmi)/2
ext_data <- mydata[(N+1):(2*N),]

#settings for independent studies
k <- 20
p <- 6
chunk <- split(1:N, ceiling(seq_along(1:N)/(N/k)))
num_study <- rep(N/k, k)

#indices, used to form X^{(k)}
A <- list(0)
A[[1]] = c(1); A[[2]] = c(1);
A[[3]] = c(1,2)
A[[4]] = c(1,3)
A[[5]] = c(1,4)
A[[6]] = c(1,5,6)
for(kk in 7:10) {   A[[kk]] = c(2,3)   }
for(kk in 11:14){   A[[kk]] = c(4,5,6)   }
for(kk in 15:18){   A[[kk]] = c(2,3,4)   }
for(kk in 19:20){   A[[kk]] = c(5,6)   }

#placeholders for constructing U and Z
d <- list(0)
S <- list(0); S_boot <- list(0)
x_study <- list(0); x_study_red <- list(0)
y_study <- list(0); y_study_ext <- array(0,c(N,k))
x_ext <- array(0, c(N,p)); x_red_ext  <- list(0)
E <- array(0, c(p+1, p+1)); 
U <- list(0); U_boot <- list(0)
beta_k <- list(0);  
for(kk in 1:k){
  d[[kk]] <- length(A[[kk]])
  beta_k[[kk]] <- rep(0, d[[kk]])
  S[[kk]] <- matrix(0, d[[kk]], d[[kk]])
  x_red_ext[[kk]] <- array(0,c(N,d[[kk]])) 
  x_study[[kk]] <- array(0, c(num_study[kk],p+1)); 
  x_study_red[[kk]] <- array(0, c(num_study[kk],d[[kk]])); 
  y_study[[k]] <- rep(0, num_study[kk])
}
U_stack <-array(0,  c(sum(unlist(d)), p+1))

#target
hat_bet <- array(0, p+1); beta_minus_k<- array(0, c(p+1, k))

w <- rep(0,k)
x_ext <- as.matrix(ext_data[,-1])
E <- (1/N)*t(cbind(1,x_ext))%*%cbind(1,x_ext)

for(kk in 1:k){
  #use external datapoints to construct U_k's
  x_red_ext[[kk]] <- cbind(1,x_ext[, A[[kk]] ])
  U[[kk]] <- solve(E[c(1,A[[kk]]+1),c(1,A[[kk]]+1) ] )%*% (E[c(1,A[[kk]]+1), ] )
  #construct study statistics: \beta^{(k)}'s and S_k's.
  x_study[[kk]] <- mydata[chunk[[kk]],-1]
  y_study[[kk]] <- mydata[chunk[[kk]],1]
  x_study_red[[kk]] <- x_study[[kk]][,(A[[kk]])]
  reg <- lm(y_study[[kk]] ~ as.matrix(x_study_red[[kk]]))
  beta_k[[kk]] <- as.vector(reg$coefficients)
  S[[kk]] <- vcov(reg)
}

beta_minus_k_i <- array(0, c(p+1,sum(unlist(d)+1)))
Var_bet_boot <- array(0, c(p+1,p+1));
U_stack <- do.call(rbind, U)
Z <- unlist(beta_k)

#construct the block diagonal matrix V=diag(S_1,\ldots, S_K)
V <- as.matrix(Matrix::bdiag(S))

#final: estimate of beta (including intercept)
hat_bet <- as.vector(solve(t(U_stack)%*%solve(V)%*%U_stack)%*%t(U_stack)%*%solve(V)%*%Z)

#settings for bootstrapping
chunklength <- sum(unlist(d)+1)
num_boot <- 1000

d_boot <- sample(1:chunklength, chunklength*num_boot, replace=TRUE)
d_boot_index <- split(d_boot, ceiling(seq_along(d_boot)/chunklength))
bet_boot_k <- rep(0,p+1)
Var_bet_boot <- array(0, c(p+1,p+1))
Z_boot <- unlist(beta_k)
Z_prime  <- sqrtm(solve(V))%*%Z_boot
U_prime  <- as.matrix(sqrtm(solve(V))%*%U_stack)

#initialize the mean and covariance matrix
mean_bet_boot <- rep(0, p+1)
for(bb in 1:num_boot){
  U_prime_boot <- U_prime[d_boot_index[[bb]],] 
  Z_prime_boot <- Z_prime[d_boot_index[[bb]],]
  
  bet_boot_k <- solve(t(U_prime[d_boot_index[[bb]],] )%*%U_prime[d_boot_index[[bb]],] )%*%t(U_prime[d_boot_index[[bb]],] )%*%Z_prime[d_boot_index[[bb]],]
  delta_boot <- bet_boot_k - mean_bet_boot
  mean_bet_boot <- mean_bet_boot + delta_boot / bb
  Var_bet_boot <- Var_bet_boot + (delta_boot %*% t(bet_boot_k - mean_bet_boot))
}
Var_bet_boot <- Var_bet_boot/(num_boot-1)

conf_int <- rbind(hat_bet[-1], round(hat_bet[-1] - qnorm(1-0.025/3)*sqrt(diag(Var_bet_boot))[-1],3), round( hat_bet[-1] +   qnorm(1-0.025/3)*sqrt(diag(Var_bet_boot))[-1],3))
colnames(conf_int)<- colnames(mydata)[-1];
rownames(conf_int) <- c("Estim","CI_LB","CI_UB")
  
#comparison to confidence interval with the entire data
#placeholders
reg_full <- lm(ggtp ~ map + bmi + smoking + drinking + blds + tc, data=mydata)
hat_bet_global <- reg_full$coef   
Var_bet_boot_global <- vcov(reg_full)
conf_int_global <- rbind( hat_bet_global[-1], 
                          round(hat_bet_global[-1] -   qnorm(1-0.025/3)*sqrt(diag(Var_bet_boot_global))[-1],3), 
                          round( hat_bet_global[-1] +   qnorm(1-0.025/3)*sqrt(diag(Var_bet_boot_global))[-1],3) )
colnames(conf_int_global)<- colnames(conf_int); rownames(conf_int_global) <- rownames(conf_int)

reg_x1 <- lm(ggtp ~ map)
hat_bet_naive <- reg_x1$coef
Var_bet_naive <- vcov(reg_x1)
conf_int_naive <- rbind(hat_bet_naive[-1], 
                        round(hat_bet_naive[-1] - qnorm(1-0.025/3)*sqrt(diag(Var_bet_naive))[-1],3), 
                  round(hat_bet_naive[-1] + qnorm(1-0.025/3)*sqrt(diag(Var_bet_naive))[-1],3))
  
conf_int_comparison <- cbind(conf_int[,1],conf_int_global[,1],conf_int_naive)
colnames(conf_int_comparison) <- c("meta", "full", "naive")
conf_int_comparison
  
#comparisons of the confidence intervals
ci_naive <- conf_int_comparison[-1,3]
ci_full <-   conf_int_comparison[-1,2]
ci_proposed <-  conf_int_comparison[-1,1]
  
#combine data into one data frame
data_combined <- data.frame(Group = factor(c("Proposed", "Naive", "Full"), levels = c("Proposed", "Naive", "Full")),
                            Mean = c(conf_int_comparison[1, 1], conf_int_comparison[1, 3], conf_int_comparison[1, 2]),
                            Lower_CI = c(ci_proposed[1], ci_naive[1], ci_full[1]),
                            Upper_CI = c(ci_proposed[2], ci_naive[2], ci_full[2]))

#making plots for comparision of confidence intervals
pdf("ci-comparison.pdf", width = 10, height = 3.5)
par(mfrow = c(1, 1))

data <- data.frame(Group = c("Naive", "Full", "Meta"),
                   Mean = c(conf_int_comparison[1, 3], conf_int_comparison[1, 2], conf_int_comparison[1, 1]),
                   Lower_CI = c(ci_naive[1], ci_full[1],ci_proposed[1]),
                   Upper_CI = c(ci_naive[2], ci_full[2],ci_proposed[2]))

plot(data$Mean,
     1:nrow(data),
     xlim = range(c(data$Lower_CI, data$Upper_CI)),
     ylim = c(0.5, nrow(data) + 0.5),
     xaxt = "n",
     yaxt = "n",
     xlab = expression(beta[map]), cex.lab=1.5,
     ylab = "",
     pch = 16,
     cex = 1.5,
     col = "black")

arrows(data$Lower_CI, 1:nrow(data), data$Upper_CI, 1:nrow(data),
       angle = 90, code = 3, length = 0.1, col = "black")

axis(2, at = 1:nrow(data), labels = data$Group, cex.axis = 1); axis(1, cex.axis = 1)
abline(h = 1:nrow(data), col = "lightgray", lty = "dotted")
dev.off()
