packages <- c("latex2exp", "rstudioapi")

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
load("gendata-ci-ma-fe.Rdata")

rho_0.3 <-  results_50000[[1]][[1]][[1]];   
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

###1. create figure 1
###1-1. create pdf for subfigure (a) 
pdf_filename_panel1 <- paste0("panel1.pdf")
pdf(pdf_filename_panel1, width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )
true_value <- bet_estim["True Value", ]
estimates <- bet_estim[-1, ]

plot(NULL, xlim=c(1, num_beta), xlab="", ylab="", xaxt='n', cex=1.4, cex.lab=1.4, ylim=c(-5.5,4))
axis(1, at = 1:num_beta, cex.axis=1.4, labels = beta_expressions)
points(1:num_beta, as.numeric(true_value),  pch=3, col=colors[1])

# Plot estimates for each rho
for (i in 1:nrow(estimates)) {
  points(1:num_beta, as.numeric(estimates[i, ]),  pch=pch_values[i], col=colors[2])
}
legend("bottomright", legend=c("True",TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), 
                               TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')), pch=c(3,pch_values), col=c("red",rep("black",4)))
dev.off()

###1-2. create pdf for subfigure (b) (Variance and MSE)
pdf_filename_panel2 <- paste0("panel2.pdf")
pdf(pdf_filename_panel2 , width = 5, height = 5)
par(mar = c(4, 3, 3, 1.5) )
variance <- bet_var
mse <- bet_mse

plot(NULL, xlim=c(1, num_beta), ylim=c(0,max(mse[1:4,])*1.5), xlab="", ylab="", xaxt='n', cex=1.4,  cex.lab=1.4)
axis(1, at = 1:num_beta, cex.axis=1.4, labels = beta_expressions)

# Plot variance and MSE for each rho
for (i in 1:4){
  points(1:num_beta, as.numeric(variance[i, ]), col=colors[1], pch=pch_values_filled[i], cex=0.8)
  points(1:num_beta, as.numeric(mse[i, ]), col=colors[2], pch=pch_values[i], cex=1.6)
}

legend("topleft", 
       legend=c(TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')), 
       pch=pch_values)
dev.off()


###1-3 Create pdf for subfigure (c)
pdf_filename_panel3 <- paste0("panel3.pdf")
pdf(pdf_filename_panel3 , width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )

coverage_real <- bet_cov[1:4,]
coverage_oracle <- bet_cov[5:8,]

plot(NULL, xlim=c(1, num_beta), ylim=c(0.92, 0.98), xlab="", ylab="", main="", xaxt='n', cex=1.4,  cex.lab=1.4)
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


hat_bet_rho0.3 <- rho_0.3$Estimates
hat_bet_rho0.6 <- rho_0.6$Estimates
hat_bet_rho0.8 <- rho_0.8$Estimates
hat_bet_rho0.9 <- rho_0.9$Estimates

bet<- true_bet
###1-4. create pdf for subfigure (d)
pdf_filename_panel4 <- paste0("panel4.pdf")
pdf(pdf_filename_panel4 , width = 5, height = 5)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1) + 0.1, oma = c(2.5, 1, 1.5, 0))

plot(density(as.numeric(hat_bet_rho0.3[,2])), xlab="", main=TeX('$\\rho_X=\\,0.3$'), cex.main=1.4)
abline(v=bet[1], col="red")
plot(density(as.numeric(hat_bet_rho0.6[,2])), xlab="", main=TeX('$\\rho_X=\\,0.6$'), cex.main=1.4)
abline(v=bet[1], col="red")
plot(density(as.numeric(hat_bet_rho0.8[,2])), xlab="", main=TeX('$\\rho_X=\\,0.8$'), cex.main=1.4)
abline(v=bet[1], col="red")
plot(density(as.numeric(hat_bet_rho0.9[,2])), xlab="", main=TeX('$\\rho_X=\\,0.9$'), cex.main=1.4)
abline(v=bet[1], col="red")

dev.off()