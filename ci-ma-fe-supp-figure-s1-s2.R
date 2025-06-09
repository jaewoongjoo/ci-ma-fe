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

rho_0   <-  results_50000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]];
rho_0.1 <-  results_50000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]];   
rho_0.2 <-  results_50000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]];   
rho_0.3 <-  results_50000[[1]][[1]][[1]][[1]][[1]][[1]][[2]];
rho_0.4 <-  results_50000[[1]][[1]][[1]][[1]][[1]][[2]];
rho_0.5 <-  results_50000[[1]][[1]][[1]][[1]][[2]];
rho_0.6 <-  results_50000[[1]][[1]][[1]][[2]];
rho_0.7 <-  results_50000[[1]][[1]][[2]];
rho_0.8 <-  results_50000[[1]][[2]];
rho_0.9 <-  results_50000[[2]];  
  
rho_0_m10000     <-  results_10000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]];   
rho_0.1_m10000   <-  results_10000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]];   
rho_0.2_m10000   <-  results_10000[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]];
rho_0.3_m10000   <-  results_10000[[1]][[1]][[1]][[1]][[1]][[1]][[2]];
rho_0.4_m10000   <-  results_10000[[1]][[1]][[1]][[1]][[1]][[2]];
rho_0.5_m10000   <-  results_10000[[1]][[1]][[1]][[1]][[2]];
rho_0.6_m10000   <-  results_10000[[1]][[1]][[1]][[2]];
rho_0.7_m10000   <-  results_10000[[1]][[1]][[2]];
rho_0.8_m10000   <-  results_10000[[1]][[2]];
rho_0.9_m10000   <-  results_10000[[2]];  

pch_values_experiment <- 2
pch_values_experiment_filled <- c(15,16)

###1. create figure S-1

coverage_0.3 <- rbind(rho_0.3_m10000$Real[6,], rho_0.3$Real[6,], rho_0.3$Oracle[6,])
coverage_0.6 <- rbind(rho_0.6_m10000$Real[6,], rho_0.6$Real[6,], rho_0.6$Oracle[6,])
coverage_0.7 <- rbind(rho_0.7_m10000$Real[6,], rho_0.7$Real[6,], rho_0.7$Oracle[6,])
coverage_0.8 <- rbind(rho_0.8_m10000$Real[6,], rho_0.8$Real[6,], rho_0.8$Oracle[6,])
coverage_0.9 <- rbind(rho_0.9_m10000$Real[6,], rho_0.9$Real[6,], rho_0.9$Oracle[6,])

rownames(coverage_0.3) <- c("rho0_Real_m10000", "rho0_Real_m50000", "rho0_Oracle")
rownames(coverage_0.6) = rownames(coverage_0.3); rownames(coverage_0.7) = rownames(coverage_0.3); 
rownames(coverage_0.8) = rownames(coverage_0.3); rownames(coverage_0.9) = rownames(coverage_0.3)

colnames(coverage_0.3) <- c("bet1","bet2","bet3","bet4","bet5","bet6","bet7","bet8")
colnames(coverage_0.6) = colnames(coverage_0.3); colnames(coverage_0.7) = colnames(coverage_0.3); 
colnames(coverage_0.8) = colnames(coverage_0.3); colnames(coverage_0.9) = colnames(coverage_0.3)

pdf("panel5.pdf",width = 10, height = 10)
par(mfrow=c(2,3))
par(mar = c(4, 3, 3, 1) )

###rho=0.3
plot(NULL, xlim=c(1, 8), ylim=c(0.89, 0.97), xlab="", ylab="", main=TeX('$\\rho_X=\\,0.3$'), xaxt='n', cex=1.4,  cex.lab=1.4, cex.main=1.8)
axis(2, at = seq(from=0.89, to=0.97, by= 0.01) )
axis(1, at = 1:num_beta, cex.axis=1.6, 
     labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")
abline(h=0.93, lty=3, col="blue")
for (i in 1:2){
  points(1:num_beta, as.numeric(coverage_0.3[i, ]), col=colors[1],   pch=pch_values_experiment_filled[i], cex=1.5)
}
points(1:num_beta, as.numeric(coverage_0.3[3, ]), col=colors[2],   pch=pch_values_experiment, cex=1.5)

###rho=0.6
plot(NULL, xlim=c(1, 8), ylim=c(0.89, 0.97), xlab="", ylab="", main=TeX('$\\rho_X=\\,0.6$'), xaxt='n', cex=1.4,  cex.lab=1.4, cex.main=1.8)
axis(2, at = seq(from=0.89, to=0.97, by= 0.01) )
axis(1, at = 1:num_beta, cex.axis=1.6, 
     labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")
abline(h=0.93, lty=3, col="blue")
for (i in 1:2){
  points(1:num_beta, as.numeric(coverage_0.6[i, ]), col=colors[1],   pch=pch_values_experiment_filled[i], cex=1.5)
}
points(1:num_beta, as.numeric(coverage_0.6[3, ]), col=colors[2],   pch=pch_values_experiment, cex=1.5)

###rho=0.7
plot(NULL, xlim=c(1, 8), ylim=c(0.89, 0.97), xlab="", ylab="", main=TeX('$\\rho_X=\\,0.7$'), xaxt='n', cex=1.4,  cex.lab=1.4, cex.main=1.8)
axis(2, at = seq(from=0.89, to=0.97, by= 0.01) )
axis(1, at = 1:num_beta, cex.axis=1.6, 
     labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")
abline(h=0.93, lty=3, col="blue")
for (i in 1:2){
  points(1:num_beta, as.numeric(coverage_0.7[i, ]), col=colors[1],   pch=pch_values_experiment_filled[i], cex=1.5)
}
points(1:num_beta, as.numeric(coverage_0.7[3, ]), col=colors[2],   pch=pch_values_experiment, cex=1.5)

###rho=0.8
plot(NULL, xlim=c(1, 8), ylim=c(0.89, 0.97), xlab="", ylab="", main=TeX('$\\rho_X=\\,0.8$'), xaxt='n', cex=1.4,  cex.lab=1.4, cex.main=1.8)
axis(2, at = seq(from=0.89, to=0.97, by= 0.01) )
axis(1, at = 1:num_beta, cex.axis=1.6, 
     labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")
abline(h=0.93, lty=3, col="blue")
for (i in 1:2){
  points(1:num_beta, as.numeric(coverage_0.8[i, ]), col=colors[1],   pch=pch_values_experiment_filled[i], cex=1.5)
}
points(1:num_beta, as.numeric(coverage_0.8[3, ]), col=colors[2],   pch=pch_values_experiment, cex=1.5)

###rho=0.9
plot(NULL, xlim=c(1, 8), ylim=c(0.89, 0.97), xlab="", ylab="", main=TeX('$\\rho_X=\\,0.9$'), xaxt='n', cex=1.4,  cex.lab=1.4, cex.main=1.8)
axis(2, at = seq(from=0.89, to=0.97, by= 0.01) )
axis(1, at = 1:num_beta, cex.axis=1.6, 
     labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")
abline(h=0.93, lty=3, col="blue")
for (i in 1:2){
  points(1:num_beta, as.numeric(coverage_0.9[i, ]), col=colors[1],   pch=pch_values_experiment_filled[i], cex=1.5)
}
points(1:num_beta, as.numeric(coverage_0.9[3, ]), col=colors[2],   pch=pch_values_experiment, cex=1.5)

###legend
plot.new()
legend("bottomleft", inset=c(0.01, 0.3),
       legend=c("10,000", "50,000", TeX('$\\infty$\\,(oracle)')),
       pch=c(15, 16, 2), 
       bty="n", 
       col=c("red", "red", "black"),
       cex=2,  # Adjust this value to control the size of the legend
       title="Size of external data set",
       adj = c(0, 0.5),  # Aligns text to the left
       xpd=TRUE)
dev.off()

###2. create figure S-2
bet_all_estim <- rbind(rho_0.3$Real[1,], rho_0.3$Real[2,], rho_0.6$Real[2,], rho_0.8$Real[2,], rho_0.9$Real[2,])

bet_all_var <- rbind(rho_0.3$Real[4,], rho_0.6$Real[4,], rho_0.8$Real[4,], rho_0.9$Real[4,],
                     rho_0.3$Rootn[4,], rho_0.6$Rootn[4,], rho_0.8$Rootn[4,], rho_0.9$Rootn[4,],
                     rho_0.3$Oracle[4,], rho_0.6$Oracle[4,], rho_0.8$Oracle[4,], rho_0.9$Oracle[4,])

bet_all_mse <- rbind(rho_0.3$Real[5,], rho_0.6$Real[5,], rho_0.8$Real[5,], rho_0.9$Real[5,],
                     rho_0.3$Rootn[5,], rho_0.6$Rootn[5,], rho_0.8$Rootn[5,], rho_0.9$Rootn[5,],
                     rho_0.3$Oracle[5,], rho_0.6$Oracle[5,], rho_0.8$Oracle[5,], rho_0.9$Oracle[5,])

bet_all_cov <- rbind(rho_0.3$Real[6,], rho_0.6$Real[6,], rho_0.8$Real[6,], rho_0.9$Real[6,],
                     rho_0.3$Rootn[6,], rho_0.6$Rootn[6,], rho_0.8$Rootn[6,], rho_0.9$Rootn[6,],
                     rho_0.3$Oracle[6,], rho_0.6$Oracle[6,], rho_0.8$Oracle[6,], rho_0.9$Oracle[6,])

rownames(bet_all_estim) <- c("True Value",  "rho0.3",  "rho0.6", "rho0.8", "rho0.9")
colnames(bet_all_estim) <- c("beta_1", "beta_2", "beta_3", "beta_4", "beta_5", "beta_6", "beta_7", "beta_8")
colnames(bet_all_cov) <- colnames(bet_all_estim); colnames(bet_all_var) <- colnames(bet_all_estim); colnames(bet_all_mse)<- colnames(bet_all_estim);
rownames(bet_all_var) <- c("rho0.3_Real", "rho0.6_Real",  "rho0.8_Real", "rho0.9_Real",
                           "rho0.3_Rootn", "rho0.6_Rootn", "rho0.8_Rootn", "rho0.9_Rootn",
                           "rho0.3_Oracle", "rho0.6_Oracle", "rho0.8_Oracle", "rho0.9_Oracle") ; 
rownames(bet_all_cov) <- rownames(bet_all_var); 
rownames(bet_all_mse) <- rownames(bet_all_var); 

colors <- c("red","black")
pch_values <- c(0,1,2,5)
pch_values_filled <- c(15, 16, 17, 18) 

pdf("panel3_rootn.pdf",width = 5, height = 5)
par(mar = c(4, 3, 3, 1) )
# Extract data for the plot
coverage_real <- bet_all_cov[1:4,]
coverage_rootn <- bet_all_cov[5:8,]
# Plot setup
plot(NULL, xlim=c(1, 8), ylim=c(0.92, 0.98), xlab="", ylab="", main="", xaxt='n', cex=1.4,  cex.lab=1.4)
axis(1, at = 1:num_beta, cex.axis=1.4, labels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]),
                                                  expression(beta[5]),expression(beta[6]),expression(beta[7]),expression(beta[8])))
abline(h=0.95, lty=2, col="red")

for (i in 1:4){
  points(1:num_beta, as.numeric(coverage_real[i, ]), col=colors[1],   pch=pch_values_filled[i], cex=1)
  points(1:num_beta, as.numeric(coverage_rootn[i, ]), col='gray50', pch=pch_values[i], cex=1.6)
}

legend("topleft", 
       legend=c(TeX('$\\rho_X=\\,0.3$'), TeX('$\\rho_X=\\,0.6$'), TeX('$\\rho_X=\\,0.8$'), TeX('$\\rho_X=\\,0.9$')),
       pch=pch_values)
dev.off()

