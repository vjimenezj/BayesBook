library(foreach)
library(doParallel)
library(iterators)
library(beepr)
library(rstan)
rstan_options(auto_write = TRUE)

# Local directory
dir_local <- Sys.getenv("BAYESBOOK_PATH")
if ((dir_local) == "") {
  dir_local <- "/data_lab_MAP/vjimenezj/HMDATR"
}
print(paste0("Local directory: ", dir_local))
source(paste0(dir_local, "/r_scripts/stan_utility.R"))
source(paste0(dir_local, "/r_scripts/plotting_inferences.R"))

# Names of the models to simulate and infer
model <- "/Poisson-lognormal_norm_lg"
            
# models <- c("/Poisson-lognormal_norm_lg", 
#             "/Poisson-lognormal", 
#             "/NB")
# title <- c("Poisson(beta,error)-Normalized-EffectiveLength", 
#             "Poisson(beta,error)", 
#             "NB(beta,error)")

# Other directories for reading, loading and writing files
dir_infer <- paste0(dir_local,"/inference_models")
dir_simu <- paste0(dir_local,"/simulation_models")

# Create results directory
dir.create(paste0(dir_local,"/results"))

dir_results <- paste0(dir_local,"/results", model, "_", Sys.Date())
i <- 1
while (file.exists(paste0(dir_results, "__", i))){
  i <- i + 1
}
dir_results <- paste0(dir_results, "__", i)
dir.create(dir_results)
dir.create(paste0(dir_results,"/prior_ensemble"))
dir.create(paste0(dir_results,"/inference_prior"))

# Setting aesthetic variables
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

# Bayesian ensemble simulation
G <- 100
S <- 6
design <- c(0, 0, 0, 1, 1, 1)
log_l_g <- log(rnorm(G, 1000, 10))
simu_data <- list("G" = G, "S" = S, "design" = design, "log_l_g" = log_l_g)
fit_ensemble <- stan(file=paste0(dir_simu, model, "_simulation.stan"),
            data=simu_data, iter=1000, warmup=0, chains=1,
            seed=4838282, algorithm="Fixed_param")
beep()

# Saving ensemble fit to load results
saveRDS(fit_ensemble, paste0(dir_results, "/prior_ensemble/fit_ensemble.rds"))
# # Load ensemble fit
# fit_ensemble <- readRDS(paste0(dir_results, "/prior_ensemble/fit_ensemble.rds"))

# Saving variables in the model to load them separately
simu_ensemble <- extract(fit_ensemble, permuted = FALSE)
expression <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE,
                      pars = c("expression"), include = TRUE)
#expression_real <- matrix(expression[1, 1, ], nrow = G, ncol = S, byrow = FALSE)
simu_expression <- vec_into_mat(expression)


simu_expression <- extract(fit_ensemble)$expression
R<-1000

# # Prior checks: Representing the simulated data as histograms
# 
# # Selecting the iteration with highest expressed gene in order to create
# # breaks for all histogram plots (range will always be big enough)
# max_dim <- which(log(simu_expression+1)==max(log(simu_expression+1)), arr.ind = TRUE)
# breaks_expr <- seq(from = 0, to = round(log(simu_expression[max_dim]+1), digits = 1) + 0.1, length.out = 61)
# B <- length(breaks_expr) - 1
# idx_expr <- rep(1:B, each = 2)
# x_expr <- rep(breaks_expr[2:(length(breaks_expr)-1)], each=2)
# x_expr <- append(x_expr, breaks_expr[length(breaks_expr)])
# x_expr <- append(breaks_expr[1], x_expr)
# 
# # hist(log(simu_expression[1, , ]+1), breaks=breaks_expr)$counts
# counts <- sapply(1:R, function(r) hist(log(simu_expression[r, , ]+1), breaks=breaks_expr, plot=FALSE)$counts)
# pad_counts <- sapply(1:R, function(r) do.call(cbind, lapply(idx_expr, function(n) counts[n, r])))
# 
# # # First simulation
# # c_superfine <- c("#8F272755")
# # plot(1, type="n", main="Prior Predictive Histogram Samples",
# #      xlim=c(x_expr[1], x_expr[length(x_expr)]), xlab="y", ylim=c(0, max(pad_counts)), ylab="")
# # lines(x_expr, pad_counts[,1], col=c_superfine, lw=2)
# # 
# # # First three simulations
# # plot(1, type="n", main="Prior Predictive Histogram Samples",
# #      xlim=c(x_expr[1], x_expr[length(x_expr)]), xlab="y", ylim=c(0, max(pad_counts)), ylab="")
# # lines(x_expr, pad_counts[,1], col=c_superfine, lw=2)
# # lines(x_expr, pad_counts[,2], col=c_superfine, lw=2)
# # lines(x_expr, pad_counts[,3], col=c_superfine, lw=2)
# # 
# # # First ten samples
# # plot(1, type="n", main="Prior Predictive Histogram Samples",
# #      xlim=c(x_expr[1], x_expr[length(x_expr)]), xlab="y", ylim=c(0, max(pad_counts)), ylab="")
# # for (r in 1:10)
# #   lines(x_expr, pad_counts[,r], col=c_superfine, lw=2)
# 
# # Plotting the quantiles of the histograms for all the simulations
# # (bad for multimodal distributions and correlations)
# probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# cred <- sapply(1:B, function(b) quantile(counts[b,], probs=probs))
# pad_cred <- do.call(cbind, lapply(idx_expr, function(n) cred[1:9, n]))
# 
# counts_obs <- hist(log(expression_real[c(1:G),]+1), breaks=breaks_expr, plot=FALSE)$counts
# pad_obs <- rep(counts_obs, each = 2)
# 
# pdf(paste0(dir_results, "/Prior_predictive_distribution.pdf"), width = 6, height = 6)
# plot(1, type="n", main="Prior Predictive Distribution",
#      xlim=c(x_expr[1], x_expr[length(x_expr)]), xlab="y", ylim=c(0, max(pad_obs)), ylab="")
# 
# polygon(c(x_expr, rev(x_expr)), c(pad_cred[1,], rev(pad_cred[9,])),
#         col = c_light, border = NA)
# polygon(c(x_expr, rev(x_expr)), c(pad_cred[2,], rev(pad_cred[8,])),
#         col = c_light_highlight, border = NA)
# polygon(c(x_expr, rev(x_expr)), c(pad_cred[3,], rev(pad_cred[7,])),
#         col = c_mid, border = NA)
# polygon(c(x_expr, rev(x_expr)), c(pad_cred[4,], rev(pad_cred[6,])),
#         col = c_mid_highlight, border = NA)
# lines(x_expr, pad_cred[5,], col=c_dark, lwd=2)
# lines(x_expr, pad_obs, col = "blue", lwd=2)
# dev.off()
# 
# # Estimated tail probability
# over10_simu <- length(simu_expression[log(simu_expression+1) > 10]) / length(simu_expression)
# over10_real <- length(expression_real[log(expression_real+1) > 10]) / length(expression_real)
# print(paste0("Percentage of genes with log expression over 10 in the simulation is ", 
#              over10_simu,
#              " whereas in the real dataset is ",
#              over10_real))



# Fitting simulated ensemble
print("Sampling ensemble performed. Starting trycatch")
tryCatch({
  registerDoParallel(makeCluster(detectCores()))
  
  simu_expression2 <- aperm(simu_expression, c(2, 3, 1))
  simu_expression_linear <- array(simu_expression2, dim = c(G * S, R))
  simu_expression_linear <- t(simu_expression_linear)
  
  simu_alpha2 <- aperm(simu_alpha, c(2, 3, 1))
  simu_alpha_linear <- array(simu_alpha2, dim = c(G * (C + 1), R))
  simu_alpha_linear <- t(simu_alpha_linear)
  
  simu_list <- data.matrix(cbind(simu_expression_linear, simu_alpha_linear, simu_phi, c(1:R)))
  
  # Compile the posterior fit model
  fit_model = stan_model(file=paste0(dir_model, "/", model, ".stan"))
  file=
  
  ensemble_output <- foreach(simu=iter(simu_list, by='row'),
                             .combine='cbind',
                             .packages = 'rstan',
                             .errorhandling = "remove",
                             .verbose = TRUE) %dopar% {
                               simu_expression_loop <- simu[1:(G*S)]
                               expression <- array(simu_expression_loop, dim = c(G, S))
                               simu_alpha_loop <- simu[(G * S + 1):(G * (S + C + 1))]
                               simu_phi_loop <- simu[(G * (S + C + 1)) + 1];
                               iter <- simu[(G * (S + C + 1)) + 2]
                               
                               # Fit the simulated observation
                               G <- dim(simu_expression)[2]
                               S <- dim(simu_expression)[3]
                               C <- dim(design)[1]
                               input_data <- list("G" = G, "S" = S, "C" = C, "expression" = expression, "design" = design)

                               a <- Sys.time()
                               print(paste0("Inititation of iteration ", iter, ": ", a))
                               capture.output(fit <- sampling(fit_model, data=input_data, seed=4938483))
                               b <- Sys.time()
                               print(paste0("End of iteration ", iter, ":", b))
                               print(paste0("Total time iteration ", iter, ": ", b-a))
                               
                               saveRDS(fit, paste0(dir_results, "/inference_prior/inference", iter, ".rds"))
                               
                               
                               # Compute diagnostics
                               source(paste0(dir_local, "/BayesWorkflow/stan_utility.R"))
                               
                               warning_code <- check_all_diagnostics(fit, quiet=TRUE)
                               
                               ##### Variable phi #####
                               # Compute rank of prior draw with respect to thinned posterior draws
                               sbc_rank_phi <- sum(simu_phi_loop < extract(fit)$phi[seq(1, 4000 - 8, 8)])
                               
                               # Compute posterior sensitivities
                               s_phi <- summary(fit, probs = c(), pars='phi')$summary
                               post_mean_phi <- s_phi[,1]
                               post_sd_phi <- s_phi[,3]
                               
                               prior_sd_phi <- sd(simu_phi)
                               
                               z_score_phi <- (post_mean_phi - simu_phi_loop) / post_sd_phi
                               contraction_phi <- 1 - (post_sd_phi / prior_sd_phi)**2
                               
                               ##### Variable alpha, first dimension (gene1 as example) #####
                               # Compute rank of prior draw with respect to thinned posterior draws
                               sbc_rank_alpha_1.1 <- sum(simu_alpha_loop[1] < extract(fit)$alpha[seq(1, 4000 - 8, 8), 1, 1])
                          
                               # Compute posterior sensitivities
                               s_alpha <- summary(fit, probs = c(), pars='alpha')$summary
                               post_mean_alpha_1.1 <- s_alpha[1, 1]
                               post_sd_alpha_1.1 <- s_alpha[1, 3]
                               
                               prior_sd_alpha_1.1 <- sd(simu_alpha[ , 1, 1])
                               
                               z_score_alpha_1.1 <- (post_mean_alpha_1.1 - simu_alpha_loop[1]) / post_sd_alpha_1.1
                               contraction_alpha_1.1 <- 1 - (post_sd_alpha_1.1 / prior_sd_alpha_1.1)**2
                               
                               ##### Variable alpha, second dimension (gene1 as example) #####
                               # Compute rank of prior draw with respect to thinned posterior draws
                               sbc_rank_alpha_1.2 <- sum(simu_alpha_loop[G + 1] < extract(fit)$alpha[seq(1, 4000 - 8, 8), 1, 2])
                               
                               # Compute posterior sensitivities
                               post_mean_alpha_1.2 <- s_alpha[4, 1]
                               post_sd_alpha_1.2 <- s_alpha[4, 3]
                               
                               prior_sd_alpha_1.2 <- sd(simu_alpha[ , 1, 2])
                               
                               z_score_alpha_1.2 <- (post_mean_alpha_1.2 - simu_alpha_loop[G + 1]) / post_sd_alpha_1.2
                               contraction_alpha_1.2 <- 1 - (post_sd_alpha_1.2 / prior_sd_alpha_1.2)**2
                               
                               ##### Variable alpha, third dimension (gene1 as example) #####
                               # Compute rank of prior draw with respect to thinned posterior draws
                               sbc_rank_alpha_1.3 <- sum(simu_alpha_loop[2 * G + 1] < extract(fit)$alpha[seq(1, 4000 - 8, 8), 1, 3])
                               
                               # Compute posterior sensitivities
                               post_mean_alpha_1.3 <- s_alpha[7, 1]
                               post_sd_alpha_1.3 <- s_alpha[7, 3]
                               
                               prior_sd_alpha_1.3 <- sd(simu_alpha[ , 1, 3])
                               
                               z_score_alpha_1.3 <- (post_mean_alpha_1.3 - simu_alpha_loop[2 * G + 1]) / post_sd_alpha_1.3
                               contraction_alpha_1.3 <- 1 - (post_sd_alpha_1.3 / prior_sd_alpha_1.3)**2
                               
                               # Saving interesting values
                               c(warning_code,
                                 sbc_rank_phi, z_score_phi, contraction_phi,
                                 sbc_rank_alpha_1.1, z_score_alpha_1.1, contraction_alpha_1.1,
                                 sbc_rank_alpha_1.2, z_score_alpha_1.2, contraction_alpha_1.2,
                                 sbc_rank_alpha_1.3, z_score_alpha_1.3, contraction_alpha_1.3)
                             }
}, finally={ stopImplicitCluster() })
b <- Sys.time()
b
b-a


saveRDS(ensemble_output, paste0(dir_results, "/inference_prior/ensemble_output.rds"))
#ensemble_output <- readRDS(paste0(dir, "/inference_prior/ensemble_output.rds"))

# Checking results for algorithmic calibration
str(ensemble_output)
class(ensemble_output)
head(ensemble_output)
head(ensemble_output[,1])

# 1:Checking diagnostics
warning_code <- ensemble_output[1,]
if (sum(warning_code) != 0) {
  cat("Some simulated posterior fits in the Bayesian ensemble encountered problems!\n")
  for (r in 1:R) {
    if (warning_code[r] != 0) {
      cat(sprintf('Replication %s of %s\n', r, R))
      parse_warning_code(warning_code[r])
      cat(sprintf('Simulated lambda = %s\n', simu_lambdas[r]))
      cat(" \n")
    }
  }
} else {
  cat("No posterior fits in the Bayesian ensemble encountered problems!\n")
}

# 2:Checking simulation-based calibration output
pdf(paste0(dir_results, "/SBC.pdf"), width = 6, height = 6)
par(mfrow=c(2,2))
##### Phi #####
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, plot=FALSE)
plot(sbc_hist, main="SBC for phi", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

##### Alpha: Dimension 1 #####
sbc_rank <- ensemble_output[5,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, plot=FALSE)
plot(sbc_hist, main="SBC for alpha (dimension 1)", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

##### Alpha: Dimension 2 #####
sbc_rank <- ensemble_output[8,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, plot=FALSE)
plot(sbc_hist, main="SBC for alpha (dimension 2)", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

##### Alpha: Dimension 3 #####
sbc_rank <- ensemble_output[11,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, plot=FALSE)
plot(sbc_hist, main="SBC for alpha (dimension 3)", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

dev.off()



# 3:Inferential calibration
pdf(paste0(dir_results, "/Inferential_calibration.pdf"), width = 6, height = 6)
par(mfrow=c(2,2))
##### Phi #####
z_score <- ensemble_output[3,]
contraction <- ensemble_output[4,]

plot(contraction, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Phi",
     xlim=c(0, 1), xlab="Posterior Contraction",
     ylim=c(-5, 5), ylab="Posterior z-Score")

##### Alpha (dimension 1) #####
z_score <- ensemble_output[6,]
contraction <- ensemble_output[7,]

plot(contraction, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Alpha (dimension 1)",
     xlim=c(0, 1), xlab="Posterior Contraction",
     ylim=c(-5, 5), ylab="Posterior z-Score")

##### Alpha (dimension 2) #####
z_score <- ensemble_output[9,]
contraction <- ensemble_output[10,]

plot(contraction, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Alpha (dimension 2)",
     xlim=c(0, 1), xlab="Posterior Contraction",
     ylim=c(-5, 5), ylab="Posterior z-Score")

##### Alpha (dimension 3) #####
z_score <- ensemble_output[12,]
contraction <- ensemble_output[13,]

plot(contraction, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Alpha (dimension 3)",
     xlim=c(0, 1), xlab="Posterior Contraction",
     ylim=c(-5, 5), ylab="Posterior z-Score")

dev.off()


# Output to print in log file
# Date
print(paste0("Execution time of the script: ", 
             Sys.time()))
# Computation time
print(paste0("Running time for the script: ", 
             proc.time()))

# Show currect R directory being executed for running the script
print(paste0("R directory being executed for running the script: ", 
             file.path(R.home("bin"), "R")))
# Arguments passed into Rscript
print(paste("Arguments passed to script:", args, sep = " "))
# Directories
print(paste0("Local directory: ", dir_local, 
             ", Model directory: ", dir_model,
             ", Data diretory: ", dir_data,
             ", Results directory: ", dir_results))

