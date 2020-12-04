library(rstan)
rstan_options(auto_write = TRUE)


a <- Sys.time()
# Local directory
dir_local <- Sys.getenv("HMDATR_PATH")
if ((dir_local) == "") {
  dir_local <- "/data3/vjimenezj/HMDATR"
}
print(paste0("Local directory: ", dir_local))

# Arguments passed into Rscript
args = commandArgs(trailingOnly=TRUE)
model <- "NB-fixedCov"
design <- rbind(c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 1, 1, 0, 0, 1, 1))
dir_data <- paste0(dir_local,"/DEGbayesData/expression_real.rds")
Iterations <- 2000




# Other directories for reading, loading and writing files
dir_model <- paste0(dir_local,"/Models/", model)
util <- new.env()
source(paste0(dir_local, "/BayesWorkflow/stan_utility.R"))


# Create results directory
dir.create(paste0(dir_local,"/FittedModels"))
dir.create(paste0(dir_local,"/FittedModels/", model))

dir_results <- paste0(dir_local,"/FittedModels/", model, 
                      "/Estimation_", model, "_", Sys.Date() )
i <- 1
while (file.exists(paste0(dir_results, "__", i))){
  i <- i + 1
}
dir_results <- paste0(dir_results, "__", i)
dir.create(dir_results)
sessionInfo()

# Loading real data
expression_real <- readRDS(dir_data)


# 4:Fitting the real model
G <- dim(expression_real)[1]
S <- dim(expression_real)[2]
design <- rbind(c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 1, 1, 0, 0, 1, 1))
C <- dim(design)[1]
input_data <- list("G" = G, "S" = S, "C" = C, 
                   "expression" = expression_real,
                   "design" = design)

#writeLines(readLines(paste0(dir_model, "/", model, "_post_pred.stan")))
options(mc.cores = 4)
fit <- stan(file=paste0(dir_model, "/", model, "_post_pred.stan"),
            data=input_data, seed=4938483, iter = Iterations, 
            chains = 4, control = list(max_treedepth = 9))
saveRDS(fit, paste0(dir_results, "/fit.rds"))

# Diagnose posterior fits
check_all_diagnostics(fit)

params = extract(fit)

pdf(paste0(dir_results, "/Histogram_phi.pdf"), width = 5, height = 6)
hist(params$phi, main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
dev.off()

pdf(paste0(dir_results, "/Histogram_alpha_dim1.pdf"), width = 5, height = 6)
hist(params$alpha[ , 1, 1], main="", xlab="alpha (dimension 1)", yaxt='n', 
     ylab="", col=c_dark, border=c_dark_highlight)
dev.off()

pdf(paste0(dir_results, "/Histogram_alpha_dim2.pdf"), width = 5, height = 6)
hist(params$alpha[ , 1, 2], main="", xlab="alpha (dimension 2)", yaxt='n', 
     ylab="", col=c_dark, border=c_dark_highlight)
dev.off()

pdf(paste0(dir_results, "/Histogram_alpha_dim3.pdf"), width = 5, height = 6)
hist(params$alpha[ , 1, 3], main="", xlab="alpha (dimension 3)", yaxt='n', 
     ylab="", col=c_dark, border=c_dark_highlight)
dev.off()


##### Posterior retrodictive checks
# Representing the simulated data as histograms

# Selecting the iteration with highest expressed gene in order to create
# breaks for all histogram plots (range will always be big enough)
max_dim <- which(log(params$expression_pred+1)==max(log(params$expression_pred+1)),
                 arr.ind = TRUE)[1]

#hist(log(params$expression_pred[max_dim, , ]+1), breaks=60)

breaks_expr <- hist(log(params$expression_pred[max_dim, , ]+1),
                    breaks=60, plot=FALSE)$breaks

B <- length(breaks_expr) - 1
idx_expr <- rep(1:B, each = 2)
x_expr <- rep(breaks_expr[2:(length(breaks_expr)-1)], each=2)
x_expr <- append(x_expr, breaks_expr[length(breaks_expr)])
x_expr <- append(breaks_expr[1], x_expr)

# hist(log(simu_expression[1, , ]+1), breaks=breaks_expr)$counts
counts <- sapply(1:(4*Iterations), 
                 function(r) hist(log(params$expression_pred[r, , ]+1),
                                  breaks=breaks_expr, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B, function(b) quantile(counts[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx_expr, function(n) cred[1:9, n]))

counts_obs <- hist(log(expression_real[c(1:G),]+1), 
                   breaks=breaks_expr, plot=FALSE)$counts
pad_obs <- rep(counts_obs, each = 2)

pdf(paste0(dir_results, "/Posterior_retrodictive_checks.pdf"), width = 6, height = 6)
plot(1, type="n", main="Posterior Retrodictive Distribution",
     xlim=c(-0.5, B + 0.5), xlab="y",
     ylim=c(0, max(c(counts_obs, cred[9,]))), ylab="")

polygon(c(x_expr, rev(x_expr)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x_expr, rev(x_expr)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x_expr, rev(x_expr)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x_expr, rev(x_expr)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x_expr, pad_cred[5,], col=c_dark, lwd=2)

#lines(x_expr, pad_obs, col="white", lty=1, lw=2.5)
lines(x_expr, pad_obs, col="black", lty=1, lw=2)
dev.off()

b <- Sys.time()
print(paste0("Total execution time: ", b - a))






