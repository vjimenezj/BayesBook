library(rstan)
library(ggplot2)
library(beepr)

options(scipen = 999)

# Names of the models to simulate and infer
# model <- "/Poisson-lognormal_lg"
# model <- "/Poisson-lognormal_norm_lg"
# model <- "/Poisson-lognormal"
model <- "/NB"
# model_name <- "Poisson(beta,error)-Normalized-EffectiveLength"
# model_name <- "Poisson(beta,error)-EffectiveLength"
# model_name <- "Poisson(beta,error)"
model_name <- "NB(beta,error)"

# Local directory
dir_local <- Sys.getenv("BAYESBOOK_PATH")
if ((dir_local) == "") {
  dir_local <- "/data3/vjimenezj/BayesBook"
}
print(paste0("Local directory: ", dir_local))

# Model directory
dir_models_sim <- paste0(dir_local, "/simulation_models")
dir_models_infer <- paste0(dir_local, "/inference_models")

# Loading plotting inferences functions
source(paste0(dir_local, "/r_scripts/plotting_inferences.R"))

# Simulating and inferring for each of the models
# Defining the data for sampling from the model
G <- 100
R <- 10
# design <- rbind(c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 1, 1, 0, 0, 1, 1))
design <- c(0, 0, 0, 1, 1, 1)
# S <- dim(design)[2]
# C <- dim(design)[1]
S <- 6
C <- 2
log_l_g <- log(rnorm(G, 1000, 10))
# Data list for sample generation
simu_data <- list("G" = G, "S" = S, "design" = design, "log_l_g" = log_l_g)
# Sampling with Stan
fit_ensemble <- stan(file=paste0(dir_models_sim, model, "_simulation.stan"),
                     data=simu_data, iter=1, warmup=0, chains=1,
                     refresh=500, seed=4838282, algorithm="Fixed_param")

# Extracting samples and true parameter values to infer with Stan and for initialization
expression <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE,
                           pars = c("expression"), include = TRUE)
#expression_real <- matrix(expression[1, 1, ], nrow = G, ncol = S, byrow = FALSE)
expression_real <- vec_into_mat(expression)
input_data <- list("G" = G, "S" = S, "expression" = expression_real, 
                   "design" = design, "log_l_g" = log_l_g)
options(mc.cores = 4)
mc.cores = parallel::detectCores()
nchains <- 4
iter_per_chain <- 2000
fit <- stan(file=paste0(dir_models_infer, model, "_inference.stan"),
             data=input_data, seed=493848, iter = 2 * iter_per_chain, 
             refresh = 400, chains = nchains, 
             control = list(max_treedepth = 10))
beep()


inference_summary <- summary(fit, probs = c(0.025, 0.975))$summary
simu_ensemble <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE)
iteration <- 1
chain <- 1

par_names_means <- c("mu_alpha", "sigma_alpha", "alpha")
par_names_changes <- "beta"
par_names_errors <- "error"
#par_names_norm <- "log_norm_factors"
par_names <- list(par_names_means, par_names_changes, 
                  par_names_errors, par_names_norm)

dir.create(paste0(dir_local, "/results/"))
for (parameters in par_names) {
  title <- paste0("Simulation_", model_name, "-Inference_", model_name,
                  "-Samples_", S, "-", paste(parameters, collapse = ","))
  
  p <- plotting_credibility(parameters, simu_ensemble = simu_ensemble,
                            iteration = 1, chain = 1, title = title,
                            inference_summary = inference_summary)
  print(p)
  ggsave(paste0(dir_local, "/results/", title, ".png"))
}

if (substr(model_name, start = 1, stop = 2) == "NB") {
  phi_infer <- as.data.frame(extract(fit)$phi)
  colnames(phi_infer) <- "Phi"
  phi_title <- paste0("Simulation_", model_name, "-Inference_", model_name, "-Samples_", S, "-Phi=1")
  
  p <- ggplot(phi_infer, aes(Phi)) + 
   # geom_histogram(aes(y=..density..), binwidth = 0.01) +
    geom_density(alpha=.5, fill="#00FF00") +
    geom_vline(xintercept = 1, color = "blue", size=1.5) +
    theme_classic()
  
  p
  ggsave(paste0(dir_local, "/results/", phi_title, ".png"))
}
