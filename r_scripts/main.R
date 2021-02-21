library(rstan)
library(ggplot2)

options(scipen = 999)

# Local directory (If the .Renviron file has not been created in the local repository with the location of
# the repository, this local variable has to be changed into the directory of the repository)
dir_local <- '/data_lab_MAP/vjimenezj/BayesBook/'
# dir_local <- Sys.getenv("BAYESBOOK_PATH")
# dir_local <- "local_directory_of_BookBayes_repository"

# Loading plotting inferences functions
source(paste0(dir_local, "/r_scripts/functions_and_plotting.R"))

# Simulating and inferring for each of the models
# Defining the data for sampling from the model prior
G <- 100
R <- 10
design <- c(0, 0, 0, 1, 1, 1)
S <- 6
C <- 2
eff_length <- rnorm(G, 1000, 10)
phi <- 1
# Data list for sample generation
simu_data <- list("G" = G, "S" = S, "design" = design, "eff_length" = eff_length, "phi" = phi)
# Sampling with Stan
fit_ensemble <- stan(file=paste0(dir_local, "/simulation_models/NB_simulation.stan"),
                     data=simu_data, iter=1, warmup=0, chains=1,
                     refresh=500, seed=4838282, algorithm="Fixed_param")

# Extracting samples and true parameter values to infer with Stan and for initialization
expression <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE,
                           pars = c("expression"), include = TRUE)
expression_real <- vec_into_mat(expression)
input_data <- list("G" = G, "S" = S, "expression" = expression_real, 
                   "design" = design, "eff_length" = eff_length)
options(mc.cores = 4)
mc.cores = parallel::detectCores()
nchains <- 4
iter_per_chain <- 2000
fit <- stan(file=paste0(dir_local, "/inference_models/NB_inference.stan"),
             data=input_data, seed=493848, iter = 2 * iter_per_chain, 
             refresh = 400, chains = nchains, 
             control = list(max_treedepth = 10))

# figures for a comparison among the inferred parameters and its simulated values
inference_summary <- summary(fit, probs = c(0.025, 0.975))$summary
simu_ensemble <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE)
iteration <- 1
chain <- 1
par_names_means <- c("mu_alpha", "sigma_alpha", "alpha")
par_names_changes <- "beta"
par_names_norm <- "log_norm_factors"
par_names <- list(par_names_means, par_names_changes, par_names_norm)

dir.create(paste0(dir_local, "/figures/"))
i <- 1
for (parameters in par_names) {
  title <- paste0("Simulation and Inference from NB model - ", paste(parameters, collapse = ","))
  p <- plotting_credibility(parameters, simu_ensemble = simu_ensemble,
                            iteration = 1, chain = 1, title = title,
                            inference_summary = inference_summary)
  assign(paste0('p',i), p)
  print(p)
  ggsave(paste0(dir_local, "/figures/", title, ".tiff"), dpi = 200)
  i <- i+1
}

# Histogram for phi posterior distribution
phi_infer <- as.data.frame(extract(fit)$phi)
colnames(phi_infer) <- "Phi"
phi_title <- paste0("Simulation_and_Inference_from_NB_model-Phi=1")
  
p <- ggplot(phi_infer, aes(Phi)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.01) +
 # geom_density(alpha=.5, fill="#00FF00") +
  geom_vline(xintercept = 1, color = "blue", size=1.5) +
  theme_classic()
p
ggsave(paste0(dir_local, "/figures/", phi_title, ".tiff"), dpi = 200)

library(ggpubr)
ggarrange(p1, p2, p3, p, 
          labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2)
ggsave(paste0(dir_local, "/figures/Figure3.tiff"), dpi = 200)

