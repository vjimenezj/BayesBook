library(rstan)
library(ggplot2)

options(scipen = 999)

# Names of the models to simulate and infer
models <- c("/Poisson-lognormal_norm_lg", 
            "/Poisson-lognormal", 
            "/NB")
titles <- c("Poisson(beta,error)-Normalized-EffectiveLength", 
            "Poisson(beta,error)", 
            "NB(beta,error)")


# Local directory
dir_local <- Sys.getenv("BAYESBOOK_PATH")
if ((dir_local) == "") {
  dir_local <- "/data_lab_MAP/vjimenezj/BayesBook"
}
print(paste0("Local directory: ", dir_local))

# Model directory
dir_models_sim <- paste0(dir_local, "/simulation_models")
dir_models_infer <- paste0(dir_local, "/inference_models")

# Loading plotting inferences functions
source(paste0(dir_local, "/r_scripts/plotting_inferences.R"))

# Simulating and inferring for each of the models
for (j in c(1:length(models))) {
  # Defining the data for sampling from the model
  G <- 100
  R <- 10
  # design <- rbind(c(0, 0, 0, 0, 1, 1, 1, 1), c(0, 0, 1, 1, 0, 0, 1, 1))
  design <- c(0, 0, 0, 0, 1, 1, 1, 1)
  # S <- dim(design)[2]
  # C <- dim(design)[1]
  S <- 8
  C <- 2
  log_l_g <- log(rnorm(G, 1000, 10))
  # Data list for sample generation
  simu_data <- list("G" = G, "S" = S, "design" = design, "log_l_g" = log_l_g)
  # Sampling with Stan
  fit_ensemble <- stan(file=paste0(dir_models_sim, models[j], ".stan"),
                       data=simu_data, iter=1, warmup=0, chains=1,
                       refresh=500, seed=4838282, algorithm="Fixed_param")
  
  # Extracting samples and true parameter values to infer with Stan and for initialization
  expression <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE,
                             pars = c("expression"), include = TRUE)
  #expression_real <- matrix(expression[1, 1, ], nrow = G, ncol = S, byrow = FALSE)
  expression_real <- vec_into_mat(expression)

  for (k in c(1:length(model2inf))) {
    input_data <- list("G" = G, "S" = S, "expression" = expression_real, "design" = design)
    options(mc.cores = 4)
    mc.cores = parallel::detectCores()
    nchains <- 4
    iter_per_chain <- 2000
    fit <- stan(file=paste0(dir_models_infer, models[k], "_inference.stan"),
                 data=input_data, seed=493848, iter = 2 * iter_per_chain, 
                 refresh = 400, chains = nchains, 
                 control = list(max_treedepth = 10))
    
    inference_summary <- summary(fit, probs = c(0.025, 0.975))$summary
    simu_ensemble <- extract(fit_ensemble, permuted = FALSE, inc_warmup = FALSE)
    iteration <- 1
    chain <- 1
    
    par_names_means <- c("mu_theta", "sigma_theta", "theta")
    par_names_changes <- "beta"
    par_names_errors <- "error"
    par_names <- list(par_names_means, par_names_changes, par_names_errors)
    
    
    for (parameters in par_names) {
      title <- paste0("Simulation_", titles[j], "-Inference_", titles[k],
                      "-Samples_", S, "-", paste(parameters, collapse = ","))
      
      p <- plotting_credibility(parameters, simu_ensemble = simu_ensemble,
                                iteration = 1, chain = 1, title = title,
                                inference_summary = inference_summary)
      print(p)
      ggsave(paste0(dir_local, "/Exploration/Plots/", title, ".png"))
    }
    
    if (substr(title_inf[j], start = 1, stop = 2) == "NB") {
      phi_infer <- extract(fit)$phi
      phi_title <- paste0("Simulation_", title_gen[k], "-Inference_", title_inf[j], "-Samples_", S, "-Phi=1")
      print(hist(phi_infer, main = phi_title))
      print(abline(v = 1))
      ggsave(paste0(dir_local, "/Exploration/Plots/", phi_title, ".png"))
    }
  }
}

