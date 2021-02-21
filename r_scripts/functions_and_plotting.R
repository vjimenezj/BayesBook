# Function to select fit1_summary rows based on parameter row names
select_summary <- function(fit_summary, par_names) {
  parameters <- gsub("\\[.*?\\]", "", rownames(fit_summary))
  summary_pars <- fit_summary[parameters %in% par_names, ]
  return(summary_pars)
}

# Select real_values iteration used for inference and parameters of interest using dimnames
select_simu <- function(simu, iteration, chain, par_names) {
  parameters <- gsub("\\[.*?\\]", "", dimnames(simu)[[3]])
  simu_pars <- simu[iteration, chain, parameters %in% par_names]
  return(simu_pars)
}

# Function to shape list into matrix of parameters
vec_into_mat <- function(vec) {
  indexes <- gsub("[^0-9,]", "",  dimnames(vec)[[3]])
  indexes <- strsplit(indexes, ",")
  indexes  <-  as.data.frame(matrix(as.numeric(unlist(indexes)), ncol=length(unlist(indexes[1])), byrow = TRUE))
  mat <- matrix(vec, nrow = max(indexes[,1]), ncol = max(indexes[,2]), byrow = FALSE)
  return(mat)
} 

# Function to compare credibility intervals for selected parameters with real values of those parameters
plotting_credibility <- function (par_names, simu_ensemble, iteration, chain,inference_summary, title) {
  inference_summary <- inference_summary[, c("mean", "2.5%", "97.5%")]
  colnames(inference_summary) <- c("estimation_model", "conf_lower_model", "conf_upper_model")
  inf_sum_pars <- select_summary(inference_summary, par_names)
  simu_pars <- select_simu(simu_ensemble, iteration, chain, par_names)
  colors <- gsub("\\[.*?\\]", "", rownames(inf_sum_pars))
  inf_sum_pars <- inf_sum_pars[match(names(simu_pars), rownames(inf_sum_pars)),]
  if (sum(rownames(inf_sum_pars) != names(simu_pars)) == 0) {
    data <- as.data.frame(cbind(inf_sum_pars, simu_pars, colors))
  } else {print("Error in plot")}
  p <- ggplot(data, aes(x = as.numeric(simu_pars), y = as.numeric(estimation_model))) + 
    geom_point(size = 1, aes(color = colors)) +
    geom_errorbar(aes(ymax = as.numeric(conf_upper_model), ymin = as.numeric(conf_lower_model))) +
    geom_abline(slope = 1, intercept = 0) + 
    xlab("Real values") + 
    ylab("Estimated values") # + 
    # ggtitle(title)
  return(p)
}

print("Plotting functions loaded")
