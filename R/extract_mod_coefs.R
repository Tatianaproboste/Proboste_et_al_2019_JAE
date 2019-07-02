# Extract host intercepts and coefficients from lmer models
extract_mod_coefs = function(model, parameters){
  library(dplyr)
  
  host.averages <- data.frame(coef(model)$Name)
  host.averages$Name <- rownames(host.averages)
  host.averages %>%
    dplyr::mutate(Intercept = X.Intercept.) %>%
    dplyr::select(Name, Intercept) -> host.averages
  
  # Extract regression coefficients
  mod.coefs <- data.frame(summary(model)$coefficients)[-1,] 
  mod.coefs$Parameter <- parameters
  mod.coefs %>%
    dplyr::select(Parameter, Estimate) -> mod.coefs
  rownames(mod.coefs) <- NULL
  
  return(list(Host.averages = host.averages,
              Coefficients = mod.coefs))
}

# Extract coefficients and significance levels from MRQAP models
extract_MRQAP_coefs = function(model, y){
  
  output <- matrix(NA, 1, 3)
  output[1, 1] <- round(model$coefficients[2], 4)
  output[1, 2] <- sum((model$fitted.values - mean(y)) ^ 2) / sum((y - mean(y)) ^ 2)
    output[1, 3] <- if(model$pgreqabs[2] > 0.05){
    'not.significant'
  } else{
    'significant'
  }
  colnames(output) <- c('Coefficient', 'R', 'Significance')
  
  return(output)
}

# Extract coefficients and significance levels from MRM
extract_MRM_coefs = function(model){

  output <- matrix(NA, 1, 3)
  output[1, 1] <- round(model$coef[2, 1], 4)
  output[1, 2] <- model$r.squared[1]
  output[1, 3] <- if(model$r.squared[2] > 0.05) {
    'not.significant'
  } else{
    'significant'
  }
  colnames(output) <- c('Coefficient', 'R', 'Significance')
  
  return(output)
}
