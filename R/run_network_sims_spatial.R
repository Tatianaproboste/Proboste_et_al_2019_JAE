
## Pipeline script


run_network_sims_spatial = function(data.ecoli, data.social, 
                            n.simulations, n.cores, weighted.host = FALSE, 
                            spatial = FALSE, spatial_2014 = NULL, spatial_2015 = NULL ) {
  
  
  library(parallel)
  library(pbapply)
  
  # If n.cores is missing, set default number of cores for parallel processing
  if(missing(n.cores)){
    n.cores <- detectCores() - 1
  }
  
  #### If n.cores > 1, check parallel library loading ####
  if(n.cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n.cores)
    setDefaultCluster(cl)
    
    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(asnipe)), silent = TRUE)
    
    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {
      
      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/asnipe$", "", system.file(package = "asnipe"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))
      
      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(asnipe)), silent = TRUE)
      
      if(class(test_load2) == "try-error"){
        
        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(asnipe)), silent = TRUE)
        
        if(class(test_load3) == "try-error"){
          
          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)
          
        } else {
          parallel_compliant <- TRUE
        }
        
      } else {
        parallel_compliant <- TRUE
      }
      
    } else {
      parallel_compliant <- TRUE
    }
  } else {
    #If n.cores = 1, set parallel_compliant to FALSE
    parallel_compliant <- FALSE
    warning('Parallel loading not supported, calculations may crash!')
  }
  
  if(parallel_compliant){
    # If parallel loading is supported, export data, libraries and functions to each cluster
    #Export necessary data
    
      clusterExport (NULL, c('data.ecoli', 'data.social', 'spatial_2014',
                             'spatial_2015',
                             'n.simulations', 'weighted.host'),
                     envir = environment())
    
    
    #Export necessary functions
    clusterExport(NULL, c('get_ecoli_matrix', 'get_host_matrix',
                          'match_ecoli_matrices',
                          'extract_mod_coefs',
                          'compute_metrics','extract_MRQAP_coefs',
                          'extract_MRM_coefs'),
                  envir = environment())
    
    #Export necessary libraries
    clusterEvalQ(cl, library(asnipe))
    clusterEvalQ(cl, library(lme4))
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(ecodist))
    clusterEvalQ(cl, library(sna))
    
    
    #### After exporting data, functions and packages, run the simulations in parallel ####
    simulations <- pbapply::pblapply(seq_len(n.simulations), function(x){
      
      #### Create ecoli strain sharing matrices using a sampled cutoff ####
      # Randomly sample an ecoli similarity cutoff value
      ecoli.cutoff <- sample(c(0.90, 0.94), 1, FALSE)
      
      # Create strain sharing matrices for each year
      ecoli.2014 <- get_ecoli_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2014, 
                                     cutoff = ecoli.cutoff)
      
      ecoli.2015 <- get_ecoli_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2015, 
                                     cutoff = ecoli.cutoff)
      
      # Create matrices that only contain hosts sampled in both years
      # for ecoli community comparisons
      ecoli.annual.comms <- match_ecoli_matrices(matrix.2014 = ecoli.2014, 
                                                 matrix.2015 = ecoli.2015)
      
      
      #### Create host social association matrices using a sampled cutoff and 
      # association index ####
      # If using spatial, get the distance matrices
     
      
        #distance matrix
        dist.2014 <- as.matrix(spatial_2014)
        dist.2015 <- as.matrix(spatial_2015)
        
        
        # Randomly sample a host association threshold cutoff from a uniform 
        # distribution U[0, 0.15] 
        # Randomly sample an association index value
        host.cutoff <- sample(seq(0, 0.15, length.out = 100), 1, FALSE)
        host.index <- sample(c("HWI", "SRI"), 1, FALSE)
      
        
       
        
        
        host.2014 <- get_host_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2014,
                                     cutoff = host.cutoff,
                                     association.index = host.index)
        
        host.2015 <- get_host_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2015,
                                     cutoff = host.cutoff,
                                     association.index = host.index)
      
      
      host.annual.comms <- match_ecoli_matrices(matrix.2014 = host.2014, 
                                                matrix.2015 = host.2015) 
      
      #### Perform sanity checks to ensure that all matrices have matching names ####
      if(!all.equal(colnames(ecoli.2014), colnames(host.2014))){
        stop('Error in name matching of ecoli and host 2014 matrices')
      }
      
      #if(!all.equal(colnames(distance_matrix_2014), colnames(host.2014))){
      #  stop('Error in name matching of ecoli and host 2014 matrices')
      #}
      
      if(!all.equal(colnames(ecoli.2015), colnames(host.2015))){
        stop('Error in name matching of ecoli and host 2015 matrices')
      }
      
      #### Compute network metrics for each year ####
      metrics.2014 <- compute_metrics(matrix.ecoli = ecoli.2014, 
                                      matrix.social = host.2014,
                                      year = 2014, weighted.host = weighted.host)
      
      metrics.2015 <- compute_metrics(matrix.ecoli = ecoli.2015, 
                                      matrix.social = host.2015,
                                      year = 2015, weighted.host = weighted.host)
      all.metrics <- rbind(metrics.2014, metrics.2015)
      rm(metrics.2014, metrics.2015)
      
      # Merge in the sex information of hosts
      data.social$sex <- ifelse(data.social$ReproStatus == "Male", "Male", "Female")
      names.sex <- unique(data.social[c("Name", "sex")])
      
      all.metrics <- merge(all.metrics, names.sex, by = 'Name')
      all.metrics$sex.fac <- as.factor(all.metrics$sex)
      all.metrics$year.fac <- as.factor(all.metrics$Year)
      rm(names.sex)
      
      #### Run linear mixed models and extract coefficients ####
      # Run degree centrality model
      degree.mod <- lme4::lmer(SOC.degree ~ scale(ECOLI.degree)*year.fac + 
                                 sex.fac + (1|Name), data = all.metrics, REML = F)
      
      # Extract model results
      degree.results <- extract_mod_coefs(model = degree.mod, 
                                          parameters <- c('Ecoli.degree','Year.2015',
                                                          'Sex.Male',
                                                          'Ecoli.degree*Year.2015'))
      rm(degree.mod)
      
      # Run betwenness model and extract results
      btw.mod <- lme4::lmer(SOC.betweenness ~ scale(ECOLI.betweenness)*year.fac + 
                              sex.fac + (1|Name), data = all.metrics, REML = F)
      
      betweeness.results <- extract_mod_coefs(model = btw.mod, 
                                              parameters <- c('Ecoli.betweenness','Year.2015',
                                                              'Sex.Male',
                                                              'Ecoli.betweenness*Year.2015'))
      rm(btw.mod)
      
      # Run eigenvector centrality model and extract results
      egv.mod <- lme4::lmer(SOC.evcent ~ scale(ECOLI.evcent)*year.fac + 
                              sex.fac + (1|Name), data = all.metrics, REML = F)
      
      eigencent.results <- extract_mod_coefs(model = egv.mod, 
                                             parameters <- c('Ecoli.evcent','Year.2015',
                                                             'Sex.Male',
                                                             'Ecoli.evcent*Year.2015'))
      rm(egv.mod)
      
      #### Compare network 'community' structures using MRQAP and MRM ####
      # Each model is run and then the 95% confidence interval of coefficients
      # is evaluated to test for 'significance'
      netlm.2014 <- sna::netlm(host.2014, ecoli.2014,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.2014.results <- extract_MRQAP_coefs(netlm.2014, host.2014)
      rm(netlm.2014)
      
      mrm.2014 <- ecodist::MRM(dist(host.2014) ~ dist(ecoli.2014), 
                               nperm = 100)
      MRM.2014.results <- extract_MRM_coefs(mrm.2014)
      rm(mrm.2014)
      
      netlm.2015 <- sna::netlm(host.2015, ecoli.2015,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.2015.results <- extract_MRQAP_coefs(netlm.2015, host.2015)
      
      mrm.2015 <- ecodist::MRM(dist(host.2015) ~ dist(ecoli.2015), 
                               nperm = 100)
      MRM.2015.results <- extract_MRM_coefs(mrm.2015)
      rm(mrm.2015)
      
      # Spatial / Ecoli
      netlm.2014.dist <- sna::netlm(dist.2014, ecoli.2014,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 1000)
      
      MRQAP.2014.results.dist <- extract_MRQAP_coefs(netlm.2014.dist, dist.2014)
      rm(netlm.2014.dist)
      
      
      netlm.2015.dist <- sna::netlm(dist.2015, ecoli.2015,
                                    mode = "graph",
                                    nullhyp = "qapspp",
                                    reps = 1000)
      
      MRQAP.2015.results.dist <- extract_MRQAP_coefs(netlm.2015.dist, dist.2015)
      rm(netlm.2015.dist)
      
      
      mrm.2014.dist <- ecodist::MRM(dist(dist.2014) ~ dist(ecoli.2014), 
                               nperm = 100)
      MRM.2014.results.dist <- extract_MRM_coefs(mrm.2014.dist)
      rm(mrm.2014.dist)
      
      mrm.2015.dist <- ecodist::MRM(dist(dist.2015) ~ dist(ecoli.2015), 
                               nperm = 100)
      MRM.2015.results.dist <- extract_MRM_coefs(mrm.2015.dist)
      rm(mrm.2015)
      
     
      ###
      
      netlm.ecoli <- sna::netlogit(ecoli.annual.comms$matrix.2014.sub,
                                   ecoli.annual.comms$matrix.2015.sub,
                                   mode = 'graph',
                                   reps = 100)
      MRQAP.ecoli.results <- extract_MRQAP_coefs(netlm.ecoli, ecoli.annual.comms$matrix.2014.sub)
      
      mrm.ecoli <- ecodist::MRM(dist(ecoli.annual.comms$matrix.2014.sub) ~ 
                                  dist(ecoli.annual.comms$matrix.2015.sub), 
                                nperm = 100)
      MRM.ecoli.results <- extract_MRM_coefs(mrm.ecoli)
      rm(mrm.ecoli)
      
      
      
      #host / host
      netlm.host <- sna::netlm(host.annual.comms$matrix.2014.sub,
                               host.annual.comms$matrix.2015.sub,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.host.results <- extract_MRQAP_coefs(netlm.host, host.annual.comms$matrix.2014.sub)
      
      mrm.host <- ecodist::MRM(dist(host.annual.comms$matrix.2014.sub) ~ 
                                 dist(host.annual.comms$matrix.2015.sub), 
                               nperm = 100)
      MRM.host.results <- extract_MRM_coefs(mrm.host)
      
      
      return(list(ecoli.cutoff = ecoli.cutoff,
                  host.cutoff = host.cutoff,
                  host.index = host.index,
                  degree.results = degree.results,
                  betweeness.results = betweeness.results,
                  eigencent.results = eigencent.results,
                  MRQAP.2014.results = MRQAP.2014.results,
                  MRM.2014.results = MRM.2014.results,
                  MRQAP.2015.results = MRQAP.2015.results,
                  MRM.2015.results = MRM.2015.results,
                  MRQAP.ecoli.results = MRQAP.ecoli.results,
                  MRM.ecoli.results = MRM.ecoli.results,
                  MRQAP.host.results = MRQAP.host.results,
                  MRM.host.results = MRM.host.results,
                  MRQAP.2014.results.dist = MRQAP.2014.results.dist,
                  MRQAP.2015.results.dist = MRQAP.2015.results.dist,
                  MRM.2014.results.dist = MRM.2014.results.dist,
                  MRM.2015.results.dist = MRM.2015.results.dist 
      ))
    }, cl = cl)
    stopCluster(cl)
    
  } else {
    # If parallel loading not supported, use lapply instead (will take a long time!)
    simulations <- lapply(seq_len(n.simulations), function(x){
      
      #### Report progress of the simulation exercise
      cat('Processing iteration', x, 'of', n.simulations, '...\n')
      
      #### Create ecoli strain sharing matrices using a sampled cutoff ####
      # Randomly sample an ecoli similarity cutoff value
      ecoli.cutoff <- sample(c(0.90, 0.94), 1, FALSE)
      
      # Create strain sharing matrices for each year
      ecoli.2014 <- get_ecoli_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2014, 
                                     cutoff = ecoli.cutoff)
      
      ecoli.2015 <- get_ecoli_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2015, 
                                     cutoff = ecoli.cutoff)
      
      # Create matrices that only contain hosts sampled in both years
      # for ecoli community comparisons
      ecoli.annual.comms <- match_ecoli_matrices(matrix.2014 = ecoli.2014, 
                                                 matrix.2015 = ecoli.2015)
      
      
      
      #distance matrix
      dist.2014 <- as.matrix(spatial_2014)
      dist.2015 <- as.matrix(spatial_2015)
      
      
      #### Create host social association matrices using a sampled cutoff and 
      # association index ####
      host.cutoff <- sample(seq(0, 0.15, length.out = 100), 1, FALSE) 
      host.index <- sample(c("HWI", "SRI"), 1, FALSE)
      
      
        # Randomly sample a host association threshold cutoff from a uniform 
        # distribution U[0, 0.15] 
        # Randomly sample an association index value
        
        
        host.2014 <- get_host_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2014,
                                     cutoff = host.cutoff,
                                     association.index = host.index)
        
        host.2015 <- get_host_matrix(data.ecoli = data.ecoli, 
                                     data.social = data.social, 
                                     year = 2015,
                                     cutoff = host.cutoff,
                                     association.index = host.index)
      
      
      host.annual.comms <- match_ecoli_matrices(matrix.2014 = host.2014, 
                                                matrix.2015 = host.2015) 
      
      
      #### Perform sanity checks to ensure that all matrices have matching names ####
      if(!all.equal(colnames(ecoli.2014), colnames(host.2014))){
        stop('Error in name matching of ecoli and host 2014 matrices')
      }
      
      if(!all.equal(colnames(ecoli.2015), colnames(host.2015))){
        stop('Error in name matching of ecoli and host 2015 matrices')
      }
      
      #### Compute network metrics for each year ####
      metrics.2014 <- compute_metrics(matrix.ecoli = ecoli.2014, 
                                      matrix.social = host.2014,
                                      year = 2014, weighted.host = weighted.host)
      
      metrics.2015 <- compute_metrics(matrix.ecoli = ecoli.2015, 
                                      matrix.social = host.2015,
                                      year = 2015, weighted.host = weighted.host)
      all.metrics <- rbind(metrics.2014, metrics.2015)
      rm(metrics.2014, metrics.2015)
      
      # Merge in the sex information of hosts
      data.social$sex <- ifelse(data.social$ReproStatus == "Male", "Male", "Female")
      names.sex <- unique(data.social[c("Name", "sex")])
      
      all.metrics <- merge(all.metrics, names.sex, by = 'Name')
      all.metrics$sex.fac <- as.factor(all.metrics$sex)
      all.metrics$year.fac <- as.factor(all.metrics$Year)
      rm(names.sex)
      
      #### Run linear mixed models and extract coefficients ####
      # Run degree centrality model
      degree.mod <- lme4::lmer(SOC.degree ~ scale(ECOLI.degree)*year.fac + 
                                 sex.fac + (1|Name), 
                               data = all.metrics, REML = F)
      
      # Extract model results
      degree.results <- extract_mod_coefs(model = degree.mod, 
                                          parameters <- c('Ecoli.degree','Year.2015',
                                                          'Sex.Male',
                                                          'Ecoli.degree*Year.2015'))
      rm(degree.mod)
      
      # Run betwenness model and extract results
      btw.mod <- lme4::lmer(SOC.betweenness ~ scale(ECOLI.betweenness)*year.fac + 
                              sex.fac + (1|Name), data = all.metrics, REML = F)
      
      betweeness.results <- extract_mod_coefs(model = btw.mod, 
                                              parameters <- c('Ecoli.betweenness','Year.2015',
                                                              'Sex.Male',
                                                              'Ecoli.betweenness*Year.2015'))
      rm(btw.mod)
      
      # Run eigenvector centrality model and extract results
      egv.mod <- lme4::lmer(SOC.evcent ~ scale(ECOLI.evcent)*year.fac + 
                              sex.fac + (1|Name), data = all.metrics, REML = F)
      
      eigencent.results <- extract_mod_coefs(model = egv.mod, 
                                             parameters <- c('Ecoli.evcent','Year.2015',
                                                             'Sex.Male',
                                                             'Ecoli.evcent*Year.2015'))
      rm(egv.mod)
      
      #### Compare network 'community' structures using MRQAP and MRM ####
      # Each model is run and then the 95% confidence interval of coefficients
      # is evaluated to test for 'significance'
      netlm.2014 <- sna::netlm(host.2014, ecoli.2014,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.2014.results <- extract_MRQAP_coefs(netlm.2014, host.2014)
      rm(netlm.2014)
      
      mrm.2014 <- ecodist::MRM(dist(host.2014) ~ dist(ecoli.2014), 
                               nperm = 100)
      MRM.2014.results <- extract_MRM_coefs(mrm.2014)
      rm(mrm.2014)
      
      netlm.2015 <- sna::netlm(host.2015, ecoli.2015,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.2015.results <- extract_MRQAP_coefs(netlm.2015, host.2015)
      
      mrm.2015 <- ecodist::MRM(dist(host.2015) ~ dist(ecoli.2015), 
                               nperm = 100)
      MRM.2015.results <- extract_MRM_coefs(mrm.2015)
      rm(mrm.2015)
      
      
      netlm.ecoli <- sna::netlogit(ecoli.annual.comms$matrix.2014.sub,
                                   ecoli.annual.comms$matrix.2015.sub,
                                   mode = 'graph',
                                   reps = 100)
      MRQAP.ecoli.results <- extract_MRQAP_coefs(netlm.ecoli, ecoli.annual.comms$matrix.2014.sub)
      
      mrm.ecoli <- ecodist::MRM(dist(ecoli.annual.comms$matrix.2014.sub) ~ dist(ecoli.annual.comms$matrix.2015.sub), 
                                nperm = 100)
      MRM.ecoli.results <- extract_MRM_coefs(mrm.ecoli)
      rm(mrm.ecoli)
      
      
      # spatial / Ecoli
      #MRQAP
      #2014
      netlm.2014.dist <- sna::netlm(dist.2014, ecoli.2014,
                                    mode = "graph",
                                    nullhyp = "qapspp",
                                    reps = 1000)
      
      MRQAP.2014.results.dist <- extract_MRQAP_coefs(netlm.2014.dist, dist.2014)
      rm(netlm.2014.dist)
      
      #2015
      netlm.2015.dist <- sna::netlm(dist.2015, ecoli.2015,
                                    mode = "graph",
                                    nullhyp = "qapspp",
                                    reps = 1000)
      
      MRQAP.2015.results.dist <- extract_MRQAP_coefs(netlm.2015.dist, dist.2015)
      rm(netlm.2015.dist)
      
      #MRM
      mrm.2014.dist <- ecodist::MRM(dist(dist.2014) ~ dist(ecoli.2014), 
                                    nperm = 100)
      MRM.2014.results.dist <- extract_MRM_coefs(mrm.2014.dist)
      rm(mrm.2014.dist)
      
      mrm.2015.dist <- ecodist::MRM(dist(dist.2015) ~ dist(ecoli.2015), 
                                    nperm = 100)
      MRM.2015.results.dist <- extract_MRM_coefs(mrm.2015.dist)
      rm(mrm.2015)
      
     
      
      
      #host / hots
      netlm.host <- sna::netlm(host.annual.comms$matrix.2014.sub,
                               host.annual.comms$matrix.2015.sub,
                               mode = "graph",
                               nullhyp = "qapspp",
                               reps = 100)
      MRQAP.host.results <- extract_MRQAP_coefs(netlm.host, host.annual.comms$matrix.2014.sub)
      
      mrm.host <- ecodist::MRM(dist(host.annual.comms$matrix.2014.sub) ~ 
                                 dist(host.annual.comms$matrix.2015.sub), 
                               nperm = 100)
      MRM.host.results <- extract_MRM_coefs(mrm.host)
      
      return(list(ecoli.cutoff = ecoli.cutoff,
                  host.cutoff = host.cutoff,
                  host.index = host.index,
                  degree.results = degree.results,
                  betweeness.results = betweeness.results,
                  eigencent.results = eigencent.results,
                  MRQAP.2014.results = MRQAP.2014.results,
                  MRM.2014.results = MRM.2014.results,
                  MRQAP.2015.results = MRQAP.2015.results,
                  MRM.2015.results = MRM.2015.results,
                  MRQAP.ecoli.results = MRQAP.ecoli.results,
                  MRM.ecoli.results = MRM.ecoli.results,
                  MRQAP.host.results = MRQAP.host.results,
                  MRM.host.results =  MRM.host.results,
                  MRQAP.2014.results.dist = MRQAP.2014.results.dist,
                  MRQAP.2015.results.dist = MRQAP.2015.results.dist,
                  MRM.2014.results.dist = MRM.2014.results.dist,
                  MRM.2015.results.dist = MRM.2015.results.dist 
      ))
    })
    
  }
  
  #### Processing simulation results ####
  library(dplyr)
  
  #1.Degree
  degree.coef.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'degree.results'),
                                               'Coefficients'))
  degree.coef.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'degree.results'),
                                                  'Coefficients')) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Lower.coef.95 = quantile(Estimate, probs = 0.025),
                     Median.coef = quantile(Estimate, probs = 0.5),
                     Upper.coef = quantile(Estimate, probs = 0.975), 
                     Lower.coef.85 = quantile(Estimate, probs = 0.075),
                     Median.coef.85 = quantile(Estimate, probs = 0.5),
                     Upper.coef.85 = quantile(Estimate, probs = 0.925))
  
  degree.intercept.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'degree.results'),
                                                    'Host.averages'))
  degree.intercept.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'degree.results'),
                                                       'Host.averages')) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(Lower.intercept = quantile(Intercept, probs = 0.025),
                     Median.intercept = quantile(Intercept, probs = 0.5),
                     Upper.intercept = quantile(Intercept, probs = 0.975))
  
  #2. Betweeness
  betweeness.coef.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'betweeness.results'),
                                                   'Coefficients'))
  betweeness.coef.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'betweeness.results'),
                                                      'Coefficients')) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Lower.coef = quantile(Estimate, probs = 0.025),
                     Median.coef = quantile(Estimate, probs = 0.5),
                     Upper.coef = quantile(Estimate, probs = 0.975),
                     Lower.coef.85 = quantile(Estimate, probs = 0.075),
                     Median.coef.85 = quantile(Estimate, probs = 0.5),
                     Upper.coef.85 = quantile(Estimate, probs = 0.925))
  
  betweeness.intercept.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'betweeness.results'),
                                                        'Host.averages'))
  betweeness.intercept.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'betweeness.results'),
                                                           'Host.averages')) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(Lower.intercept = quantile(Intercept, probs = 0.025),
                     Median.intercept = quantile(Intercept, probs = 0.5),
                     Upper.intercept = quantile(Intercept, probs = 0.975))
  #3. eigencent
  eigencent.coef.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'eigencent.results'),
                                                  'Coefficients'))
  eigencent.coef.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'eigencent.results'),
                                                     'Coefficients')) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Lower.coef = quantile(Estimate, probs = 0.025),
                     Median.coef = quantile(Estimate, probs = 0.5),
                     Upper.coef = quantile(Estimate, probs = 0.975), 
                     Lower.coef.85 = quantile(Estimate, probs = 0.075),
                     Median.coef.85 = quantile(Estimate, probs = 0.5),
                     Upper.coef.85 = quantile(Estimate, probs = 0.925))
  
  eigencent.intercept.raw <- do.call(rbind, purrr::map(purrr::map(simulations, 'eigencent.results'),
                                                       'Host.averages'))
  eigencent.intercept.summary = do.call(rbind, purrr::map(purrr::map(simulations, 'eigencent.results'),
                                                          'Host.averages')) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(Lower.intercept = quantile(Intercept, probs = 0.025),
                     Median.intercept = quantile(Intercept, probs = 0.5),
                     Upper.intercept = quantile(Intercept, probs = 0.975))
  #4. MRQAP 2014
  MRQAP.2014.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2014.results')))
  MRQAP.2014.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2014.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  #.4.1. MRQAP 2014 (ecoli & distance)
  MRQAP.2014.dist.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2014.results.dist')))
  MRQAP.2014.dist.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2014.results.dist'))) %>%
  dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
               R = as.numeric(as.character(R))) %>%
  dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                 Median.coef = quantile(Coefficient, probs = 0.5),
                 Upper.coef = quantile(Coefficient, probs = 0.975),
                 Lower.r= quantile(R, probs = 0.025),
                 Median.r = quantile(R, probs = 0.5),
                 Upper.r = quantile(R, probs = 0.975),
                 Prop.significant = length(Significance[which(Significance == 'significant')])/
                 length(Significance))
  
  
  #.4.2. MRQAP 2015 (ecoli & distance)
  MRQAP.2015.dist.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2015.results.dist')))
  MRQAP.2015.dist.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2015.results.dist'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  
  #4.3. MRM 2014 (ecoli & distance)
  MRM.2014.dist.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2014.results.dist')))
  MRM.2014.dist.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2014.results.dist'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  #4.4 MRM 2015 (ecoli & distance)
  MRM.2015.dist.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2015.results.dist')))
  MRM.2015.dist.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2015.results.dist'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  
  #5. MRM 2014
  MRM.2014.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2014.results')))
  MRM.2014.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2014.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  #6. MRQAP 2015
  MRQAP.2015.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2015.results')))
  MRQAP.2015.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.2015.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  #7. MRM 2015
  MRM.2015.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2015.results')))
  MRM.2015.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.2015.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  #8. MRQAP Ecoli
  MRQAP.ecoli.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.ecoli.results')))
  MRQAP.ecoli.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.ecoli.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  #9 MRM Ecoli
  MRM.ecoli.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.ecoli.results')))
  MRM.ecoli.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.ecoli.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  #10 MRQAP host
  MRQAP.host.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.host.results')))
  MRQAP.host.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRQAP.host.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  #11 MRM host 
  MRM.host.raw <- data.frame(do.call(rbind, purrr::map(simulations, 'MRM.host.results')))
  MRM.host.summary = data.frame(do.call(rbind, purrr::map(simulations, 'MRM.host.results'))) %>%
    dplyr::mutate(Coefficient = as.numeric(as.character(Coefficient)),
                  R = as.numeric(as.character(R))) %>%
    dplyr::summarise(Lower.coef = quantile(Coefficient, probs = 0.025),
                     Median.coef = quantile(Coefficient, probs = 0.5),
                     Upper.coef = quantile(Coefficient, probs = 0.975),
                     Lower.r= quantile(R, probs = 0.025),
                     Median.r = quantile(R, probs = 0.5),
                     Upper.r = quantile(R, probs = 0.975),
                     Prop.significant = length(Significance[which(Significance == 'significant')])/
                       length(Significance))
  
  #### Extract statistics on types of networks that were used ####
  ecoli.cutoff.raw <- unlist(purrr::map(simulations, 'ecoli.cutoff'))
  ecoli.cutoff.summary <- matrix(NA, 1, 2)
  ecoli.cutoff.summary[1, 1] <- length(ecoli.cutoff.raw[which(ecoli.cutoff.raw == 0.90)])
  ecoli.cutoff.summary[1, 2] <- length(ecoli.cutoff.raw[which(ecoli.cutoff.raw == 0.94)])
  colnames(ecoli.cutoff.summary) <- c('Number.sims.at.0.90', 'Number.sims.at.0.94')
  
  host.index.raw <- unlist(purrr::map(simulations, 'host.index'))
  host.index.summary <- matrix(NA, 1, 2)
  host.index.summary[1, 1] <- length(host.index.raw[which(host.index.raw == 'HWI')])
  host.index.summary[1, 2] <- length(host.index.raw[which(host.index.raw == 'SRI')])
  colnames(host.index.summary) <- c('Number.sims.using.HWI', 'Number.sims.using.SRI')
  
  host.cutoff.raw <- unlist(purrr::map(simulations, 'host.cutoff'))
  host.cutoff.summary <- quantile(host.cutoff.raw, probs = c(0.025, 0.5, 0.975))
  
  return(list(degree.coef.summary = degree.coef.summary,
              degree.intercept.summary = degree.intercept.summary,
              betweeness.coef.summary = betweeness.coef.summary,
              betweeness.intercept.summary = betweeness.intercept.summary,
              eigencent.coef.summary = eigencent.coef.summary,
              eigencent.intercept.summary = eigencent.intercept.summary,
              MRQAP.2014.summary = MRQAP.2014.summary,
              MRQAP.2014.dist.summary = MRQAP.2014.dist.summary,
              MRQAP.2015.dist.summary = MRQAP.2015.dist.summary,
              MRM.2014.dist.summary = MRM.2014.dist.summary,
              MRM.2015.dist.summary = MRM.2015.dist.summary,
              MRM.2014.summary = MRM.2014.summary,
              MRQAP.2015.summary = MRQAP.2015.summary,
              MRM.2015.summary = MRM.2015.summary,
              MRQAP.ecoli.summary = MRQAP.ecoli.summary,
              MRM.ecoli.summary = MRM.ecoli.summary,
              ecoli.cutoff.summary = ecoli.cutoff.summary,
              host.index.summary = host.index.summary,
              host.cutoff.summary = host.cutoff.summary,
              degree.coef.raw = degree.coef.raw,
              degree.intercept.raw = degree.intercept.raw,
              betweeness.coef.raw = betweeness.coef.raw, 
              betweeness.intercept.raw = betweeness.intercept.raw,
              eigencent.coef.raw = eigencent.coef.raw,
              eigencent.intercept.raw = eigencent.intercept.raw,
              MRQAP.2014.raw = MRQAP.2014.raw,
              MRQAP.2014.dist.raw = MRQAP.2014.dist.raw,
              MRQAP.2015.dist.raw = MRQAP.2015.dist.raw,
              MRM.2014.dist.raw = MRM.2014.dist.raw,
              MRM.2015.dist.raw = MRM.2015.dist.raw,
              MRM.2014.raw = MRM.2014.raw,
              MRQAP.2015.raw = MRQAP.2015.raw, 
              MRM.2015.raw = MRM.2015.raw, 
              MRQAP.ecoli.raw = MRQAP.ecoli.raw,
              MRM.ecoli.raw = MRM.ecoli.raw,
              ecoli.cutoff.raw = ecoli.cutoff.raw,
              host.index.raw = host.index.raw,
              host.cutoff.raw = host.cutoff.raw,
              MRQAP.host.summary = MRQAP.host.summary, 
              MRM.host.summary = MRM.host.summary ))
}



