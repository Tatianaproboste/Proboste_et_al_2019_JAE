

compute_metrics = function(matrix.ecoli, matrix.social, year, weighted.host = FALSE){
  
  #### Make sure year and association.index specifications are allowed ####
  if(!(year %in% c(2014, 2015))){
    # if year isn't legal, issue an error
    stop('Year must either be a numeric value matching either 2014 or 2015')
  }
  
  ## Compute degree centrality for each network
  # gmode "graph" indicates that edges are undirected  
  s.strength <- sna::degree(matrix.social, ignore.eval = weighted.host, gmode = "graph") 
  e.strength <- sna::degree(matrix.ecoli,  gmode = "graph", ignore.eval = T)

  
  s.names <- colnames(matrix.social)
  e.names <- colnames(matrix.ecoli)
  
  ## Compute flow betweenness
  s.betweenness <- sna::flowbet(matrix.social, gmode = 'graph', cmode="normflow", 
                                rescale = T)
  e.betweenness <- sna::flowbet(matrix.ecoli, gmode = 'graph',
                                rescale = T)
  
  ## Compute eigen centrality
  s.evcent <- sna::evcent(matrix.social, gmode = 'graph')
  e.evcent <- sna::evcent(matrix.ecoli, gmode = 'graph')
  
  t <- data.frame(Name = s.names, SOC.degree = s.strength, 
                 SOC.betweenness = s.betweenness,
                 SOC.evcent = s.evcent, 
                 ECOLI.degree = e.strength, 
                 ECOLI.betweenness = e.betweenness, 
                 ECOLI.evcent = e.evcent,
                 Year = year)
  
  return(t)
  
}