get_host_matrix = function(data.ecoli, data.social, year, cutoff, association.index){
  
  #### Make sure year and association.index specifications are allowed ####
  if(!(year %in% c(2014, 2015))){
    # if year isn't legal, issue an error
    stop('Year must either be a numeric value matching either 2014 or 2015')
  }
  
  if(!(association.index %in% c("HWI", "SRI"))){
    stop('association.index must be a character string matching either "HWI" or "SRI"')
  }
  
  if(!is.numeric(cutoff)){
    stop('cutoff must be a positive numeric value between 0 and 0.15')
  }
  
  #### Subset the host and ecoli datasets for the appropriate year ####
  if(year == 2014){
    d.soc <- subset(data.social, Year == 2014)
    d.coli <- subset(data.ecoli, year == 2014)
    
  } else{
    d.soc <- subset(data.social, Year == 2015)
    d.coli <- subset(data.ecoli, year == 2015)
  }
  
  #### Clean the data bases with only the individuals that appear in both datasets ####
  d.soc <- d.soc[d.soc$Name %in% d.coli$Name, ]
  d.coli <- d.coli[d.coli$Name %in% d.soc$Name, ]
  
  #### Get the sex of each individual ####
  d.soc$sex <- ifelse(d.soc$ReproStatus == "Male", "Male", "Female")
  names_sex <- unique(d.soc[c("Name", "sex")])
  
  # Put names in one column and uniquegroup in another column
  d.soc$uniquegroup <- paste(d.soc$SurveyNumber, d.soc$GroupNumber)
  survey_modified <- data.frame(d.soc$Name, d.soc$uniquegroup)
  names(survey_modified) <- c("Name", "Group")
  survey_modified <- survey_modified[order(survey_modified$Name), ]
  
  #### Creating the host social assication matrix ####
  # Get group by individual format
  s.gbi <- asnipe::get_group_by_individual(survey_modified,data_format="individuals")
  
  # Create an association matrix from the data based on the specified association.index value
  s.hwi <- asnipe::get_network(s.gbi, data_format = "GBI",
                               association_index = association.index)
  
  # Create vectors of host names
  host <- d.soc$Name
  n.host <- length(unique(host)) 
  .host.name <- unique(sort(host))
  
  # Construct and fill the symmetric association matrix
  social.matrix <- matrix(NA, nrow = n.host, ncol = n.host)
  colnames(social.matrix) <- .host.name
  rownames(social.matrix) <- .host.name
  
  indivs.ecoli <- sort(unique(d.coli$Name))
  indivs.social <- colnames(s.hwi)
  indivs.social[!(indivs.ecoli %in% indivs.social)]
  
  ### seeing what's in indivs and not in indivs_social
  for(i in 1:length(indivs.ecoli)){
    for(j in 1:length(indivs.ecoli)){
      indiv1 <- indivs.ecoli[i]
      indiv2 <- indivs.ecoli[j]
      indiv1_numb <- which(colnames(s.hwi) == indiv1)
      indiv1_numb <- as.numeric(indiv1_numb)
      indiv2_numb <- which(colnames(s.hwi) == indiv2)
      indiv2_numb <- as.numeric(indiv2_numb)
      hwi <- s.hwi[indiv1_numb, indiv2_numb]
      social.matrix[i, j] <- hwi
    }
  }
  

  
  ## Specify the cutoff value, below which individuals will not be considered as 'associating'
  for(i in 1:length(indivs.social)){
    for(j in 1:length(indivs.social)){ 
      
      if(s.hwi[i, j] < cutoff){
        social.matrix[i, j] <- 0
        
      }}}
  
  return(social.matrix)
  
}
