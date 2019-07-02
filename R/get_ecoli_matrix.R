get_ecoli_matrix = function(data.ecoli, data.social, year, cutoff){
  
  #### Make sure year and cutoff specifications are allowed ####
  if(!(year %in% c(2014, 2015))){
    # if year isn't legal, issue an error
    stop('Year must be a numeric value matching either 2014 or 2015')
  }
  
  if(!(cutoff %in% c(0.90, 0.94))){
    stop('cutoff must be a numeric value matching either 0.90 or 0.94')
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
  
  #### Create vectors of host.A and strain names ####
  host <- d.coli$Name
  
  # Strain names are assigned based on the specified similarity cutoff value
  if(cutoff == 0.90){
    strain <- d.coli$S3
  }
  
  if(cutoff == 0.94){
    strain <- d.coli$S1
  }
  
  # Sort host.A names to make identical row and column orders in the host.AXhost.A matrix
  n.host <- length(unique(host)) 
  .host.name <- unique(sort(host))
  
  #### Creating the ecoli strain sharing matrix ####
  # Create an empty host.AXhost.A matrix to fill with a for loop
  e.coli.matrix <- matrix(NA, nrow = n.host, ncol = n.host)
  
  # For loop: find unique strains for each host.A pair & calculate how many are shared
  for (i in 1:n.host){
    for (l in 1:n.host){
      host1.strain <- strain[which(host == .host.name[i])]
      host2.strain <- strain[which(host == .host.name[l])]
      e.coli.matrix[i, l] <- length(which(host1.strain %in% host2.strain))
    }
  }
  
  # Give names to rows and columns
  colnames(e.coli.matrix) <- rownames(e.coli.matrix) <- .host.name
  
  # diagonal = 0
  diag(e.coli.matrix) <- 0
  
  # weithed network 
  #e.coli.matrix <- e.coli.matrix # this doesn't work
  
  #binary Network
  e.coli.matrix[which(e.coli.matrix > 0)] <- 1 
  
  
  return(e.coli.matrix)
  
}