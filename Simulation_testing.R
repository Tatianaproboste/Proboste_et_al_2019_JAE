# You shouldn't need to set a working directory as it will search in the source file directory
# 1. Load the data from the 'raw_data' folder
getwd()
gc()
rm(list = ls())
.rs.restartR()

#setwd("m:/Users/uqtprobo/Documents/Network_Simulations/Network_Simulations")
setwd("C:/Users/uqtprobo/Dropbox/PhD/Chapter_1/R/Network_Simulations/Network_Simulations")

data.e <- read.csv ("./raw_data/V1_S3.csv", header = T,
                    stringsAsFactors = FALSE)
data.s <- read.csv("./raw_data/Survey_2014_2016_clipped.csv", header = T,
                   stringsAsFactors = FALSE)
# Run simulation using the spatial data 
data.d.2014 <- as.data.frame(read.csv("./raw_data/distances_centroids_2014.csv"))
data.d.2014$X = NULL

data.d.2015 <- read.csv("./raw_data/distances_centroids_2015.csv")
data.d.2015$X=NULL


# 2. Source all functions needed to perform analyses. 
# In this case, we have stored all of these functions in the folder 'R', 
# and can source them all at once with a simple function
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}

sourceDir('./R')

# 3. Run simulations by 
#(a) allowing cutoff and social index values to vary,
# (b) constructing matrices and calculating network metrics, 
# and (c) running regression models. 
# This function produces summary statistics of all the necessary regressions, 
# but will take a while to run (calculation of network metrics is a bit slow). 
# !! NOTE, the following packages should be installed !!: 
# ecodist, pbapply, lme4, asnipe, sna, dplyr, purrr
#
# test <- run_network_sims(data.ecoli = data.e,
#                 data.social = data.s,
#                 n.simulations = 1,
#                 n.cores = 4,
#                 weighted.host = TRUE,
#                 spatial = TRUE,
#                 spatial_2014 = data.d.2014,
#                 spatial_2015 = data.d.2015)


test <- run_network_sims_spatial(data.ecoli = data.e,
                         data.social = data.s,
                         n.simulations = 1000,
                         n.cores = 10,
                         weighted.host = TRUE,
                         spatial = TRUE,
                         spatial_2014 = data.d.2014,
                         spatial_2015 = data.d.2015)


#save(test, file="./test_results/coli_space_1000.rda")

load("./test_results/test.unweigthed.v3.rda")
load("./test_results/test.weigthed.v3.rda")
test <-test.unweigthed.v3
test.w <-test.weigthed.v3


### What host cutoffs were tested?
test$host.cutoff.summary

### What host association indices were tested?
test$host.index.summary

### What ecoli cutoffs were tested?
test$ecoli.cutoff.summary



###___________________________Analysis 

# Table 2
#Results of tests for structural similarities between hosts' social networks and 
#commensal bacterial networks in 2014 and 2015. Coefficients represent 95% quantiles 
#of relationships between commensal bacterial network structure and hosts' social network
#structure, calculated using a multiple regression quadratic assignments procedure (MRQAP)
#and a multiple regression on distance matrices procedure (MRM) with 100 permutations.
#Estimates of coefficients, r squared value and the proportion of significant
#tests were generated from 1,000 simulations. Significant proportion is given in bold.

test$MRQAP.2014.summary
test$MRQAP.2015.summary
test$MRM.2014.summary
test$MRM.2015.summary

# Table 3
#Results of tests for structural similarities between bacterial commensal networks
#across years and hosts' social networks across years. Coefficients represent 
#relationships between networks' structures between the years, calculated using a
#multiple regression quadratic assignments procedure (MRQAP) and a multiple regression 
#on distance matrices procedure (MRM) with 100 permutations. Estimates of coefficients,
#r squared and the proportion of significant test were generated from 1,000 simulations.

test$MRQAP.ecoli.summary
test$MRQAP.host.summary
test$MRM.ecoli.summary
test$MRM.host.summary

# Individual-based comparison
test$degree.coef.summary


###__________________________ Suplementary Material 

# Table S1. 
#a)
#Comparison between structural similarities between hosts' social 
#networks (weighted and unweighted social network) and commensal bacterial networks
#in 2014 and 2015. Coefficients represent 95% quantiles of relationships between 
#commensal bacterial network structure and hosts' social network structure, calculated
#using: a) multiple regression quadratic assignments procedure (MRQAP); and b) multiple 
#regression on distance matrices procedure (MRM) with 100 permutations. Estimates
#of coefficients, r squared value and the proportion of tests that were significant were
#generated from 1,000 simulations. Significant proportion is given in bold.

#unweigthed 2014 and 2015
test$MRQAP.2014.summary
test$MRQAP.2015.summary

#weigthed 2014 and 2015
test.w$MRQAP.2014.summary
test.w$MRQAP.2015.summary

#b)
test$MRM.2014.summary
test$MRM.2015.summary

test.w$MRM.2014.summary
test.w$MRM.2015.summary


#Table S2. 
#Estimates of fixed effects from the linear mixed model to predict
#social network a) betweenness and b) eigenvector centrality. Individual ID was 
#included as a random effect. Regression estimates for each of the predictor
#variables varied substantially across simulations with 95% and 
#85% quantiles of regression coefficients overlapping zero for nearly all predictors

#a)
test$betweeness.coef.summary

#b)
test$eigencent.coef.summary


#Table S3. 
#Tests of community relationships between host social networks and 
#commensal bacterial networks in 2014 and 2015. Coefficients represent relationships 
#between commensal bacterial network structure and host social network structure, 
#calculated using a multiple regression quadratic assignments procedure (MRQAP) 
#with 1,000 permutations.




#Table S4. 
#Summary of network-level metrics for a commensal bacterial network 
#constructed using a 90% and 94% similarity cut-off and social network build
#based on HWI index and a threshold of 0.05 for 2014 and 2015. Density and global
#clustering have been calculated using function in "sna" package and assortativity 
#and individual degree has been calculated using function 
#in "igraph" package. We show the mean degree and the stander deviation (SD).  

#upload the matrix
e.coli.matrix.2014.0.90 <- read.csv("./raw_data/e_coli_matrix_2014_090.csv")
e.coli.matrix.2014.0.90 <- e.coli.matrix.2014.0.90[ ,2:30]
e.coli.matrix.2015.0.90 <- read.csv("./raw_data/e_coli_matrix_2015_090.csv")
e.coli.matrix.2015.0.90 <- e.coli.matrix.2015.0.90[, 3:63]
e.coli.matrix.2014.0.94 <- read.csv("./raw_data/e_coli_matrix_2014_094.csv")
e.coli.matrix.2014.0.94 <- e.coli.matrix.2014.0.94 [,2:30]
e.coli.matrix.2015.0.94 <- read.csv("./raw_data/e_coli_matrix_2015_094.csv")
e.coli.matrix.2015.0.94 <- e.coli.matrix.2015.0.94 [,2:62]
social.matrix.2014<- read.csv("./raw_data/social.matrix.2014.csv")
social.matrix.2014 <- social.matrix.2014[,2:30]
social.matrix.2015 <- read.csv("./raw_data/social.matrix.2015.csv")
social.matrix.2015 <- social.matrix.2015[,2:62]

#generated the networks
g.2014.0.90<-graph_from_adjacency_matrix (as.matrix(e.coli.matrix.2014.0.90),
                                          weighted = NULL, mode ="undirected", diag = FALSE)

g.2015.0.90<-graph_from_adjacency_matrix (as.matrix(e.coli.matrix.2015.0.90),
                                          weighted = NULL, mode ="undirected", diag = FALSE) # E coli

g.2014.0.94<-graph_from_adjacency_matrix (as.matrix(e.coli.matrix.2014.0.94),
                                          weighted = NULL, mode ="undirected", diag = FALSE)

g.2015.0.94<-graph_from_adjacency_matrix (as.matrix(e.coli.matrix.2015.0.94),
                                          weighted = NULL, mode ="undirected", diag = FALSE) # E coli

s.2014<-graph_from_adjacency_matrix (as.matrix(social.matrix.2014),
                                     weighted = TRUE, mode ="undirected", diag = FALSE)

s.2015<-graph_from_adjacency_matrix (as.matrix(social.matrix.2015),
                                     weighted = TRUE, mode ="undirected", diag = FALSE)

# Ecoli 2014 90%
density.2014.90 <- sna::gden(e.coli.matrix.2014.0.90)
dens.2014.90 <- igraph::edge_density(g.2014.0.90)
transitivity.2014.90 <- sna::gtrans(e.coli.matrix.2014.0.90) # global clustering coefficient (transitivity):
assortativity.2014.90 <- igraph::assortativity_degree(g.2014.0.90, directed = F)
deg.2014.90<- igraph::degree(g.2014.0.90)
summary(deg.2014.90)
sd(deg,deg.2014.90)
plot(g.2014.0.90)

#ecoli 2015 90%
density.2015.90 <- sna::gden(e.coli.matrix.2015.0.90)
dens.2015.90 <- igraph::edge_density(g.2015.0.90)
transitivity.2015.90 <- sna::gtrans(e.coli.matrix.2015.0.90)
assortativity.2015.90 <- igraph::assortativity_degree(g.2015.0.90, directed = F)
plot(g.2015.0.90)
deg.2015.90<- igraph::degree(g.2015.0.90)
summary(deg.2015.90)
sd(deg.2015.90)

#ecoli 2014 94%
density.2014.94 <- sna::gden(e.coli.matrix.2014.0.94)
dens.2014.94 <- igraph::edge_density(g.2014.0.94)
transitivity.2014.94 <- sna::gtrans(e.coli.matrix.2014.0.94)
assortativity.2014.94 <- igraph::assortativity_degree(g.2014.0.94, directed = F)
plot(g.2014.0.94)
deg.2014.94<- igraph::degree(g.2014.0.94)
summary(deg.2014.94)
sd(deg.2014.94)

#ecoli 2015 94%
density.2015.94 <- sna::gden(e.coli.matrix.2015.0.94)
dens.2015.94 <- igraph::edge_density(g.2015.0.94)
transitivity.2015.94 <- sna::gtrans(e.coli.matrix.2015.0.94)
assortativity.2015.94 <- igraph::assortativity_degree(g.2015.0.94, directed = F)
plot(g.2015.0.94)
deg.2015.94<- igraph::degree(g.2015.0.94)
summary(deg.2015.94)
sd(deg.2015.94)

#social 2014
s.density.2014 <- sna::gden(social.matrix.2014, mode = "graph", ignore.eval = T)  # revisar si egnore.eval esta bien siendo "True"
s.dens.2014 <- igraph::edge_density(s.2014)
s.transitivity.2014 <- sna::gtrans(social.matrix.2014)
s.mean_p.2014 <- igraph::mean_distance(s.2014)
assortativity.s.2014 <- igraph::assortativity_degree(s.2014, directed = F)
plot(s.2014)
deg.social.2014 <- igraph::degree(s.2014)
summary(deg.social.2014)
sd(deg.social.2014)

#social 2015
s.density.2015 <- sna::gden(social.matrix.2015, mode = "graph", ignore.eval = T)  # revisar si egnore.eval esta bien siendo "True"
s.dens.2015 <- igraph::edge_density(s.2015)
s.transitivity.2015 <- sna::gtrans(social.matrix.2015)
s.mean_p.2015 <- igraph::mean_distance(s.2015)
assortativity.s.2015 <- igraph::assortativity_degree(s.2015, directed = F)
plot(s.2015)
deg <- sna::degree(social.matrix.2015)
deg.social.2015 <- igraph::degree(s.2015)
summary(deg.social.2015)
sd(deg.social.2015)


# for the Figures, we used the code "Figure" which create the matrix in the format 
# needed to work with Gephi.


#Figure S5. 
#Degree coefficient for commensal bacterial
#network built based on a 94% similarity and 90% similarity 
#generated using 1,000 simulations.


library(dplyr)
library(ggplot2)
plot.data <- data.frame(Ecoli.cutoff = as.factor(test$ecoli.cutoff.raw),
                        Coef = test$degree.coef.raw %>%
                          filter(Parameter=='Ecoli.degree'))

e94 <- subset(plot.data, Ecoli.cutoff==0.94)
e90 <- subset(plot.data, Ecoli.cutoff==0.90)

hist(e94$Coef.Estimate, main = "94% similarity", ylab = "Frequency", xlab = "Degree Coefficient")
hist(e90$Coef.Estimate, main = "94% similarity", ylab = "Frequency", xlab = "Degree Coefficient")


#Table S6. 
#Results of tests for structural similarities between hosts' spatial 
#networks and commensal bacterial networks in 2014 and 2015. Coefficients represent
#95% quantiles of relationships between commensal bacterial network structure and hosts' 
#social network structure, calculated using a multiple regression quadratic assignments
#procedure (MRQAP) with 100 permutations. Estimates of coefficients, r squared value and
#the proportion of significant tests were generated from 100 simulations. Significant
#proportion is given in bold.


#check spatial matrix data
test$MRQAP.2014.dist.summary
test$MRQAP.2015.dist.summary

test$MRM.2014.dist.summary
test$MRM.2015.dist.summary

# Check a few results, such as 95% CIs of coefficients for the degree regression
test$degree.coef.summary
test.unweigthed.v3$degree.coef.summary

# or for the eigencentrality regression
test.unweigthed.v3$eigencent.coef.summary

# also check some of the matrix regression results 
# (prop.significant indicates the proportion of total models 
# that were statistically 'significant')
test.unweighted$MRQAP.2015.raw
test.unweigthed.v3$MRQAP.2015.summary

test$MRM.2014.summary

test$MRM.2014.summary
test$MRQAP.2015.summary
test$MRM.2015.summary

test.unweigthed.v3$MRM.ecoli.summary
test$MRQAP.ecoli.summary

test$MRQAP.host.summary
test$MRM.host.summary

test$host.index.raw

# What host cutoffs were tested?
test$host.cutoff.summary

#What host association indices were tested?
test$host.index.summary

# What ecoli cutoffs were tested?
test$ecoli.cutoff.summary

