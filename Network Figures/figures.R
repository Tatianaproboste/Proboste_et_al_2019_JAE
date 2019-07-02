#####################################################
# code to create csv files to build the netwoks using Gephi. 

setwd("./Network_Simulations")

data.e <- read.csv ("./raw_data/V1_S3.csv", header = T,
                    stringsAsFactors = FALSE)
data.s <- read.csv("./raw_data/Survey_2014_2016_clipped.csv", header = T,
                   stringsAsFactors = FALSE)

Abrev <- read.csv("./raw_data/Abrev.csv", header = T,
                   stringsAsFactors = FALSE)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir('./Figures')


E.coli_matrix_2015_94 <- ecoli_matrix(data.ecoli=data.e, 
                        data.social=data.s, 
                        year= 2015, 
                        cutoff= 0.94)

Social_matrix_2015_005 <- social_matrix (data.ecoli= data.e,
                                       data.social=data.s,
                                       year=2015,
                                       cutoff=0.05, 
                                       association.index="HWI")


gephi <- data_gephi(ecoli=E.coli_matrix_2015_94,
                   social=Social_matrix_2015_005)

#write.csv(gephi$S.edges, file = "./Figures/GephiData/2015/social_edges_2015_005.csv")
#write.csv(gephi$S.node.table, file = "./Figures/GephiData/2015/social_node_2015_005.csv")
#write.csv(gephi$E.edges, file = "./Figures/GephiData/2014/ecoli_edges_2014_94.csv")
#write.csv(gephi$E.node.table, file = "./Figures/GephiData/2014/ecoli_node_2014_94.csv")

