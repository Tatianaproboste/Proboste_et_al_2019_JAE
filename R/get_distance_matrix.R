
get_distance_matrix = function(data.social, data.ecoli, year, cutoff, Area){

  
#### 1. First part of Digiroo file ####
library(rgeos)
library(raster)
library(dismo)
library(RgoogleMaps)
library(spatstat)
library(spdep)
library(adehabitatHR)


setwd("C:/Users/uqtprobo/Dropbox/PhD/Chapter_1/R/Network_Simulations/Network_Simulations")
path <- ("C:/Users/uqtprobo/Dropbox/PhD/Chapter_1/R/Network_Simulations/Network_Simulations/Spatia_data")

  # Specify projections used
WGSproj <- "+init=epsg:4326"
GDAproj <- "+init=epsg:28355"
mercproj <-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
  
  
library(rgdal)
StudyArea <- readOGR(path, "Study area forest Erase EB")
StudyAreaP <- spTransform(StudyArea,CRS(GDAproj))

#data
data.social$Year <- as.numeric(as.character(data.social$Year))
survey_year <- data.social


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

d.soc <- d.soc[d.soc$Name %in% d.coli$Name, ]
d.coli <- d.coli[d.coli$Name %in% d.soc$Name, ]

#### Create vectors of host.A and strain names ####
host <- d.coli$Name
n.host <- length(unique(host)) 
.host.name <- unique(sort(host))
social.matrix <- matrix(NA, nrow = n.host, ncol = n.host)

# Give names to rows and columns
colnames(social.matrix) <- rownames(social.matrix) <- .host.name

if(year == 2014){
 survey_year <- data.social[which(data.social$Year==2014),]
 
  } else{
 survey_year <- data.social[which(data.social$Year==2015),]
  }


#### Clean the data bases with only the individuals that appear in both datasets ####
survey_year <- survey_year[survey_year$Name %in% data.ecoli$Name,]

survey_year$Latitude <- as.numeric(as.character(survey_year$Latitude))
survey_year$Longitude <- as.numeric(as.character(survey_year$Longitude))

survey_NA_remove <- survey_year[which(!is.na(survey_year$Latitude)),]
survey_NA_remove <- survey_NA_remove[which(!is.na(survey_NA_remove$Longitude)),]

Roo <- data.frame(Name=survey_NA_remove$Name,Latitude=survey_NA_remove$Latitude,Longitude=survey_NA_remove$Longitude,SurveyNumber=survey_NA_remove$SurveyNumber)
Roo$Latitude <- as.numeric(as.character(Roo$Latitude))
Roo$Longitude <- as.numeric(as.character(Roo$Longitude))


locs <- Roo[,c("Longitude","Latitude")]
locssp<- SpatialPoints(locs,CRS(WGSproj))
locsspdf<- Roo[,c("Name","Longitude","Latitude")]
locsspdf<- Roo
coordinates(locsspdf)<- ~Longitude+Latitude
proj4string(locsspdf)<- CRS(GDAproj)

###Load ArcGIS created .shp polygon and reproject in correct coordinate system
library(rgdal)
StudyAreaP <- spTransform(Area,CRS(GDAproj))

#3 Generate buffer of study area and extract only those points falling within the buffer

#Firstly we need to tell R there are holes in the polygon and assign a comment 
#saying which hole belongs to which polygon

library(rgeos)
lapply(slot(StudyAreaP, "polygons"), comment)
any(c(sapply(slot(StudyAreaP, "polygons"),
             function(x) sapply(slot(x, "Polygons"), slot, "hole"))))
slot(createSPComment(StudyAreaP), "polygons")

StudyAreaP1 <- StudyAreaP
Pls <- slot(StudyAreaP1, "polygons")
pls <- slot(Pls[[1]], "Polygons")
pls1 <- lapply(pls, function(p) {slot(p, "hole") <- FALSE; return(p)})
slot(Pls[[1]], "Polygons") <- pls1
slot(StudyAreaP1, "polygons") <- Pls
any(c(sapply(slot(StudyAreaP1, "polygons"),
             function(x) sapply(slot(x, "Polygons"), slot, "hole"))))

lapply(slot(createSPComment(StudyAreaP1), "polygons"), comment)
gIsValid(StudyAreaP1, reason=TRUE)
Pls <- slot(StudyAreaP1, "polygons")
lapply(slot(StudyAreaP1, "polygons"), comment)
gIsValid(StudyAreaP1, reason=TRUE)

##Ensure buffer is not too large
StudyAreaPB <- gBuffer(StudyAreaP1,width=10)

if(gIsEmpty(StudyAreaPB)){
  stop("[RD] Buffer size too large")
}else{
  StudyAreaB <- spTransform(StudyAreaPB,CRS(WGSproj))
}


#Extract points falling within Study Area Polygon
locsspdf2<- Roo
coordinates(locsspdf2)<- ~Longitude+Latitude
proj4string(locsspdf2)<- CRS(WGSproj)
NewRoo_test2 <- locsspdf2[!is.na(over(locsspdf2, StudyAreaB)),]
NewRoo <- data.frame(NewRoo_test2)

# get rid of multiple sightings in a survey
NewRoo2 <- NewRoo[,c("Name","SurveyNumber")]
NewRoo2 <- NewRoo[!duplicated(NewRoo2),]


#### 2. Get homerange overlap ####
# Select for Individuals with more than 5 sightings

Roof <- tapply(rep(1,nrow(NewRoo2)),NewRoo2$Name,sum)
AllRoofreq <- data.frame(Name=names(Roof),Frequency=as.vector(Roof),
                         stringsAsFactors=FALSE)
AllRoofreq <- AllRoofreq[order(AllRoofreq[,2],decreasing = TRUE),]

RooDatTable <- data.frame(X=1:nrow(AllRoofreq),AllRoofreq)


# Here we can filter to keep only individulas where we have more 8+ locations, but for our analysis we 
# are keeping all the individuals.

NoRoos <- as.character(AllRoofreq[which(AllRoofreq[,2]>8),1]) # only those individuals where we have 9+ locations
NewRoo2 <- NewRoo2[as.character(NewRoo2[,"Name"]) %in% NoRoos,]
NewRoo2$Name <- as.character(NewRoo2$Name)
tapply(rep(1,nrow(NewRoo2)),NewRoo2$Name,sum) #All data


# Estimate UDs for selection of individuals
# Ensure data is in the correct coordinate system
NewRoosub <- NewRoo2[c("Longitude","Latitude","Name")]
coordinates(NewRoosub) <- c("Longitude","Latitude")
NewRoosub$Name <- as.character(NewRoosub$Name)
proj4string(NewRoosub) <- CRS(WGSproj)

#Convert all objects to a projected coordinate system (helps UD calculations)
NewRoosubP <- spTransform(NewRoosub,CRS(GDAproj))

## Estimation of UD for the 5 animals
ud <- kernelUD(NewRoosubP, h="href",same4all=TRUE) # works


## Calculation of the 95 percent home range in hectares
(ver95 <- getverticeshr(ud, 95, unin = c("m"),unout = c("ha")))

library(rgeos)
trueCentroids = gCentroid(ver95,byid=T)
plot(ver95)
points(coordinates(ver95),pch=1)
points(trueCentroids,pch=2,col="red")

# transform centroids
#proj4string(trueCentroids) <- CRS(WGSproj)
trueCentroids <- spTransform(trueCentroids,CRS(GDAproj))

# get distances between centroids for all individuals
distances_centroids <- gDistance(trueCentroids,trueCentroids,byid=T)


# Make sure that the distance_centroids matrix has the same names social matrix
indivs_distances_centroids <- rownames(distances_centroids)
length(indivs_distances_centroids)
d.dist <- as.data.frame(distances_centroids)


#invividuals 

indivs_social_matrix <- .host.name

distances_centroids <- data.frame(matrix(NA,nrow=length(indivs_social_matrix),ncol=length(indivs_social_matrix)))
colnames(distances_centroids) <- indivs_social_matrix
rownames(distances_centroids) <- indivs_social_matrix

for(i in 1:length(indivs_distances_centroids)){
  for(j in 1:length(indivs_distances_centroids)){
    indiv1 <- indivs_distances_centroids[i]
    indiv2 <- indivs_distances_centroids[j]
    distance <- d.dist[i,j]
    
    indiv1_number <- which(colnames(distances_centroids)==indiv1)
    indiv2_number <- which(colnames(distances_centroids)==indiv2)
    distances_centroids[indiv1_number,indiv2_number] <- distance
  }
}

return(distances_centroids)
}

#write.csv(distances_centroids, file="./Spatial_matrix_analysis/distance.centroids_2014.csv")

