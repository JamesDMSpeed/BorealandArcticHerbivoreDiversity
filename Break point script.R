
#Calculating Eucledian distnace to borders

#Load packages
require(raster) #Spatial
require(rgdal)#Spatial
require(sp)#Spatial
require(rasterVis)#Spatial
require(picante)#Diversity analysis
require(vegan)#Diversity analysis
require(ape)#View dendrograms
require(NbClust)#Find optimum number of clusters
require(mapdata)#Background world map
require(maptools)#Background world map
require(rgeos)#Crop biomes map
require(reshape2)#Data manipulation
require(dplyr)#find shapefile intercept


#########################Spearating the biomes
alltundra<-ecoregions[ecoregions$BIOME==11,]
tundra<-crop(alltundra,extent(-180,180,40,90))#Remove southern hemisphere
#Project
tundrapp<-spTransform(tundra,polarproj)



allboreal<-ecoregions[ecoregions$BIOME==6,]
boreal<-crop(allboreal,extent(-180,180,40,90))#Remove southern hemisphere
#Project
borealpp<-spTransform(boreal,polarproj)


plot(tundrapp)
plot(borealpp)


writeOGR(tundrapp,dsn="M:\\privat\\MSc_paper_1\\Tundra_Boreal shp files", layer= "tundrashp", driver = "ESRI Shapefile")
writeOGR(borealpp,dsn="M:\\privat\\MSc_paper_1\\Tundra_Boreal shp files", layer= "borealshp", driver = "ESRI Shapefile")

#########################Dissolving the shape files to find border between them

#onetundra<-unionSpatialPolygons(tundrapp,'tundra',)

onetundra<-gUnionCascaded(tundrapp, id = NULL)

oneboreal<-gUnionCascaded(borealpp, id = NULL)



#########################Finding the border

border<-gIntersection(onetundra, oneboreal, byid=FALSE, id=NULL,
              drop_lower_td=FALSE, unaryUnion_if_byid_false=TRUE, checkValidity=NULL)


#export border for checking. need to turn "LargeSpatial Line object into spatial dataframe
#border2<-SpatialLinesDataFrame(border, data.frame(id=1:length(border)))
#writeOGR(border2,dsn="M:\\privat\\MSc_paper_1\\Tundra_Boreal shp files", layer= "border", driver = "ESRI Shapefile")
#WORKS!

######################Distance to border

#raster to points
srp = as(speciesrichness,"SpatialPoints")

#distance
srd = gDistance(srp, border, byid=TRUE)

srdmin = apply(srd,2,min)

rastersrmin<-rasterize(srdmin, sr, field, fun='last', background=NA, mask=FALSE, update=FALSE, updateValue='all', na.rm=TRUE)











