#Analysing boreal and arctic herbivore diversity

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
require(FD)


# Set up ------------------------------------------------------------------

#Polar projection
polarproj<-'+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 '

#Country outlines
noreco_countries <- c('NO', 'SE', 'FI','RU','CA','IS','GL','SJ','MN','JP') 
noreco_shp1 <- do.call("bind", lapply(noreco_countries, function(x)  raster::getData('GADM', country=x, level=0)))
#Only Alaska from USA #had to add raster::getData due to error in reading the argument US
us<-raster::getData('GADM',country='US',level=1)
alaska<-us[us$NAME_1=='Alaska',]
#Bind
noreco_shp<-bind(noreco_shp1,alaska)
#Project
noreco_shppp<-spTransform(noreco_shp,CRS=crs(polarproj))

#plot(noreco_shppp)


#Biomes
ecoregions<-readOGR('Biomes','wwf_terr_ecos')
tundraboreal<-ecoregions[ecoregions$BIOME==6|ecoregions$BIOME==11,]
northernecosystems<-crop(tundraboreal,extent(-180,180,40,90))#Remove southern hemisphere
#Project
northernecosystemspp<-spTransform(northernecosystems,polarproj)
#plot(northernecosystemspp)

#Regions
regions<-readOGR('CMEC regions & realms','Regions')
regionspp<-spTransform(regions,polarproj)
regionAB<-crop(regionspp,northernecosystemspp)

#Productivity http://www.ntsg.umt.edu/project/modis/mod17.php
globnpp_url<-('https://ntnu.box.com/shared/static/yo8wmtj09k1vqop2gpstlsv889wff2mp.tif')
download.file(globnpp_url,'GlobNPP_AnnMean00_15.tif')
globnpp<-raster('GlobNPP_AnnMean00_15.tif')
#Rescale
globnpps<-globnpp*0.1#http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/GeoTIFF/MOD17A3/readme.txt 
plot(globnpps)#g/m2/yr

# Species data ------------------------------------------------------------

filelist<-list.files('SpeciesDistributionData',full.names = T)
filelist
#Stack up
herbivore_dataset<-stack(filelist)

#Tidy names
names(herbivore_dataset)<-substring(names(herbivore_dataset),3)
names(herbivore_dataset)

plot(herbivore_dataset[[1]])

#Project
herbivore_dataset<-projectRaster(herbivore_dataset,crs=polarproj)
herbivore_dataset


#---Different Extent Issue fixed by stacking stacks-not classy but works---

filelist2<-list.files('DifferentExtent',full.names = T)
filelist2
#Stack up
herbivore_dataset2<-stack(filelist2)

#Tidy names
names(herbivore_dataset2)<-substring(names(herbivore_dataset2),3)
names(herbivore_dataset2)

plot(herbivore_dataset2[[1]])

#Project
herbivore_dataset2<-projectRaster(herbivore_dataset2,crs=polarproj)
herbivore_dataset2

#setExtent(herbivore_dataset2,herbivore_dataset)

#Stack the two
herbivore_dataset3<-stack(herbivore_dataset,crop(herbivore_dataset2,herbivore_dataset))
#Working rasterstack = herbivore_dataset3
writeRaster(herbivore_dataset3,'FinalRanges/',bylayer=T,suffix=names(herbivore_dataset3))
lf<-list.files('FinalRanges/',full.names = T,pattern = 'grd')
herbivore_dataset3<-stack(lf)

#Trait data
traits<-read.csv('FunctionalClassification/TraitTableFeb2019.csv',header=T)
#Make column with species names formatted same as raster names
traits$SpeciesNames<-gsub(" ",".",traits$Binomial)
#Check these
traits$SpeciesNames%in%names(herbivore_dataset3)
names(herbivore_dataset3)%in%traits$SpeciesNames

#Make a stack for bodymass
bodymassstack<-herbivore_dataset3
graminoidstack<-herbivore_dataset3
shrubstack<-herbivore_dataset3
forbstack<-herbivore_dataset3
mossstack<-herbivore_dataset3
lichenstack<-herbivore_dataset3
littersizestack<-herbivore_dataset3
#Replace with trait values
for (i in 1:nlayers(herbivore_dataset3)){
bodymassstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$body_mass[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
graminoidstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_graminoid[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
shrubstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_shrub[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
forbstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_forb[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
mossstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_moss[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
lichenstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_lichen[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
littersizestack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$Litter_clutch_size[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
}

meantraitstack<-stack(calc(bodymassstack,mean,na.rm=T),
                      calc(graminoidstack,mean,na.rm=T),
                      calc(shrubstack,mean,na.rm=T),
                      calc(forbstack,mean,na.rm=T),
                      calc(mossstack,mean,na.rm=T),
                      calc(lichenstack,mean,na.rm=T),
                      calc(littersizestack,mean,na.rm=T))
names(meantraitstack)<-c('Body Mass','Graminoids_mean','Shrubs_mean','Forbs_mean','Bryophytes_mean','Lichens_mean','Litter/Clutch size')
plot(meantraitstack)

#Relative diet contributions since diets do not sum to the same across species
#Total diet
traits$TotalDiet<-traits$diet_item_forb+traits$diet_item_graminoid+traits$diet_item_shrub+traits$diet_item_moss+traits$diet_item_lichen

rel_graminoidstack<-herbivore_dataset3
rel_shrubstack<-herbivore_dataset3
rel_forbstack<-herbivore_dataset3
rel_mossstack<-herbivore_dataset3
rel_lichenstack<-herbivore_dataset3

for (i in 1:nlayers(herbivore_dataset3)){
rel_graminoidstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_graminoid[traits$SpeciesNames==names(herbivore_dataset3[[i]])]/traits$TotalDiet[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
rel_shrubstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_shrub[traits$SpeciesNames==names(herbivore_dataset3[[i]])]/traits$TotalDiet[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
rel_forbstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_forb[traits$SpeciesNames==names(herbivore_dataset3[[i]])]/traits$TotalDiet[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
rel_mossstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_moss[traits$SpeciesNames==names(herbivore_dataset3[[i]])]/traits$TotalDiet[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
rel_lichenstack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,traits$diet_item_lichen[traits$SpeciesNames==names(herbivore_dataset3[[i]])]/traits$TotalDiet[traits$SpeciesNames==names(herbivore_dataset3[[i]])]))
} 

#Accumulated diet types
graminoidsum<-calc(rel_graminoidstack,sum,na.rm=T)
forbsum<-calc(rel_forbstack,sum,na.rm=T)
shrubsum<-calc(rel_shrubstack,sum,na.rm=T)
mosssum<-calc(rel_mossstack,sum,na.rm=T)
lichensum<-calc(rel_lichenstack,sum,na.rm=T)

dietstacksum<-stack(graminoidsum,forbsum,shrubsum,mosssum,lichensum)
names(dietstacksum)<-c('Graminoid_sum','Forbs_sum','Shrubs_sum','Moss_sum','Lichens_sum')
levelplot(dietstacksum,scales=list(draw=F),names.attr=names(dietstacksum),par.settings='YlOrRdTheme')+
  layer(sp.polygons(bPolslaea))+
  layer(sp.polygons(arczones_laea[arczones_laea$ZONE_==0,],lty=2))


#Functional classification coordinates
famdcoord<-read.csv('FunctionalClassification/FAMD_coords.csv',header=T)
famdcoord$X<-gsub(" ",".",famdcoord$X)
pc1stack<-herbivore_dataset3
pc2stack<-herbivore_dataset3
for (i in 1:nlayers(herbivore_dataset3)){
  pc1stack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,famdcoord$Dim.1[famdcoord$X==names(herbivore_dataset3[[i]])]))
  pc2stack[[i]]<-reclassify(herbivore_dataset3[[i]],cbind(1,famdcoord$Dim.2[famdcoord$X==names(herbivore_dataset3[[i]])]))
}

pcstack<-stack(calc(pc1stack,mean,na.rm=T),
                      calc(pc2stack,mean,na.rm=T))
names(pcstack)<-c('PC1','PC2')
plot(pcstack)


#Distance to Biome Boundary
distBiomeBound<-raster('Biomes/DistnaceBiomeBoundary')
diverge0 <- function(p, ramp) {
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(rev(brewer.pal(11, ramp))))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
  p
}
pBB<-levelplot(distBiomeBound,margin=F)
diverge0(pBB,'RdBu')

#Vegetation Cover #See here https://github.com/JamesDMSpeed/AB_FunctionalTraits/blob/master/HerbivoreTraitData.R
shrubcov<-raster('VegCover/ShrubCover.tif')
treecov<-raster('VegCover/treeCover.tif')
vegcov<-stack(shrubcov,treecov)
plot(vegcov)


#Stack all up
AllVars<-stack(distBiomeBound,vegcov,pcstack,meantraitstack,dietstacksum)
names(AllVars)[1]<-'DistanceBiomeBoundary'
names(AllVars)

pairs(AllVars[[c(1:6,12)]])
pairs(AllVars[[c(1:5,7:11)]])
AllVars_ex<-extract(AllVars,1:ncell(AllVars))

#Functional diversity package
traitselection<-traits[,c(12,15,18,21,24,27,30,33,36,42,45,48,51,54,57,60,63)]#Not human managed
rownames(traitselection)<-traits$SpeciesNames
traitdf<-traitselection[order(row.names(traitselection)),]#Order alphabetically
speciesdat<-extract(herbivore_dataset3,1:ncell(herbivore_dataset3),df=T)
speciesdata<-speciesdat[,2:ncol(speciesdat)]#Removing ID column
speciesmat<-speciesdata[,order(names(speciesdata))]#Order alphabeticaaly to match traits
speciesCommunities<-speciesmat[rowSums(speciesmat,na.rm=T)>0,]#Remove empty rows
speciescoms<-speciesCommunities[,colSums(speciesCommunities,na.rm=T)>0]#Remove empty columns
traitdf1<-traitdf[rownames(traitdf)%in%names(speciescoms),]#Remove empty species from trait data too

funcdiv<-dbFD(traitdf1,speciescoms,corr='lingoes')
