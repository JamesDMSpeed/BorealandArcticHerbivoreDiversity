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
require(FD)#FD indices
require(piecewiseSEM)#SEM
library(segmented)#Segmented regs



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

# filelist<-list.files('SpeciesDistributionData',full.names = T)
# filelist
# #Stack up
# herbivore_dataset<-stack(filelist)
# 
# #Tidy names
# names(herbivore_dataset)<-substring(names(herbivore_dataset),3)
# names(herbivore_dataset)
# 
# plot(herbivore_dataset[[1]])
# 
# #Project
# herbivore_dataset<-projectRaster(herbivore_dataset,crs=polarproj)
# herbivore_dataset
# 
# 
# #---Different Extent Issue fixed by stacking stacks-not classy but works---
# 
# filelist2<-list.files('DifferentExtent',full.names = T)
# filelist2
# #Stack up
# herbivore_dataset2<-stack(filelist2)
# 
# #Tidy names
# names(herbivore_dataset2)<-substring(names(herbivore_dataset2),3)
# names(herbivore_dataset2)
# 
# plot(herbivore_dataset2[[1]])
# 
# #Project
# herbivore_dataset2<-projectRaster(herbivore_dataset2,crs=polarproj)
# herbivore_dataset2
# 
# #setExtent(herbivore_dataset2,herbivore_dataset)
# 
# #Stack the two
# herbivore_dataset3<-stack(herbivore_dataset,crop(herbivore_dataset2,herbivore_dataset))
# #Working rasterstack = herbivore_dataset3
# writeRaster(herbivore_dataset3,'FinalRanges/',bylayer=T,suffix=names(herbivore_dataset3))
# 


lf<-list.files('FinalRanges/',full.names = T,pattern = 'grd')
herbivore_dataset3<-stack(lf)

#Trait data
traits<-read.csv('FunctionalClassification/traits_cor.csv',header=T)
#Make column with species names formatted same as raster names
traits$SpeciesNames<-gsub(" ",".",traits$Binomial)
#Check these
traits$SpeciesNames%in%names(herbivore_dataset3)
names(herbivore_dataset3)%in%traits$SpeciesNames

#Make a stack for variables
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

#SD trait stacks
sdtraitstack<-stack(calc(bodymassstack,sd,na.rm=T),
                      calc(graminoidstack,sd,na.rm=T),
                      calc(shrubstack,sd,na.rm=T),
                      calc(forbstack,sd,na.rm=T),
                      calc(mossstack,sd,na.rm=T),
                      calc(lichenstack,sd,na.rm=T),
                      calc(littersizestack,sd,na.rm=T),
                      calc(pc1stack,sd,na.rm=T),
                      calc(pc2stack,sd,na.rm=T))
names(sdtraitstack)<-c('Body Mass_sd','Graminoids_sd','Shrubs_sd','Forbs_sd','Bryophytes_sd','Lichens_sd','Litter/Clutch size_sd','PC1_sd','PC2_sd')
plot(sdtraitstack)

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


#Functional diversity package
#traitselection<-traits[,c(12,15,18,21,24,27,30,33,36,42,45,48,51,54,57,60,63)]#Not human managed
traitselection<-traits[,c(7:16,18:25)]
rownames(traitselection)<-traits$SpeciesNames
traitdf<-traitselection[order(row.names(traitselection)),]#Order alphabetically
speciesdat<-extract(herbivore_dataset3,1:ncell(herbivore_dataset3),df=T)
#speciesdata<-speciesdat[,2:ncol(speciesdat)]#Removing ID column
speciesmat<-speciesdat[,order(names(speciesdat)[2:length(names(speciesdat))])]#Order alphabeticaaly to match traits
speciesCommunities<-speciesmat[rowSums(speciesmat[,2:ncol(speciesmat)],na.rm=T)>0,]#Remove empty rows #Omit ID column
speciescoms<-speciesCommunities[,colSums(speciesCommunities,na.rm=T)>0]#Remove empty columns
s2<-speciescoms[,2:ncol(speciescoms)]
traitdf1<-traitdf[rownames(traitdf)%in%names(s2),]#Remove empty species from trait data too

#funcdiv<-dbFD(x=traitdf1[,c(1:2,5,11:15)],a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F)
#funcdiv<-dbFD(x=traitdf1[,c(3:4)],a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F)
funcdiv<-dbFD(x=traitdf1,a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F,CWM.type = 'all')#Abundance of all ordinal traits.


#Fill raster with FD indices
dfFD<-data.frame(ID=speciescoms$ID,FDiv=funcdiv$FDiv,FRic=funcdiv$FRic,FEve=funcdiv$FEve,FDis=funcdiv$FDis,RaoQ=funcdiv$RaoQ)#Dataframe with Fdiv indices and cell numbers
dfFDCWM<-cbind(dfFD,funcdiv$CWM)
df2<-data.frame(cell=seq(1,ncell(herbivore_dataset3[[1]]),1))#Dataframe with all (inc NA) cells
m1<-merge(df2,dfFDCWM,by.x='cell',by.y='ID',all.x=T)#Merge together
dummyras<-herbivore_dataset3[[1]]
fdivras<-setValues(dummyras,m1$FDiv)
fricras<-setValues(dummyras,m1$FRic)
feveras<-setValues(dummyras,m1$FEve)
fdisras<-setValues(dummyras,m1$FDis)
raoQras<-setValues(dummyras,m1$RaoQ)
fdivstack<-stack(fdivras,fricras,feveras,fdisras,raoQras)
names(fdivstack)<-c('FDiv','FRic','FEve','FDis','RaoQ')
plot(fdivstack)

dfCWM<-data.frame(ID=speciescoms$ID,Body.Mass_cwm=funcdiv$CWM$body_mass)#Dataframe with CWM indices and cell numbers
m2<-merge(df2,dfCWM,by.x='cell',by.y='ID',all.x=T)#Merge together
Body.Mass_cwm<-setValues(dummyras,m2$Body.Mass_cwm)
#Check consistency with CWM calculated direct (outside of FD package)
plot(Body.Mass_cwm)
plot(meantraitstack$Body.Mass)



#Temperature & NDVI etc 

#Global BioClim 2.5deg
bioclimdat<-raster::getData('worldclim',var='bio',res=2.5)
bioclimlaea<-projectRaster(bioclimdat,herbivore_dataset3)
bioclimlaea_m<-mask(crop(bioclimlaea,herbivore_dataset3[[1]]),herbivore_dataset3[[1]])
plot(bioclimlaea_m[[1]])

#Productivity
globnpp_url<-('https://ntnu.box.com/shared/static/yo8wmtj09k1vqop2gpstlsv889wff2mp.tif')
download.file(globnpp_url,'Biomes/GlobNPP_AnnMean00_15.tif',mode='wb')
globnpp<-raster('Biomes/GlobNPP_AnnMean00_15.tif')
#Rescale
globnpps<-globnpp*0.1#http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/GeoTIFF/MOD17A3/readme.txt 
globnpps[globnpps==6553.5]<-NA#Set NAs
plot(globnpps)#g/m2/yr  #Note missing data
globnpp_laea<-projectRaster(globnpps,herbivore_dataset3)
globnpp_m<-mask(crop(globnpp_laea,herbivore_dataset3[[1]]),herbivore_dataset3[[1]])
plot(globnpp_m)


#Stack all up
AllVars<-stack(distBiomeBound,vegcov,pcstack,meantraitstack,dietstacksum,sdtraitstack,fdivstack,bioclimlaea_m,globnpp_m)
names(AllVars)
names(AllVars)[1]<-'DistanceBiomeBoundary'
names(AllVars)[51]<-'NPP'
AllVars$bio10<-AllVars$bio10/10
writeRaster(AllVars,'AnalysisVars/',by.layer=T,suffix=names(AllVars),overwrite=T)

AVlist<-list.files('AnalaysisVars/')
AllVars<-stack('AnalaysisVars/')
names(AllVars)

pairs(AllVars[[c(1:6,12)]])
pairs(AllVars[[c(1:5,7:11)]])

AllVars_ext<-extract(AllVars,1:ncell(AllVars),df=T)
cellsXY<-xyFromCell(AllVars,1:ncell(AllVars))
AllVars_ex<-cbind(AllVars_ext,cellsXY)
AllVars_ex1<-merge(AllVars_ex,m1,by.x='ID',by.y='cell')
Av1<-AllVars_ex1[!is.na(AllVars_ex1$PC1),]
Av1$DistanceBiomeBoundary_km<-Av1$DistanceBiomeBoundary/1000
Av1$TreeShrubCover<-Av1$treeCover+Av1$ShrubCover
Av1$Biome<-Av1$DistanceBiomeBoundary_km
Av1$Biome[Av1$DistanceBiomeBoundary_km<0]<-'Boreal'
Av1$Biome[Av1$DistanceBiomeBoundary_km>0]<-'Tundra'
names(Av1)<-sub('.x','',names(Av1))
#write.csv(Av1,'AnalysisVariables.csv')

Av1<-read.csv('AnalysisVariables.csv',header=T)
plot(Av1$DistanceBiomeBoundary_km,Av1$ShrubCover,cex=0.5,pch=16)
plot(Av1$DistanceBiomeBoundary_km,Av1$Body.Mass,cex=0.5,pch=16)

plot(Av1$DistanceBiomeBoundary_km,Av1$Shrubs_sum,cex=0.5,pch=16)
plot(Av1$ShrubCover,Av1$Shrubs_sum,cex=0.2,pch=16)
plot(Av1$ShrubCover,Av1$Shrubs_mean,cex=0.2,pch=16)


#Segmented regression
#Tree cover
lmT<-lm(treeCover~DistanceBiomeBoundary_km,data=Av1)
st1<-segmented(lmT,seg.Z = ~DistanceBiomeBoundary_km,psi=c(-1000,1000))
plot(Av1$DistanceBiomeBoundary_km,Av1$treeCover,pch=16,cex=0.1)
plot(st1,add=T,col=2)

#Tree and shrub cover
lmST<-lm(TreeShrubCover~DistanceBiomeBoundary_km,data=Av1)
sst1<-segmented(lmST,seg.Z = ~DistanceBiomeBoundary_km,psi=c(-1000,1000))
plot(Av1$DistanceBiomeBoundary_km,Av1$TreeShrubCover,pch=16,cex=0.1)
plot(sst1,add=T,col=2)

#Shrubs in diets Biome
lmS<-lm(Shrubs_sum~DistanceBiomeBoundary_km,data=Av1)
s1<-segmented.lm(lmS,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(s1)
plot(Av1$DistanceBiomeBoundary_km,Av1$Shrubs_sum,pch=16,cex=0.1)
plot.segmented(s1,add=T,col=2)

#PC1 Biome
lmPC1<-lm(PC1~DistanceBiomeBoundary_km,data=Av1)
sPC1<-segmented.lm(lmPC1,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sPC1)
plot(Av1$DistanceBiomeBoundary_km,Av1$PC1,pch=16,cex=0.1)
plot.segmented(sPC1,add=T,col=2)

#PC2 Biome
lmPC2<-lm(PC2~DistanceBiomeBoundary_km,data=Av1)
sPC2<-segmented.lm(lmPC2,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sPC2)
plot(Av1$DistanceBiomeBoundary_km,Av1$PC2,pch=16,cex=0.1)
plot.segmented(sPC2,add=T,col=2)

#Trait coords
traitcoords<-read.csv('FunctionalClassification/Trait_coords.csv')

#PC fig
tiff('Figures/PCA.tiff',height=6,width=12,units='in',res=150)
par(mfrow=c(1,2))
cols=c('darkgreen','orange3')
pcdf<-Av1[!is.na(Av1$PC1)&!is.na(Av1$DistanceBiomeBoundary_km),]
plot(famdcoord[,2:3],type='n')
text(famdcoord[,2:3],famdcoord[,1],cex=0.4)
with(Av1[Av1$DistanceBiomeBoundary_km<0,],points(PC1,PC2,cex=0.1,col=cols[1]))#Boreal
with(Av1[Av1$DistanceBiomeBoundary_km>0,],points(PC1,PC2,cex=0.1,col=cols[2]))#Tundra
#bhull<-chull(pcdf[pcdf$DistanceBiomeBoundary_km<0,5],pcdf[pcdf$DistanceBiomeBoundary_km<0,6])#Boreal
#bhull<-c(bhull,bhull[1])
#thull<-chull(pcdf[pcdf$DistanceBiomeBoundary_km>0,5],pcdf[pcdf$DistanceBiomeBoundary_km>0,6])#TUndra
#thull<-c(thull,thull[1])
#lines(pcdf[pcdf$DistanceBiomeBoundary_km<0,][bhull,5:6],col=cols[1])
#lines(pcdf[pcdf$DistanceBiomeBoundary_km>0,][thull,5:6],col=cols[2])
legend('topl',pch=16,col=cols,paste(levels(as.factor(Av1$Biome))))
#with(pcdf,dataEllipse(x=PC1,y=PC2,groups=as.factor(Biome),levels=0.95))
#text(traitcoords[,2:3],as.character(traitcoords[,1]),cex=0.8,col=2)

plot(famdcoord[,2:3],type='n',xlim=c(-4,2))
text(traitcoords[,2:3],as.character(traitcoords[,1]),cex=0.4)

dev.off()


plot(famdcoord[,3:4],type='n')
text(famdcoord[,3:4],famdcoord[,1],cex=0.4)
with(Av1[Av1$DistanceBiomeBoundary_km<0,],points(PC1,PC2,cex=0.1,col=2))
with(Av1[Av1$DistanceBiomeBoundary_km>0,],points(PC1,PC2,cex=0.1,col=3))

#Variables with high loading in res.famd axis 2 expected to show correlation with biome boundary distance
sort(res.famd$var$coord[,2])

#Segmented regression plots
#Biome dist
ndDistBiome<-data.frame(DistanceBiomeBoundary_km=seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T),max(Av1$DistanceBiomeBoundary_km,na.rm=T),length.out = 100))

par(mfrow=c(3,3))
#with(Av1,plot(DistanceBiomeBoundary_km,log(Body.Mass),cex=0.1,col=2))
lmBody.Mass<-lm(log(Body.Mass)~DistanceBiomeBoundary_km,data=Av1)
summary(lmBody.Mass)
sBody.Mass<-segmented.lm(lmBody.Mass,seg.Z=~DistanceBiomeBoundary_km)
summary(sBody.Mass)
plot(Av1$DistanceBiomeBoundary_km,log(Av1$Body.Mass),pch=16,cex=0.1)
plot.segmented(sBody.Mass,add=T,col=2)
pseg_sBody.Mass<-data.frame(predict(sBody.Mass,ndDistBiome,interval='confidence',level=0.95))
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sBody.Mass$upr,col=2)
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sBody.Mass$lwr,col=2)
abline(v=0)

#with(Av1,plot(DistanceBiomeBoundary_km,(Litter.Clutch.size),cex=0.1,col=2))
lmLitter.Clutch.size<-lm((Litter.Clutch.size)~DistanceBiomeBoundary_km,data=Av1)
summary(lmLitter.Clutch.size)
sLitter.Clutch.size<-segmented.lm(lmLitter.Clutch.size,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sLitter.Clutch.size)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Litter.Clutch.size),pch=16,cex=0.1)
plot.segmented(sLitter.Clutch.size,add=T,col=2)
pseg_sLitter.Clutch.size<-data.frame(predict(sLitter.Clutch.size,ndDistBiome,interval='confidence',level=0.95))
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sLitter.Clutch.size$upr,col=2)
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sLitter.Clutch.size$lwr,col=2)


#with(Av1,plot(DistanceBiomeBoundary_km,(Shrubs_mean),cex=0.1,col=2))
lmShrubs_mean<-lm((Shrubs_mean)~DistanceBiomeBoundary_km,data=Av1)
summary(lmShrubs_mean)
sShrubs_mean<-segmented.lm(lmShrubs_mean,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sShrubs_mean)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Shrubs_mean),pch=16,cex=0.1)
plot.segmented(sShrubs_mean,add=T,col=2)
pseg_sShrubsm<-data.frame(predict(sShrubs_mean,ndDistBiome,interval='confidence',level=0.95))
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubsm$upr,col=2)
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubsm$lwr,col=2)


#with(Av1,plot(DistanceBiomeBoundary_km,(Shrubs_sum),cex=0.1,col=2))
lmShrubs_sum<-lm((Shrubs_sum)~DistanceBiomeBoundary_km,data=Av1)
summary(lmShrubs_sum)
sShrubs_sum<-segmented.lm(lmShrubs_sum,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sShrubs_sum)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Shrubs_sum),pch=16,cex=0.1)
plot.segmented(sShrubs_sum,add=T,col=2)
pseg_sShrubs<-data.frame(predict(sShrubs_sum,ndDistBiome,interval='confidence',level=0.95))
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubs$upr,col=2)
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubs$lwr,col=2)

with(Av1,plot(DistanceBiomeBoundary_km,Population_dynamics_cyclic,cex=0.1,col=2))
lmPopDyn_cyclic<-lm(Population_dynamics_cyclic~DistanceBiomeBoundary_km,data=Av1)
summary(lmPopDyn_cyclic)
sPopDyn_cyclic<-segmented.lm(lmPopDyn_cyclic,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sPopDyn_cyclic)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Population_dynamics_cyclic),pch=16,cex=0.1)
plot.segmented(sPopDyn_cyclic,add=T,col=2)
pseg_sShrubs<-data.frame(predict(sPopDyn_cyclic,ndDistBiome,interval='confidence',level=0.95))
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubs$upr,col=2)
lines(ndDistBiome$DistanceBiomeBoundary_km,pseg_sShrubs$lwr,col=2)



##
#FD as sd
par(mfrow=c(1,2))
with(Av1,plot(DistanceBiomeBoundary_km,Body.Mass,cex=0.1,col=2))
with(Av1,plot(DistanceBiomeBoundary_km,Body.Mass_sd,cex=0.1,col=2))

par(mfrow=c(1,2))
with(Av1,plot(DistanceBiomeBoundary_km,PC1,cex=0.1,col=2))
with(Av1,plot(DistanceBiomeBoundary_km,PC1_sd,cex=0.1,col=2))
par(mfrow=c(1,2))
with(Av1,plot(DistanceBiomeBoundary_km,PC2,cex=0.1,col=2))
with(Av1,plot(DistanceBiomeBoundary_km,PC2_sd,cex=0.1,col=2))


par(mfrow=c(1,3))
with(Av1,plot(DistanceBiomeBoundary_km,Shrubs_sum,cex=0.1,col=2))
with(Av1,plot(DistanceBiomeBoundary_km,Shrubs_mean,cex=0.1,col=2))
with(Av1,plot(DistanceBiomeBoundary_km,Shrubs_sd,cex=0.1,col=2))

par(mfrow=c(1,3))
with(Av1,plot(ShrubCover,Shrubs_sum,cex=0.1,col=2))
with(Av1,plot(ShrubCover,Shrubs_mean,cex=0.1,col=2))
with(Av1,plot(ShrubCover,Shrubs_sd,cex=0.1,col=2))

par(mfrow=c(1,3))
with(Av1,plot(TreeShrubCover,Shrubs_sum,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Shrubs_mean,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Shrubs_sd,cex=0.1,col=2))

par(mfrow=c(1,3))
with(Av1,plot(TreeShrubCover,Graminoid_sum,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Graminoids_mean,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Graminoids_sd,cex=0.1,col=2))

par(mfrow=c(1,3))
with(Av1,plot(TreeShrubCover,Forbs_sum,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Forbs_mean,cex=0.1,col=2))
with(Av1,plot(TreeShrubCover,Forbs_sd,cex=0.1,col=2))

par(mfrow=c(2,3))
with(Av1,plot(DistanceBiomeBoundary_km,FRic,cex=0.1))
with(Av1,plot(DistanceBiomeBoundary_km,FDiv,cex=0.1))
with(Av1,plot(DistanceBiomeBoundary_km,FDis,cex=0.1))
with(Av1,plot(DistanceBiomeBoundary_km,FEve,cex=0.1))
with(Av1,plot(DistanceBiomeBoundary_km,RaoQ,cex=0.1))

par(mfrow=c(2,3))
with(Av1,plot(TreeShrubCover,FRic,cex=0.1))
with(Av1,plot(TreeShrubCover,FDiv,cex=0.1))
with(Av1,plot(TreeShrubCover,FDis,cex=0.1))
with(Av1,plot(TreeShrubCover,FEve,cex=0.1))
with(Av1,plot(TreeShrubCover,RaoQ,cex=0.1))

par(mfrow=c(2,3))
with(Av1,plot(treeCover,FRic,cex=0.1))
with(Av1,plot(treeCover,FDiv,cex=0.1))
with(Av1,plot(treeCover,FDis,cex=0.1))
with(Av1,plot(treeCover,FEve,cex=0.1))
with(Av1,plot(treeCover,RaoQ,cex=0.1))

par(mfrow=c(2,3))
with(Av1,plot(ShrubCover,FRic,cex=0.1))
with(Av1,plot(ShrubCover,FDiv,cex=0.1))
with(Av1,plot(ShrubCover,FDis,cex=0.1))
with(Av1,plot(ShrubCover,FEve,cex=0.1))
with(Av1,plot(ShrubCover,RaoQ,cex=0.1))


par(mfrow=c(2,3))
with(Av1,plot(bio10,FRic,cex=0.1))
with(Av1,plot(bio10,FDiv,cex=0.1))
with(Av1,plot(bio10,FDis,cex=0.1))
with(Av1,plot(bio10,FEve,cex=0.1))
with(Av1,plot(bio10,RaoQ,cex=0.1))

par(mfrow=c(2,3))
with(Av1,plot(bio12,FRic,cex=0.1))
with(Av1,plot(bio12,FDiv,cex=0.1))
with(Av1,plot(bio12,FDis,cex=0.1))
with(Av1,plot(bio12,FEve,cex=0.1))
with(Av1,plot(bio12,RaoQ,cex=0.1))

par(mfrow=c(2,3))
with(Av1,plot(NPP,FRic,cex=0.1))
with(Av1,plot(NPP,FDiv,cex=0.1))
with(Av1,plot(NPP,FDis,cex=0.1))
with(Av1,plot(NPP,FEve,cex=0.1))
with(Av1,plot(NPP,RaoQ,cex=0.1))

with(Av1,plot(bio10,NPP))
with(Av1,plot(bio10,DistanceBiomeBoundary_km))

with(Av1,plot(treeCover,Use_of_vegetation_ground_vegetation,cex=0.1))
with(Av1,lines(loess.smooth(treeCover,Use_of_vegetation_ground_vegetation),col=2))
with(Av1,plot(treeCover,Use_of_vegetation_ground_vegetation,cex=0.1))
with(Av1,lines(loess.smooth(treeCover,Use_of_vegetation_ground_vegetation),col=2))



# Biome boundary figure ---------------------------------------------------
#Biome boundary figure
tiff('Figures/BiomeBoundary.tif',width=6, height=10,units='in',res=250)
{
cols=c('darkgreen','orange3')
par(mfcol=c(5,2))
par(oma=c(5,1,1,1))
par(mar=c(1,5,0,1))

with(Av1,plot(DistanceBiomeBoundary_km,FRic,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",las=1))
axis(1,labels=F)
legend('topr',pch=16,col=cols,legend=paste(levels(as.factor(Av1$Biome))))
lmFRic<-lm(FRic~DistanceBiomeBoundary_km,data=Av1)
summary(lmFRic)
davies.test(lmFRic,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sFRic<-segmented.lm(lmFRic,seg.Z=~DistanceBiomeBoundary_km)
summary(sFRic)
abline(v=sFRic$psi[2],lty=1,col=2)
abline(v=sFRic$psi[2]+sFRic$psi[3],lty=2,col=2)
abline(v=sFRic$psi[2]-sFRic$psi[3],lty=2,col=2)
plot.segmented(sFRic,add=T,col=2,conf.level = 0.99,shade=T)

with(Av1,plot(DistanceBiomeBoundary_km,FDiv,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",las=1))
axis(1,labels=F)
lmFDiv<-lm(FDiv~DistanceBiomeBoundary_km,data=Av1)
summary(lmFDiv)
davies.test(lmFDiv,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sFDiv<-segmented.lm(lmFDiv,seg.Z=~DistanceBiomeBoundary_km)
summary(sFDiv)
abline(v=sFDiv$psi[2],lty=1,col=2)
abline(v=sFDiv$psi[2]+sFDiv$psi[3],lty=2,col=2)
abline(v=sFDiv$psi[2]-sFDiv$psi[3],lty=2,col=2)
plot.segmented(sFDiv,add=T,col=2,conf.level = 0.99,shade=T)

with(Av1,plot(DistanceBiomeBoundary_km,FDis,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",las=1))
axis(1,labels=F)
lmFDis<-lm(FDis~DistanceBiomeBoundary_km,data=Av1)
summary(lmFDis)
davies.test(lmFDis,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sFDis<-segmented.lm(lmFDis,seg.Z=~DistanceBiomeBoundary_km)
summary(sFDis)
abline(v=sFDis$psi[2],lty=1,col=2)
abline(v=sFDis$psi[2]+sFDis$psi[3],lty=2,col=2)
abline(v=sFDis$psi[2]-sFDis$psi[3],lty=2,col=2)
plot.segmented(sFDis,add=T,col=2,conf.level = 0.99,shade=T)

with(Av1,plot(DistanceBiomeBoundary_km,FEve,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",las=1))
axis(1,labels=F)
lmFEve<-lm(FEve~DistanceBiomeBoundary_km,data=Av1)
summary(lmFEve)
davies.test(lmFEve,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
#sFEve<-segmented.lm(lmFEve,seg.Z=~DistanceBiomeBoundary_km)
#summary(sFEve)
#abline(v=sFEve$psi[2],lty=1,col=2)
#abline(v=sFEve$psi[2]+sFEve$psi[3],lty=2,col=2)
#abline(v=sFEve$psi[2]-sFEve$psi[3],lty=2,col=2)
#plot.segmented(sFEve,add=T,col=2,conf.level = 0.99,shade=T)
xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
ypt <- data.frame(predict(lmFEve, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = "confidence") )
lines(xdat,ypt$fit)
lines(ypt$lwr ~ xdat, lwd = 1.5, lty = 2,col=1)
lines(ypt$upr ~ xdat, lwd = 1.5, lty = 2,col=1)

#with(Av1,plot(DistanceBiomeBoundary_km,PC1,cex=0.1))
with(Av1,plot(DistanceBiomeBoundary_km,PC2,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",las=1))
axis(1,labels=T)
mtext('Distance to biome boundary (km)',side=1,line=3,cex=0.8)
mtext('Arctic tundra',side=1,line=-1,adj=0.9,cex=0.8)
mtext('Boreal forest',side=1,line=-1,adj=0.1,cex=0.8)
lmPC2<-lm(PC2~DistanceBiomeBoundary_km,data=Av1)
summary(lmPC2)
davies.test(lmPC2,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
#sPC2<-segmented.lm(lmPC2,seg.Z=~DistanceBiomeBoundary_km)
#summary(sPC2)
#abline(v=sPC2$psi[2],lty=1,col=2)
#abline(v=sPC2$psi[2]+sPC2$psi[3],lty=2,col=2)
#abline(v=sPC2$psi[2]-sPC2$psi[3],lty=2,col=2)
#plot.segmented(sPC2,add=T,col=2,conf.level = 0.99,shade=T)
xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
ypt <- data.frame(predict(lmPC2, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = "confidence",level=0.99) )
lines(xdat,ypt$fit)
lines(ypt$lwr ~ xdat, lwd = 1.5, lty = 2,col=1)
lines(ypt$upr ~ xdat, lwd = 1.5, lty = 2,col=1)

with(Av1,plot(DistanceBiomeBoundary_km,log(body_mass),cex=0.1,col=cols[as.factor(Biome)],xaxt='n',ylab="",las=1))
axis(1,labels=F)
title(ylab='Body mass \n(log kg)')
lmbody_mass<-lm(log(body_mass)~DistanceBiomeBoundary_km,data=Av1)
summary(lmbody_mass)
davies.test(lmbody_mass,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
#sbody_mass<-segmented.lm(lmbody_mass,seg.Z=~DistanceBiomeBoundary_km)
#summary(sbody_mass)
#abline(v=sbody_mass$psi[2],lty=1,col=2)
#abline(v=sbody_mass$psi[2]+sbody_mass$psi[3],lty=2,col=2)
#abline(v=sbody_mass$psi[2]-sbody_mass$psi[3],lty=2,col=2)
#plot.segmented(sbody_mass,add=T,col=2,conf.level = 0.99,shade=T)
xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
ypt <- data.frame(predict(lmbody_mass, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = 'confidence') )
lines(xdat,ypt$fit)
lines(xdat,ypt$lwr, lwd = 1.5, lty = 2,col=1)
lines(xdat,ypt$upr, lwd = 1.5, lty = 2,col=1)

with(Av1,plot(DistanceBiomeBoundary_km,Shrubs_mean,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",ylab="",las=1))
axis(1,labels=F)
title(ylab='Woody plant diet \n(index)')
lmShrubs_mean<-lm((Shrubs_mean)~DistanceBiomeBoundary_km,data=Av1)
summary(lmShrubs_mean)
davies.test(lmShrubs_mean,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sShrubs_mean<-segmented.lm(lmShrubs_mean,seg.Z=~DistanceBiomeBoundary_km)
summary(sShrubs_mean)
abline(v=sShrubs_mean$psi[2],lty=1,col=2)
abline(v=sShrubs_mean$psi[2]+sShrubs_mean$psi[3],lty=2,col=2)
abline(v=sShrubs_mean$psi[2]-sShrubs_mean$psi[3],lty=2,col=2)
plot.segmented(sShrubs_mean,add=T,col=2,conf.level = 0.99,shade=T)
#xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
#ypt <- data.frame(predict(lmShrubs_mean, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = 'confidence') )
#lines(xdat,ypt$fit)
#lines(xdat,ypt$lwr, lwd = 1.5, lty = 2,col=1)
#lines(xdat,ypt$upr, lwd = 1.5, lty = 2,col=1)

with(Av1,plot(DistanceBiomeBoundary_km,Use_of_vegetation_ground_vegetation,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",ylab="",las=1))
axis(1,labels=F)
title(ylab='Ground foraging \n(index)')
lmUse_of_vegetation_ground_vegetation<-lm((Use_of_vegetation_ground_vegetation)~DistanceBiomeBoundary_km,data=Av1)
summary(lmUse_of_vegetation_ground_vegetation)
davies.test(lmUse_of_vegetation_ground_vegetation,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sUse_of_vegetation_ground_vegetation<-segmented.lm(lmUse_of_vegetation_ground_vegetation,seg.Z=~DistanceBiomeBoundary_km)
summary(sUse_of_vegetation_ground_vegetation)
abline(v=sUse_of_vegetation_ground_vegetation$psi[2],lty=1,col=2)
abline(v=sUse_of_vegetation_ground_vegetation$psi[2]+sUse_of_vegetation_ground_vegetation$psi[3],lty=2,col=2)
abline(v=sUse_of_vegetation_ground_vegetation$psi[2]-sUse_of_vegetation_ground_vegetation$psi[3],lty=2,col=2)
plot.segmented(sUse_of_vegetation_ground_vegetation,add=T,col=2,conf.level = 0.99,shade=T)
#xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
#ypt <- data.frame(predict(lmUse_of_vegetation_ground_vegetation, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = 'confidence') )
#lines(xdat,ypt$fit)
#lines(xdat,ypt$lwr, lwd = 1.5, lty = 2,col=1)
#lines(xdat,ypt$upr, lwd = 1.5, lty = 2,col=1)

with(Av1,plot(DistanceBiomeBoundary_km,Population_dynamics_cyclic,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",ylab="",las=1))
axis(1,labels=F)
title(ylab='Cyclic population \n dynamics (index)')
lmPopulation_dynamics_cyclic<-lm((Population_dynamics_cyclic)~DistanceBiomeBoundary_km,data=Av1)
summary(lmPopulation_dynamics_cyclic)
davies.test(lmPopulation_dynamics_cyclic,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sPopulation_dynamics_cyclic<-segmented.lm(lmPopulation_dynamics_cyclic,seg.Z=~DistanceBiomeBoundary_km)
summary(sPopulation_dynamics_cyclic)
abline(v=sPopulation_dynamics_cyclic$psi[2],lty=1,col=2)
abline(v=sPopulation_dynamics_cyclic$psi[2]+sPopulation_dynamics_cyclic$psi[3],lty=2,col=2)
abline(v=sPopulation_dynamics_cyclic$psi[2]-sPopulation_dynamics_cyclic$psi[3],lty=2,col=2)
plot.segmented(sPopulation_dynamics_cyclic,add=T,col=2,conf.level = 0.99,shade=T)
#xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
#ypt <- data.frame(predict(lmPopulation_dynamics_cyclic, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = 'confidence') )
#lines(xdat,ypt$fit)
#lines(xdat,ypt$lwr, lwd = 1.5, lty = 2,col=1)
#lines(xdat,ypt$upr, lwd = 1.5, lty = 2,col=1)

# with(Av1,plot(DistanceBiomeBoundary_km,Belowground_feeding_Belowground_feeding_grubbing,cex=0.1,col=cols[as.factor(Biome)]))
# lmBelowground_feeding_Belowground_feeding_grubbing<-lm((Belowground_feeding_Belowground_feeding_grubbing)~DistanceBiomeBoundary_km,data=Av1)
# summary(lmBelowground_feeding_Belowground_feeding_grubbing)
# davies.test(lmBelowground_feeding_Belowground_feeding_grubbing,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
# sBelowground_feeding_Belowground_feeding_grubbing<-segmented.lm(lmBelowground_feeding_Belowground_feeding_grubbing,seg.Z=~DistanceBiomeBoundary_km)
# summary(sBelowground_feeding_Belowground_feeding_grubbing)
# abline(v=sBelowground_feeding_Belowground_feeding_grubbing$psi[2],lty=1,col=2)
# abline(v=sBelowground_feeding_Belowground_feeding_grubbing$psi[2]+sBelowground_feeding_Belowground_feeding_grubbing$psi[3],lty=2,col=2)
# abline(v=sBelowground_feeding_Belowground_feeding_grubbing$psi[2]-sBelowground_feeding_Belowground_feeding_grubbing$psi[3],lty=2,col=2)
# plot.segmented(sBelowground_feeding_Belowground_feeding_grubbing,add=T,col=2,conf.level = 0.99,shade=T)
# #xdat <- seq(min(Av1$DistanceBiomeBoundary_km,na.rm=T), max(Av1$DistanceBiomeBoundary_km,na.rm=T), length.out = 100)
# #ypt <- data.frame(predict(lmBelowground_feeding_Belowground_feeding_grubbing, newdata = data.frame(DistanceBiomeBoundary_km=xdat), interval = 'confidence') )
# #lines(xdat,ypt$fit)
# #lines(xdat,ypt$lwr, lwd = 1.5, lty = 2,col=1)
# #lines(xdat,ypt$upr, lwd = 1.5, lty = 2,col=1)
# 
# 
# with(Av1,plot(DistanceBiomeBoundary_km,Belowground_feeding_Belowground_feeding_burrowing,cex=0.1,col=cols[as.factor(Biome)]))
# lmBelowground_feeding_Belowground_feeding_burrowing<-lm((Belowground_feeding_Belowground_feeding_burrowing)~DistanceBiomeBoundary_km,data=Av1)
# summary(lmBelowground_feeding_Belowground_feeding_burrowing)
# davies.test(lmBelowground_feeding_Belowground_feeding_burrowing,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
# sBelowground_feeding_Belowground_feeding_burrowing<-segmented.lm(lmBelowground_feeding_Belowground_feeding_burrowing,seg.Z=~DistanceBiomeBoundary_km)
# summary(sBelowground_feeding_Belowground_feeding_burrowing)
# abline(v=sBelowground_feeding_Belowground_feeding_burrowing$psi[2],lty=1,col=2)
# abline(v=sBelowground_feeding_Belowground_feeding_burrowing$psi[2]+sBelowground_feeding_Belowground_feeding_burrowing$psi[3],lty=2,col=2)
# abline(v=sBelowground_feeding_Belowground_feeding_burrowing$psi[2]-sBelowground_feeding_Belowground_feeding_burrowing$psi[3],lty=2,col=2)
# plot.segmented(sBelowground_feeding_Belowground_feeding_burrowing,add=T,col=2,conf.level = 0.99,shade=T)

with(Av1,plot(DistanceBiomeBoundary_km,Belowground_feeding_Belowground_feeding_none,cex=0.1,col=cols[as.factor(Biome)],xaxt='n',xlab="",ylab="",las=1))
axis(1,labels=T)
title(ylab='Aboveground feeding \n(index)')
mtext("Distance to biome boundary (km)",side=1,line=3,cex=0.8)
mtext('Arctic tundra',side=1,line=-1,adj=0.9,cex=0.8)
mtext('Boreal forest',side=1,line=-1,adj=0.1,cex=0.8)
lmBelowground_feeding_Belowground_feeding_none<-lm((Belowground_feeding_Belowground_feeding_none)~DistanceBiomeBoundary_km,data=Av1)
summary(lmBelowground_feeding_Belowground_feeding_none)
davies.test(lmBelowground_feeding_Belowground_feeding_none,seg.Z=~DistanceBiomeBoundary_km,alternative = 'less')
sBelowground_feeding_Belowground_feeding_none<-segmented.lm(lmBelowground_feeding_Belowground_feeding_none,seg.Z=~DistanceBiomeBoundary_km)
summary(sBelowground_feeding_Belowground_feeding_none)
abline(v=sBelowground_feeding_Belowground_feeding_none$psi[2],lty=1,col=2)
abline(v=sBelowground_feeding_Belowground_feeding_none$psi[2]+sBelowground_feeding_Belowground_feeding_none$psi[3],lty=2,col=2)
abline(v=sBelowground_feeding_Belowground_feeding_none$psi[2]-sBelowground_feeding_Belowground_feeding_none$psi[3],lty=2,col=2)
plot.segmented(sBelowground_feeding_Belowground_feeding_none,add=T,col=2,conf.level = 0.99,shade=T)
}
dev.off()

# AbioticBiotic FD Figure -------------------------------------------------


par(mfcol=c(4,3))
par(mar=c(1,1,1,1))
par(oma=c(5,5,1,1))
with(Av1,plot(TreeShrubCover,FRic,cex=0.1,las=1,xlab=F,ylab='Functional richness',xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
mtext('Functional richness',side=2,line=3,cex=0.8)
axis(1,tick=T,labels = F)
legend('bottomr',pch=16,col=cols,legend=paste(levels(as.factor(Av1$Biome))))
with(Av1,plot(TreeShrubCover,FDiv,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = T)
mtext('Functional divergence',side=2,line=3,cex=0.8)
with(Av1,plot(TreeShrubCover,FDis,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = T)
mtext('Functional dissimilarity',side=2,line=3,cex=0.8)
with(Av1,plot(TreeShrubCover,FEve,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(2,las=1,labels = T)
mtext('Functional evenness',side=2,line=3,cex=0.8)
axis(1,tick=T,labels = T)
mtext('Woody plant cover (%)',side=1,line=3,cex=0.8)
with(Av1,plot(NPP,FRic,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(NPP,FDiv,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(NPP,FDis,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(NPP,FEve,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = T)
axis(2,las=1,labels = F)
mtext(expression('Net primary productivity (g m'^-2*' yr'^-1*')'),side=1,line=3,cex=0.8)
with(Av1,plot(bio10,FRic,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(bio10,FDiv,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(bio10,FDis,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = F)
axis(2,las=1,labels = F)
with(Av1,plot(bio10,FEve,cex=0.1,las=1,xlab=F,ylab=F,xaxt='n',yaxt='n',col=cols[as.factor(Biome)]))
axis(1,tick=T,labels = T)
axis(2,las=1,labels = F)
mtext(expression('Mean summer temperature ('~degree*'C)'),side=1,line=3,cex=0.8)



# SEMs --------------------------------------------------------------------


#Piecewise SEMs
#Temp, NPP and woody plant cover as drivers of trait diversity in herbivore communities
library(piecewiseSEM)
library(nlme)

#Drop NA rows
semdf<-Av1[!is.na(Av1$FRic)&!is.na(Av1$NPP)&!is.na(Av1$bio1)&!is.na(Av1$ShrubCover),]
dim(semdf)
dim(Av1)

sem_standdf<-semdf[,c(3:52,55:56)]
sem_standdf<-data.frame(scale(sem_standdf,scale = T,center=T))
sem_standdf<-cbind(sem_standdf,semdf[,53:54])#Add XY coords


model <- psem(lm(FRic ~ bio10 + NPP + TreeShrubCover, semdf), lm(TreeShrubCover ~ bio10, semdf), lm(NPP ~ bio10, semdf))
summary(model)
coefs(model, standardize = "scale")
coefs(model, standardize = "range")
a$Estimate
library(gstat)
va<-variogram(FRic ~ bio10 + NPP + TreeShrubCover,data= semdf,locations= ~x+y)
plot(va)


model_stand <- psem(lm(FRic ~ bio10 + NPP + TreeShrubCover, sem_standdf), lm(TreeShrubCover ~ bio10, sem_standdf), lm(NPP ~ bio10, sem_standdf))
summary(model_stand)
vas<-variogram(FRic ~ bio10 + NPP + TreeShrubCover,data= sem_standdf,locations= ~x+y)
plot(vas,smooth=T)


#GLS
#Fit GLS for all component models
gls_1<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= semdf,correlation=corLin(form=~x+y, nugget=T), method="ML")
gls_2<-gls(NPP ~ bio10, data= semdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_3<-gls(TreeShrubCover ~ bio10, data= semdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_list<-list(gls_1,gls_2,gls_3)
psemlist1<-as.psem(gls_list)
summary(psemlist1, .progressBar = T)

#Spatial autocor?
semdf$resids<-residuals(gls1)
vario <- Variogram(gls_1, form = ~x + y, resType = "pearson")
plot(vario,smooth=T)

#Standardised

#Find best correlation structure
gls_stand_nc<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf, method="ML")
gls_stand_exp<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_lin<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corLin(form=~x+y, nugget=T), method="ML")
gls_stand_gau<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corGaus(form=~x+y, nugget=T), method="ML")
gls_stand_rat<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corRatio(form=~x+y, nugget=T), method="ML")

AIC(gls_stand_nc,gls_stand_exp,gls_stand_gau,gls_stand_rat)#corLin convergence error so ommitted
#Exp has lowest AIC
varioNull<-Variogram(gls_stand_nc,form=~x+y, resType="pearson")
varioExp<-Variogram(gls_stand_exp,form=~x+y, resType="pearson")
par(mfrow=c(1,2))
plot(varioNull,smooth=T)
plot(varioExp,smooth=T,add=T)

#FRic
gls_stand_FRic<-gls(FRic ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_FRic<-list(gls_stand_FRic,gls_stand_2,gls_stand_3)
psem_stand_list1_FRic<-as.psem(gls_stand_list_FRic)
ps_sum_FRic<-summary(psem_stand_list1_FRic, .progressBar = T)
ps_sum_FRic
vario_FRic<-Variogram(gls_stand_FRic,form= ~x +y,resType = "pearson")
plot(vario_FRic,smooth=T)

#FEve
gls_stand_FEve<-gls(FEve ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_FEve<-list(gls_stand_FEve,gls_stand_2,gls_stand_3)
psem_stand_list1_FEve<-as.psem(gls_stand_list_FEve)
ps_sum_FEve<-summary(psem_stand_list1_FEve, .progressBar = T)
ps_sum_FEve
vario_FEve<-Variogram(gls_stand_FEve,form= ~x +y,resType = "pearson")
plot(vario_FEve,smooth=T)

#FDis
gls_stand_FDis<-gls(FDis ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_FDis<-list(gls_stand_FDis,gls_stand_2,gls_stand_3)
psem_stand_list1_FDis<-as.psem(gls_stand_list_FDis)
ps_sum_FDis<-summary(psem_stand_list1_FDis, .progressBar = T)
ps_sum_FDis
vario_FDis<-Variogram(gls_stand_FDis,form= ~x +y,resType = "pearson")
plot(vario_FDis,smooth=T)

#FDiv
gls_stand_FDiv<-gls(FDiv ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_FDiv<-list(gls_stand_FDiv,gls_stand_2,gls_stand_3)
psem_stand_list1_FDiv<-as.psem(gls_stand_list_FDiv)
ps_sum_FDiv<-summary(psem_stand_list1_FDiv, .progressBar = T)
ps_sum_FDiv
vario_FDiv<-Variogram(gls_stand_FDiv,form= ~x +y,resType = "pearson")
plot(vario_FDiv,smooth=T)

#RaoQ
gls_stand_RaoQ<-gls(RaoQ ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_RaoQ<-list(gls_stand_RaoQ,gls_stand_2,gls_stand_3)
psem_stand_list1_RaoQ<-as.psem(gls_stand_list_RaoQ)
ps_sum_RaoQ<-summary(psem_stand_list1_RaoQ, .progressBar = T)
ps_sum_RaoQ
vario_RaoQ<-Variogram(gls_stand_RaoQ,form= ~x +y,resType = "pearson")
plot(vario_RaoQ,smooth=T)

#PC1
gls_stand_PC1<-gls(PC1 ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_PC1<-list(gls_stand_PC1,gls_stand_2,gls_stand_3)
psem_stand_list1_PC1<-as.psem(gls_stand_list_PC1)
ps_sum_PC1<-summary(psem_stand_list1_PC1, .progressBar = T)
ps_sum_PC1
vario_PC1<-Variogram(gls_stand_PC1,form= ~x +y,resType = "pearson")
plot(vario_PC1,smooth=T)

#PC2
gls_stand_PC2<-gls(PC2 ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_PC2<-list(gls_stand_PC2,gls_stand_2,gls_stand_3)
psem_stand_list1_PC2<-as.psem(gls_stand_list_PC2)
ps_sum_PC2<-summary(psem_stand_list1_PC2, .progressBar = T)
ps_sum_PC2
vario_PC2<-Variogram(gls_stand_PC2,form= ~x +y,resType = "pearson")
plot(vario_PC2,smooth=T)

#Body.Mass
gls_stand_Body.Mass<-gls(Body.Mass ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Body.Mass<-list(gls_stand_Body.Mass,gls_stand_2,gls_stand_3)
psem_stand_list1_Body.Mass<-as.psem(gls_stand_list_Body.Mass)
ps_sum_Body.Mass<-summary(psem_stand_list1_Body.Mass, .progressBar = T)
ps_sum_Body.Mass
vario_Body.Mass<-Variogram(gls_stand_Body.Mass,form= ~x +y,resType = "pearson")
plot(vario_Body.Mass,smooth=T)

#Litter.Clutch.size
gls_stand_Litter.Clutch.size<-gls(Litter.Clutch.size ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Litter.Clutch.size<-list(gls_stand_Litter.Clutch.size,gls_stand_2,gls_stand_3)
psem_stand_list1_Litter.Clutch.size<-as.psem(gls_stand_list_Litter.Clutch.size)
ps_sum_Litter.Clutch.size<-summary(psem_stand_list1_Litter.Clutch.size, .progressBar = T)
ps_sum_Litter.Clutch.size
vario_Litter.Clutch.size<-Variogram(gls_stand_Litter.Clutch.size,form= ~x +y,resType = "pearson")
plot(vario_Litter.Clutch.size,smooth=T)


#Shrubs_mean
gls_stand_Shrubs_mean<-gls(Shrubs_mean ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Shrubs_mean<-list(gls_stand_Shrubs_mean,gls_stand_2,gls_stand_3)
psem_stand_list1_Shrubs_mean<-as.psem(gls_stand_list_Shrubs_mean)
ps_sum_Shrubs_mean<-summary(psem_stand_list1_Shrubs_mean, .progressBar = T)
ps_sum_Shrubs_mean
vario_Shrubs_mean<-Variogram(gls_stand_Shrubs_mean,form= ~x +y,resType = "pearson")
plot(vario_Shrubs_mean,smooth=T)

#Graminoids_mean
gls_stand_Graminoids_mean<-gls(Graminoids_mean ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Graminoids_mean<-list(gls_stand_Graminoids_mean,gls_stand_2,gls_stand_3)
psem_stand_list1_Graminoids_mean<-as.psem(gls_stand_list_Graminoids_mean)
ps_sum_Graminoids_mean<-summary(psem_stand_list1_Graminoids_mean, .progressBar = T)
ps_sum_Graminoids_mean
vario_Graminoids_mean<-Variogram(gls_stand_Graminoids_mean,form= ~x +y,resType = "pearson")
plot(vario_Graminoids_mean,smooth=T)

#Forbs_mean
gls_stand_Forbs_mean<-gls(Forbs_mean ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Forbs_mean<-list(gls_stand_Forbs_mean,gls_stand_2,gls_stand_3)
psem_stand_list1_Forbs_mean<-as.psem(gls_stand_list_Forbs_mean)
ps_sum_Forbs_mean<-summary(psem_stand_list1_Forbs_mean, .progressBar = T)
ps_sum_Forbs_mean
vario_Forbs_mean<-Variogram(gls_stand_Forbs_mean,form= ~x +y,resType = "pearson")
plot(vario_Forbs_mean,smooth=T)

#Bryophytes_mean
gls_stand_Bryophytes_mean<-gls(Bryophytes_mean ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Bryophytes_mean<-list(gls_stand_Bryophytes_mean,gls_stand_2,gls_stand_3)
psem_stand_list1_Bryophytes_mean<-as.psem(gls_stand_list_Bryophytes_mean)
ps_sum_Bryophytes_mean<-summary(psem_stand_list1_Bryophytes_mean, .progressBar = T)
ps_sum_Bryophytes_mean
vario_Bryophytes_mean<-Variogram(gls_stand_Bryophytes_mean,form= ~x +y,resType = "pearson")
plot(vario_Bryophytes_mean,smooth=T)

#Lichens_mean
gls_stand_Lichens_mean<-gls(Lichens_mean ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Lichens_mean<-list(gls_stand_Lichens_mean,gls_stand_2,gls_stand_3)
psem_stand_list1_Lichens_mean<-as.psem(gls_stand_list_Lichens_mean)
ps_sum_Lichens_mean<-summary(psem_stand_list1_Lichens_mean, .progressBar = T)
ps_sum_Lichens_mean
vario_Lichens_mean<-Variogram(gls_stand_Lichens_mean,form= ~x +y,resType = "pearson")
plot(vario_Lichens_mean,smooth=T)

#Shrubs_sum
gls_stand_Shrubs_sum<-gls(Shrubs_sum ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Shrubs_sum<-list(gls_stand_Shrubs_sum,gls_stand_2,gls_stand_3)
psem_stand_list1_Shrubs_sum<-as.psem(gls_stand_list_Shrubs_sum)
ps_sum_Shrubs_sum<-summary(psem_stand_list1_Shrubs_sum, .progressBar = T)
ps_sum_Shrubs_sum
vario_Shrubs_sum<-Variogram(gls_stand_Shrubs_sum,form= ~x +y,resType = "pearson")
plot(vario_Shrubs_sum,smooth=T)

#Graminoid_sum
gls_stand_Graminoid_sum<-gls(Graminoid_sum ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Graminoid_sum<-list(gls_stand_Graminoid_sum,gls_stand_2,gls_stand_3)
psem_stand_list1_Graminoid_sum<-as.psem(gls_stand_list_Graminoid_sum)
ps_sum_Graminoid_sum<-summary(psem_stand_list1_Graminoid_sum, .progressBar = T)
ps_sum_Graminoid_sum
vario_Graminoid_sum<-Variogram(gls_stand_Graminoid_sum,form= ~x +y,resType = "pearson")
plot(vario_Graminoid_sum,smooth=T)

#Forbs_sum
gls_stand_Forbs_sum<-gls(Forbs_sum ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Forbs_sum<-list(gls_stand_Forbs_sum,gls_stand_2,gls_stand_3)
psem_stand_list1_Forbs_sum<-as.psem(gls_stand_list_Forbs_sum)
ps_sum_Forbs_sum<-summary(psem_stand_list1_Forbs_sum, .progressBar = T)
ps_sum_Forbs_sum
vario_Forbs_sum<-Variogram(gls_stand_Forbs_sum,form= ~x +y,resType = "pearson")
plot(vario_Forbs_sum,smooth=T)

#Bryophytes_sum
gls_stand_Bryophytes_sum<-gls(Moss_sum ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Bryophytes_sum<-list(gls_stand_Bryophytes_sum,gls_stand_2,gls_stand_3)
psem_stand_list1_Bryophytes_sum<-as.psem(gls_stand_list_Bryophytes_sum)
ps_sum_Bryophytes_sum<-summary(psem_stand_list1_Bryophytes_sum, .progressBar = T)
ps_sum_Bryophytes_sum
vario_Bryophytes_sum<-Variogram(gls_stand_Bryophytes_sum,form= ~x +y,resType = "pearson")
plot(vario_Bryophytes_sum,smooth=T)

#Lichens_sum
gls_stand_Lichens_sum<-gls(Lichens_sum ~ bio10 + NPP + TreeShrubCover, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_2<-gls(NPP ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
#gls_stand_3<-gls(TreeShrubCover ~ bio10, data= sem_standdf,correlation=corExp(form=~x+y, nugget=T), method="ML")
gls_stand_list_Lichens_sum<-list(gls_stand_Lichens_sum,gls_stand_2,gls_stand_3)
psem_stand_list1_Lichens_sum<-as.psem(gls_stand_list_Lichens_sum)
ps_sum_Lichens_sum<-summary(psem_stand_list1_Lichens_sum, .progressBar = T)
ps_sum_Lichens_sum
vario_Lichens_sum<-Variogram(gls_stand_Lichens_sum,form= ~x +y,resType = "pearson")
plot(vario_Lichens_sum,smooth=T)