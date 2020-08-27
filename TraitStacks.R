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
speciesmat<-speciesdat[,order(names(speciesdata))]#Order alphabeticaaly to match traits
speciesCommunities<-speciesmat[rowSums(speciesmat[,2:ncol(speciesmat)],na.rm=T)>0,]#Remove empty rows #Omit ID column
speciescoms<-speciesCommunities[,colSums(speciesCommunities,na.rm=T)>0]#Remove empty columns
s2<-speciescoms[,2:ncol(speciescoms)]
traitdf1<-traitdf[rownames(traitdf)%in%names(s2),]#Remove empty species from trait data too

funcdiv<-dbFD(x=traitdf1[,c(1:2,5,11:15)],a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F)
funcdiv<-dbFD(x=traitdf1[,c(3:4)],a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F)
funcdiv<-dbFD(x=traitdf1,a=speciescoms[,2:ncol(speciescoms)],corr='lingoes',w.abun=F,stand.x=F)


#Fill raster with FD indices
dfFD<-data.frame(ID=speciescoms$ID,FDiv=funcdiv$FDiv,FRic=funcdiv$FRic,FEve=funcdiv$FEve,FDis=funcdiv$FDis,RaoQ=funcdiv$RaoQ)#Dataframe with Fdiv indices and cell numbers
df2<-data.frame(cell=seq(1,ncell(herbivore_dataset3[[1]]),1))#Dataframe with all (inc NA) cells
m1<-merge(df2,dfFD,by.x='cell',by.y='ID',all.x=T)#Merge together
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



#Stack all up
AllVars<-stack(distBiomeBound,vegcov,pcstack,meantraitstack,dietstacksum,sdtraitstack,fdivstack)
writeRaster(AllVars,'AnalysisVars/',by.layer=T,suffix=names(AllVars),overwrite=T)
names(AllVars)[1]<-'DistanceBiomeBoundary'
names(AllVars)

pairs(AllVars[[c(1:6,12)]])
pairs(AllVars[[c(1:5,7:11)]])

AllVars_ex<-extract(AllVars,1:ncell(AllVars),df=T)
Av1<-AllVars_ex[!is.na(AllVars_ex$PC1),]
write.csv(Av1,'AnalysisVariables.csv')
Av1$DistanceBiomeBoundary_km<-Av1$DistanceBiomeBoundary/1000
Av1$TreeShrubCover<-Av1$treeCover+Av1$ShrubCover

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


#PC fig
plot(famdcoord[,2:3],type='n')
text(famdcoord[,2:3],famdcoord[,1],cex=0.4)
with(Av1[Av1$DistanceBiomeBoundary_km<0,],points(PC1,PC2,cex=0.1,col=2))#Boreal
with(Av1[Av1$DistanceBiomeBoundary_km>0,],points(PC1,PC2,cex=0.1,col=4))#Tundra

plot(famdcoord[,3:4],type='n')
text(famdcoord[,3:4],famdcoord[,1],cex=0.4)
with(Av1[Av1$DistanceBiomeBoundary_km<0,],points(PC1,PC2,cex=0.1,col=2))
with(Av1[Av1$DistanceBiomeBoundary_km>0,],points(PC1,PC2,cex=0.1,col=3))

#Variables with high loading in res.famd axis 2 expected to show correlation with biome boundary distance
sort(res.famd$var$coord[,2])

with(Av1,plot(DistanceBiomeBoundary_km,log(Body.Mass),cex=0.1,col=2))
lmBody.Mass<-lm((Body.Mass)~DistanceBiomeBoundary_km,data=Av1)
summary(lmBody.Mass)
sBody.Mass<-segmented.lm(lmBody.Mass,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sBody.Mass)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Body.Mass),pch=16,cex=0.1)
plot.segmented(sBody.Mass,add=T,col=2)

with(Av1,plot(DistanceBiomeBoundary_km,(Litter.Clutch.size),cex=0.1,col=2))
lmLitter.Clutch.size<-lm((Litter.Clutch.size)~DistanceBiomeBoundary_km,data=Av1)
summary(lmLitter.Clutch.size)
sLitter.Clutch.size<-segmented.lm(lmLitter.Clutch.size,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sLitter.Clutch.size)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Litter.Clutch.size),pch=16,cex=0.1)
plot.segmented(sLitter.Clutch.size,add=T,col=2)

with(Av1,plot(DistanceBiomeBoundary_km,(Shrubs_mean),cex=0.1,col=2))
lmShrubs_mean<-lm((Shrubs_mean)~DistanceBiomeBoundary_km,data=Av1)
summary(lmShrubs_mean)
sShrubs_mean<-segmented.lm(lmShrubs_mean,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sShrubs_mean)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Shrubs_mean),pch=16,cex=0.1)
plot.segmented(sShrubs_mean,add=T,col=2)


with(Av1,plot(DistanceBiomeBoundary_km,(Shrubs_sum),cex=0.1,col=2))
lmShrubs_sum<-lm((Shrubs_sum)~DistanceBiomeBoundary_km,data=Av1)
summary(lmShrubs_sum)
sShrubs_sum<-segmented.lm(lmShrubs_sum,seg.Z=~DistanceBiomeBoundary_km,psi=c(-1000,1000))
summary(sShrubs_sum)
plot(Av1$DistanceBiomeBoundary_km,(Av1$Shrubs_sum),pch=16,cex=0.1)
plot.segmented(sShrubs_sum,add=T,col=2)


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


#Temperature & NDVI etc 