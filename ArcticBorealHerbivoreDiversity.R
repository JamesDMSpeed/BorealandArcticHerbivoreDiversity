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

#------Fixed Different Extent Issue-------------------------------------


#
#Here we need to remove unused species, and add those with missing data!!!
#Not complete!!

spplist<-read.csv('TraitTableFeb2019.csv')
spplist$Binomial

#Match names
#Replace dot with space to match with raster stack names
spplist1<-spplist$Binomial
spplist1<-sub(' ','.',spplist1)
spplist1%in%names(herbivore_dataset3)

#These species are ok
spplist1[which(spplist1%in%names(herbivore_dataset3))]
spplist1[which(spplist1%in%names(herbivore_dataset3)==F)]#These species are not in the spatial data 
names(herbivore_dataset3)[which(names(herbivore_dataset3)%in%spplist1==F)]#These species are in the spatial data but not spp list

#Several of these are synonyms
#Anser==Chen
#Anas==Mareca
#Urocitellus.parryii==Spermophilus.parryii

#Remove livestock & Irrelevant Herbivores from spatial data
#livestocklist<-c('Bos.taurus','Capra.aegagrus','Ovis.aries')
#herbivore_dataset3<-herbivore_dataset3[[which(names(herbivore_dataset3)%in%livestocklist==F)]]



#Two species in spatial dataset not on trait list
#Aythya.collaris (DONE)
#Ursus.arctos (DONE)

#Species on trait list, not in spatial data - check these
#Aix.galericulata (DONE)
#Allactaga.major(DONE)
#Allactaga.sibirica(DONE)
#Allocricetulus.eversmanni(DONE)
#Anas.platyrhynchos(DONE)
#Anas.rubripes"  (DONE)          
#Cervus.canadensis (DONE)
#Cervus.elaphus (DONE)
#Glaucomys.sabrinus (DONE)  
#Melanitta.americana (DONE) 
#Melanitta.deglandi (DONE) 
#Melanitta.stejnegeri(DONE)
#Micromys.minutus(DONE)
#Microtus.hyperboreus (DONE) 
#Napaeozapus.insignis   (DONE)   
#Phasianus.colchicus   (DONE)  
#Phodopus.campbelli (DONE)  
#Phodopus.sungorus (DONE)
#Pteromys.volans (DONE)
#Sciurus.carolinensis (DONE)
#Sciurus.vulgaris(DONE)
#Spatula.querquedula(DONE)
#Tamias.striatus (DONE)



#Not complete!! - check the above

#Simple biome map
simplebiome<-rasterize(northernecosystemspp,herbivore_dataset3,field='BIOME')

#Convert boreal biome to points
borealpoints<-rasterToPoints(simplebiome,fun=function(x)x==6,spatial=T)
arcticpoints<-rasterToPoints(simplebiome,fun=function(x)x==11,spatial=T)
plot(simplebiome)
points(borealpoints,cex=0.1,pch=16)
points(arcticpoints,cex=0.1,pch=16,col=2)

#Distance to boreal and arctic
distancetoBoreal<-distanceFromPoints(simplebiome,borealpoints)
distancetoArctic<-distanceFromPoints(simplebiome,arcticpoints)
plot(distancetoBoreal)
plot(distancetoArctic)

#Combine
distanceToBiomeBoundary<-distancetoBoreal
distanceToBiomeBoundary[distanceToBiomeBoundary==0]<-0-(distancetoArctic[distanceToBiomeBoundary==0])
distanceToBiomeBoundary<-mask(distanceToBiomeBoundary,simplebiome)
plot(distanceToBiomeBoundary)
writeRaster(distanceToBiomeBoundary,'Biomes/DistnaceBiomeBoundary')

pBB<-levelplot(distanceToBiomeBoundary,margin=F,scales=list(draw=F))+
  layer(sp.polygons(northernecosystemspp))
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
diverge0(pBB,'RdBu')



#Conservation of Arctic Flora and Fauna Working Group (2010) CAFF Map No.53 - Boundaries of the geographic area covered by the Arctic Biodiversity Assessment.
arczones<-readOGR('Biomes/ABA-Boundaries','Arctic_Zones')
arczones_laea<-spTransform(arczones,polarproj)
#Subarctic lower boundary
subarcbound<-arczones_laea[arczones_laea@data$Zone=='Sub arctic',]

# Species diversity analysis ----------------------------------------------


#Diversity analysis
#Calculate species richness 
sr<-calc(herbivore_dataset3,fun=sum,na.rm=T)
#Cells with sr>0
speciesrichness<-sr
speciesrichness[sr==0]<-NA
levelplot(speciesrichness,par.settings=YlOrRdTheme,margin=F)#+
  layer(sp.polygons(noreco_shppp))#Need to get better outline map at somepoint
  

#Better map
northerneco_studyregion<-rasterize(northernecosystemspp,herbivore_dataset3,field='BIOME')
ext<-as.vector(extent(projectRaster(northerneco_studyregion,crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')))
boundaries <- map('worldHires', fill=TRUE,
                 xlim=ext[1:2], ylim=ext[3:4],
                 plot=FALSE)
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                             proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
bPolslaea<-spTransform(bPols,crs(speciesrichness))

levelplot(speciesrichness,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F))+
  layer(sp.polygons(bPolslaea))+
  #layer(sp.polygons(arczones_laea,lty=2))+
  layer(sp.polygons(subarcbound))


####-----------------------------Basic Raster Descriptive stats-------------####

cellStats(speciesrichness, stat='mean', na.rm=TRUE, asSample=TRUE)
cellStats(speciesrichness, stat='max', na.rm=TRUE, asSample=TRUE)
cellStats(speciesrichness, stat='min', na.rm=TRUE, asSample=TRUE)

srfrequency <- data.frame(freq(speciesrichness))
sum(srfrequency$count)

#### to obtain number of cells, subtract NA cell count from total sum of count

# #Phylogeny --------------------------------------------------------------

#Read in phylogeny
phylogeny<-read.tree('Phylogeny/BestTree.newick')
plot(phylogeny)
phylogeny$tip.label

#Change format of tip.labels to match species data
phylogeny$tip.label<-gsub('_','.',phylogeny$tip.label)

#Match synonym names #Complete
phylogeny$tip.label[phylogeny$tip.label=="Anas.querquedula"] <- "Spatula.querquedula"
phylogeny$tip.label[phylogeny$tip.label=="Anser.cygnoides"]<-  "Anser.cygnoid"
#phylogeny$tip.label[phylogeny$tip.label=="Chen.rossii"]<- "Anser.rossii" 
#phylogeny$tip.label[phylogeny$tip.label=="Spermophilus.parryii"]<-"Urocitellus.parryii"


#Convert raster stack to community dataframe
communitydata<- getValues(herbivore_dataset3)
#Replace NA with 0
communitydata[is.na(communitydata)]<-0

#Use picante to trim community and phylogenetic data
phydata<-match.phylo.comm(phylogeny,communitydata)

###
##Check through the dropped species carefully and fix any that are errors.###
###Checked, No errors, extra species in phylogeny should be dropped.
#
##Locate spatial data for those that are missing
###

#Calculate phylogenetic diversity
phydiv<-pd(phydata$comm,phydata$phy,include.root=T)

#Rasterize this
phydivraster<-raster(speciesrichness)
phydivraster<-setValues(phydivraster,phydiv$PD)
phydivraster<-mask(phydivraster,speciesrichness,maskvalue=NA)

levelplot(phydivraster,par.settings=YlOrRdTheme,margin=F)+
  layer(sp.polygons(bPolslaea))+
  layer(sp.polygons(arczones_laea,lty=2))



#Stack together - each as a proportion of the total species richness or phylogeney branch length
diversitystack<-stack(speciesrichness/nlayers(herbivore_dataset3),phydivraster/sum(phylogeny$edge.length))
names(diversitystack)<-c('Species richness','Phylogenetic diversity')

levelplot(diversitystack,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F))+
  layer(sp.polygons(bPolslaea))

####-------------------Functional Diversity------------------------------####

#Read in Functional Dendrogram
ftt<-read.tree('functional_tree.nwk')
#plot(ftt)
ftt$tip.label

#Change format of tip.labels to match species data
ftt$tip.label<-gsub('_','.',ftt$tip.label)

#Convert raster stack to community dataframe
communitydata<- getValues(herbivore_dataset3)
#Replace NA with 0
communitydata[is.na(communitydata)]<-0

#Use vegan to trim community and functional data
fttdata<-match.phylo.comm(ftt,communitydata)

#Calculate phylogenetic diversity
fundiv<-pd(fttdata$comm,fttdata$phy,include.root=T)

#Rasterize this
fundivraster<-raster(speciesrichness)
fundivraster<-setValues(fundivraster,fundiv$PD)
fundivraster<-mask(fundivraster,speciesrichness,maskvalue=NA)

levelplot(fundivraster,par.settings=YlOrRdTheme,margin=F)+
  layer(sp.polygons(bPolslaea))

#Stack together - each as a proportion of the total species richness or phylogeney branch length
diversitystack2<-stack(speciesrichness/nlayers(herbivore_dataset3),fundivraster/sum(ftt$edge.length))
names(diversitystack2)<-c('Species richness','Functional diversity')

levelplot(diversitystack2,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F))+
  layer(sp.polygons(bPolslaea))



#If I do 3 then its a proportion of all three operation? Results look odd
#Stack together - each as a proportion of the total species richness or phylogeney branch length
diversitystack3<-stack(speciesrichness/nlayers(herbivore_dataset3),phydivraster/sum(phylogeny$edge.length),fundivraster/sum(ftt$edge.length))
names(diversitystack3)<-c('Species richness','Phylogenetic diversity','Functional diversity')

levelplot(diversitystack3,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F))+
  layer(sp.polygons(bPolslaea))+
  layer(sp.polygons(arczones_laea,lty=2))


####--------------------------------------------------------------------###

# Phylogenetic Diversity pairplots ----------------------------------------
diversitydata<-data.frame(getValues(diversitystack))

diversitydata$biomeval<-getValues(simplebiome)

with(diversitydata,plot(Species.richness,Phylogenetic.diversity,type='n')) #Very low PD in sites with low SR - only birds in these?
with(diversitydata[diversitydata$biomeval==6,],points(Species.richness,Phylogenetic.diversity,pch=16,cex=0.5,col='green3'))
with(diversitydata[diversitydata$biomeval==11,],points(Species.richness,Phylogenetic.diversity,pch=16,cex=0.5,col='tan4'))
legend('bottomr',pch=16,col=c('green3','tan4'),c('Forest','Tundra'))

#Linear model
lmPDSR<-with(diversitydata,lm(Phylogenetic.diversity~Species.richness))
summary(lmPDSR)
newdat<-data.frame(Species.richness=seq(min(diversitydata$Species.richness,na.rm=T),max(diversitydata$Species.richness,na.rm=T),length.out=100))
newdat$predpd_lin<-predict(lmPDSR,newdata=newdat)
lines(newdat$Species.richness,newdat$predpd_lin) #Linear is best

#Log model
logmPDSR<-with(diversitydata,lm(Phylogenetic.diversity~log(Species.richness)))
summary(logmPDSR)
newdat$predpd_log<-predict(logmPDSR,newdata=newdat)
lines(newdat$Species.richness,newdat$predpd_log) #Linear is best

#Phylogenetic diversity in relation to species richness
pdresiduals<-residuals(lmPDSR)
pdresidualsraster<-raster(speciesrichness)
pdresidualsraster[!is.na(speciesrichness)]<-pdresiduals
pdresidualsraster<-mask(pdresidualsraster,speciesrichness,maskvalue=NA)

#Make a function to plot diverging colour scales around 0 
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


#Plot to show where PD is greater/lower than expected given species richness
pdp<-levelplot(pdresidualsraster,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F),main='Phylogenetic diversity')+
  layer(sp.polygons(bPolslaea,lwd=0.5,col=grey(0.5)))+
  layer(sp.polygons(arczones_laea,lty=2))
diverge0(pdp,'RdBu') #Higher PD than expected in Quebec. Lower than expected in Siberia. Also W coast of Norway

meltdat<-melt(diversitydata,measure.vars=c('Species.richness','Phylogenetic.diversity'))

require(ggplot2)
vio1<-ggplot(data=meltdat[!is.na(meltdat$biomeval),],
             aes(y=value, x=as.factor(biomeval),fill=variable))+
  geom_violin()+
  scale_x_discrete(labels=c("Tundra","Boreal forest"),
                     name="Biome") +
  scale_y_continuous('Diversity (proportion of total)')+
  scale_fill_discrete(name = "Diversity")
vio1 #Tundra has low species richness but high PD

viodat<-data.frame(cbind(pdresiduals,biome=diversitydata$biomeval[!is.na(diversitydata$Species.richness)]))
vio2<-ggplot(data=viodat[!is.na(viodat$biome),],aes(x=as.factor(biome),y=pdresiduals))+geom_violin()+
  scale_x_discrete(labels=c("Tundra","Boreal forest"),
                   name="Biome") +
  scale_y_continuous('Residuals of PD~SR')
vio2

detach(package:ggplot2)#Avoiding conflict with plotting spatial objects 


####----------------Functional Diversity pairplots ----------------------#######
diversitydata2<-data.frame(getValues(diversitystack2))

diversitydata2$biomeval<-getValues(simplebiome)

with(diversitydata2,plot(Species.richness,Functional.diversity,type='n')) #Very low PD in sites with low SR - only birds in these?
with(diversitydata2[diversitydata2$biomeval==6,],points(Species.richness,Functional.diversity,pch=16,cex=0.5,col='green'))
with(diversitydata2[diversitydata2$biomeval==11,],points(Species.richness,Functional.diversity,pch=16,cex=0.5,col='tan4'))
legend('bottomr',pch=16,col=c('green3','tan4'),c('Forest','Tundra'))

#Linear model
lmFDSR<-with(diversitydata2,lm(Functional.diversity~Species.richness))
summary(lmFDSR)
newdat2<-data.frame(Species.richness=seq(min(diversitydata2$Species.richness,na.rm=T),max(diversitydata2$Species.richness,na.rm=T),length.out=100))
newdat2$predfd_lin<-predict(lmFDSR,newdata=newdat2)
lines(newdat2$Species.richness,newdat2$predfd_lin) #Logarithmic is best

#Log model
logmFDSR<-with(diversitydata2,lm(Functional.diversity~log(Species.richness)))
summary(logmFDSR)
newdat2$predfd_log<-predict(logmFDSR,newdata=newdat2)
lines(newdat2$Species.richness,newdat2$predfd_log) #Logarithmic is best

#Functional diversity in relation to species richness
fdresiduals<-residuals(lmFDSR)
fdresidualsraster<-raster(speciesrichness)
fdresidualsraster[!is.na(speciesrichness)]<-fdresiduals
fdresidualsraster<-mask(fdresidualsraster,speciesrichness,maskvalue=NA)

#Make a function to plot diverging colour scales around 0 
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


#Plot to show where FD is greater/lower than expected given species richness
pfp<-levelplot(fdresidualsraster,par.settings=YlOrRdTheme,margin=F,scales=list(draw=F),main='Functional diversity')+
  layer(sp.polygons(bPolslaea,lwd=0.5,col=grey(0.5)))+
  layer(sp.polygons(arczones_laea,lty=2))
diverge0(pfp,'RdBu') #Odd patterns, it apears that The highest ammount of FD relative to SR is in the middle of the distribution

meltdat2<-melt(diversitydata2,measure.vars=c('Species.richness','Functional.diversity'))

require(ggplot2)
vio3<-ggplot(data=meltdat2[!is.na(meltdat2$biomeval),],
             aes(y=value, x=as.factor(biomeval),fill=variable))+
  geom_violin()+
  scale_x_discrete(labels=c("Tundra","Boreal forest"),
                   name="Biome") +
  scale_y_continuous('Diversity (proportion of total)')+
  scale_fill_discrete(name = "Diversity")
vio3 #Tundra has low species richness but high PD

viodat2<-data.frame(cbind(fdresiduals,biome=diversitydata2$biomeval[!is.na(diversitydata2$Species.richness)]))
vio4<-ggplot(data=viodat[!is.na(viodat2$biome),],aes(x=as.factor(biome),y=pdresiduals))+geom_violin()+
  scale_x_discrete(labels=c("Tundra","Boreal forest"),
                   name="Biome") +
  scale_y_continuous('Residuals of FD~SR')
vio4

detach(package:ggplot2)#Avoiding conflict with plotting spatial objects 


######------------------------Latitudinal Comparisons------------------####


#---------------------------------Functional------------------------------#

Functionaldata<-data.frame(getValues(fundivraster))
Functionaldata$biomeval<-getValues(simplebiome)

colnames(Functionaldata) <- c("FD", "biomeval")

with(Functionaldata,plot(Species.richness,Functional.diversity,type='n')) #Very low PD in sites with low SR - only birds in these?
with(diversitydata2[diversitydata2$biomeval==6,],points(Species.richness,Functional.diversity,pch=16,cex=0.5,col='tan4'))
with(diversitydata2[diversitydata2$biomeval==11,],points(Species.richness,Functional.diversity,pch=16,cex=0.5,col='green3'))
legend('bottomr',pch=16,col=c('tan4','green3'),c('Tundra','Forest'))




#######-----------------------------END--------------------------------####

# Cluster analysis ------------------------------------------

###Species

#Community data
herbcomdata<-getValues(herbivore_dataset3)
herbcomdata1<-herbcomdata[rowSums(herbcomdata,na.rm=T)>1,]
herbcomdata1[is.na(herbcomdata1)]<-0
#Distance matrix
dist1<-vegdist(herbcomdata1)#Bray curtis
#Hierarchical clustering - Ward.D method
hc1<-hclust(dist1,method='ward.D')
plot(hc1)

#Cut tree (check optimal number of clusters first)
optclust<-NbClust(dist1,method='ward.D',index='cindex',min.nc=2,max.nc=10)#9
optclust
clusters<-cutree(hc1,k=10)

#Make a raster layer for cluster and populate it 
speciesclusts<-speciesrichness
speciesclusts[speciesrichness>1]<-clusters
speciesclusts[speciesrichness<=1]<-NA

#Ratify the cluster raster to plot as a factor
rat_speciesclusts<-ratify(speciesclusts)
rat <- levels(rat_speciesclusts)[[1]]
rat$clusterID<-rat$ID
levels(rat_speciesclusts)<-rat
#Set colours
mycol<-c(brewer.pal(8,'Dark2'),1,2)
#Plot with country and ecozone outlines
#polyCentroids = gCentroid(northernecosystemspp,byid=T)
levelplot(rat_speciesclusts,att='clusterID',col.regions=mycol,scales=list(draw=F),colorkey=list(title='Cluster'))+
  layer(sp.polygons(bPolslaea,col=grey(0.5)))+
  layer(sp.polygons(arczones_laea,lty=2))
#layer(sp.polygons(northernecosystemspp,lwd=1,lty=2,col='blue'))#+
#layer(sp.text(polyCentroids@coords,northernecosystemspp$ECO_NAME))


#Plot the cluster dendrogram
#hc1$labels<-rep('llllllllllllllll',times=length(hc1$height)) #Bodge a coloured bar as a label
#plot(as.phylo(hc1),tip.color=mycol[clusters],cex=0.5,no.margin=T)
#plot(hc1)


#PhyloSorrensen similarity index

#phylosor1<-function (samp, tree) 
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for phylosor calculation")
  }
  samp <- as.matrix(samp)
  s <- nrow(samp)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(samp)
  colnames(phylodist) <- rownames(samp)
  samp_comb <- matrix(NA, s * (s - 1)/2, ncol(samp))
  colnames(samp_comb) <- colnames(samp)
  i <- 1
  for (l in 1:(s - 1)) {
    print(paste('l=',l))
    for (k in (l + 1):s) {
      samp_comb[i, ] <- samp[l, ] + samp[k, ]
      i <- i + 1
    }
  }
  pdsamp <- pd(samp, tree)
  pdsamp_comb <- pd(samp_comb, tree)
  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdsamp[l, "PD"]
    for (k in (l + 1):s) {
      print(paste('k=',k))
      pdk <- pdsamp[k, "PD"]
      pdcomb <- pdsamp_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] = 2 * pdsharedlk/(pdl + pdk)
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}

#pd_sor1<-phylosor1(phydata$comm[rowSums(phydata$comm)>1,],phydata$phy)#Only for sites with >1sppp
#fd_sor1<-phylosor1(fttdata$comm[rowSums(fttdata$comm)>1,],fttdata$phy)

#pd_sor<-phylosor(phydata$comm[rowSums(phydata$comm)>1,],phydata$phy)#Only for sites with >1sppp
#fd_sor<-phylosor(fttdata$comm[rowSums(fttdata$comm)>1,],fttdata$phy)#Only for sites with >1sppp


#Convert phylogenetic similarity into distances
#pd_sor_dist<-1-pd_sor
#fd_sor_dist<-1-fd_sor

#Write to file
#write.csv(as.matrix(dist1),'SpDistances.csv')
#write.csv(as.matrix(pd_sor_dist),'PdDistances.csv')
#write.csv(as.matrix(fd_sor_dist),'FdDistances.csv')

#These files are available here https://ntnu.box.com/s/9jpqhxpwh51a8lw2ehlactx9ymtk8hkq
#To reimport distance matrices
sd2<-data.matrix(read.csv('SpDistances.csv'))
spdist<-as.dist(sd2[,2:ncol(sd2)])
pd2<-data.matrix(read.csv('PdDistances.csv'))
pddist<-as.dist(pd2[,2:ncol(pd2)])
fd2<-data.matrix(read.csv('FdDistances.csv'))
fddist<-as.dist(fd2[,2:ncol(fd2)])

spclusts<-hclust(spdist,method='ward.D')
pdclusts<-hclust(pddist,method='ward.D')                 
fdclusts<-hclust(fddist,method='ward.D')

par(mfrow=c(1,3))
plot(spclusts,ylim=c(0,1),main='Species')
plot(pdclusts,ylim=c(0,1),main='Phylogenetic')
plot(fdclusts,ylim=c(0,1),main='Functional')

# Optimal number of clusters ----------------------------------------------

require(NbClust)
indexchoice<-'cindex'
methodchoice<-'ward.D'
nbSpp<-NbClust(diss=spdist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)
nbPD<-NbClust(diss=pddist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)
nbFD<-NbClust(diss=fddist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)

nbSpp
nbPD
nbFD
clusters<-c(nbSpp$Best.nc[1],nbPD$Best.nc[1],nbFD$Best.nc[1])
clusters

#Cut tree with optimum number of clusters
spclusts<-hclust(spdist,method=methodchoice)
pdclusts<-hclust(pddist,method=methodchoice)
fdclusts<-hclust(fddist,method=methodchoice)
par(mfrow=c(1,3))
plot(spclusts)
plot(pdclusts)
plot(fdclusts)

spgroups<-cutree(spclusts,nbSpp$Best.nc[1])
pdgroups<-cutree(pdclusts,nbPD$Best.nc[1])
fdgroups<-cutree(fdclusts,nbFD$Best.nc[1])
spgroups<-cutree(spclusts,6)
pdgroups<-cutree(pdclusts,6)
fdgroups<-cutree(fdclusts,6)


#Set up vector to populate
r1<-raster(speciesrichness)
spprich<-rowSums(phydata$comm)
spprich[spprich<2]<-NA

spgroupvalues<-spprich
pdgroupvalues<-spprich
fdgroupvalues<-spprich
spgroupvalues[!is.na(spgroupvalues)]<-spgroups
pdgroupvalues[!is.na(pdgroupvalues)]<-pdgroups
fdgroupvalues[!is.na(fdgroupvalues)]<-fdgroups

sp_sorras<-setValues(r1,spgroupvalues)
pd_sorras<-setValues(r1,pdgroupvalues)
fd_sorras<-setValues(r1,fdgroupvalues)

np<-SpatialPoints(cbind(0,0))
clustpal<-rasterTheme(region=(c(brewer.pal(8,'Dark2'),'black','red','blue','green')))
clustpal<-rasterTheme(region=(c(brewer.pal(6,'Dark2'))))
levelplot(stack(sp_sorras,pd_sorras,fd_sorras),par.settings=clustpal,colorkey=F,names.attr=c('Species','Phylogenetic','Functional'),scales=(list(draw=F)))+
    layer(sp.polygons(arczones_laea,lty=2))+
    layer(sp.points(np,col=1))+
    layer(sp.polygons(bPolslaea))+
    layer(sp.polygons(regionAB))

#Plot showing biome and regional boundaries
r1<-speciesrichness
r1[r1>0]<-1
r1[r1==0]<-NA
levelplot(r1,margin=F,scales=list(draw=F))+
  layer(sp.polygons(arczones_laea,col='red'))+
  layer(sp.polygons(regionAB,col='blue'))+
  layer(sp.polygons(noreco_shppp))+
  layer(sp.points(np))


#Productivity
globnpps[globnpps>5000]<-NA
nppherb<-projectRaster(globnpps,herbivore_dataset3)
plot(nppherb)
