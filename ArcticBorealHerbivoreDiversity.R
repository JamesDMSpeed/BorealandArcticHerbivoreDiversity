#Analysing boreal and arctic herbivore diversity

#Load packages
require(raster) #Spatial
require(rgdal)#Spatial
require(sp)#Spatial
require(rasterVis)#Spatial
require(picante)#Diversity analysis
require(vegan)#Diversity analysis
require(ape)#View dendrograms

# Set up ------------------------------------------------------------------

#Polar projection
polarproj<-'+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 '

#Country outlines
noreco_countries <- c('NO', 'SE', 'FI','RU','CA','IS','GL','SJ','MN','JP') 
noreco_shp1 <- do.call("bind", lapply(noreco_countries, function(x)  raster::getData('GADM', country=x, level=0)))
#Only Alaska from USA
us<-getData('GADM',country='US',level=1)
alaska<-us[us$NAME_1=='Alaska',]
#Bind
noreco_shp<-bind(noreco_shp1,alaska)
#Project
noreco_shppp<-spTransform(noreco_shp,CRS=crs(polarproj))

plot(noreco_shppp)


#Biomes
ecoregions<-readOGR('Biomes','wwf_terr_ecos')
tundraboreal<-ecoregions[ecoregions$BIOME==6|ecoregions$BIOME==11,]
northernecosystems<-crop(tundraboreal,extent(-180,180,40,90))#Remove southern hemisphere
#Project
northernecosystemspp<-spTransform(northernecosystems,polarproj)
plot(northernecosystemspp)



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


#
#
#Here we need to remove unused species, and add those with missing data!!!
#Not complete!!

spplist<-read.csv('TraitSpreadSheetJan2019.csv')
spplist$Binomial

#Match names
#Replace dot with space to match with raster stack names
spplist1<-spplist$Binomial
spplist1<-sub(' ','.',spplist1)
spplist1%in%names(herbivore_dataset)

#These species are ok
spplist1[which(spplist1%in%names(herbivore_dataset))]
spplist1[which(spplist1%in%names(herbivore_dataset)==F)]#These species are not in the spatial data 
names(herbivore_dataset)[which(names(herbivore_dataset)%in%spplist1==F)]#These species are in the spatial data but not spp list

#Several of these are synonyms
#Anser==Chen
#Anas==Mareca
#Urocitellus.parryii==Spermophilus.parryii

#Remove livestock from spatial data
livestocklist<-c('Bos.taurus','Capra.aegagrus','Ovis.aries')
herbivore_dataset<-herbivore_dataset[[which(names(herbivore_dataset)%in%livestocklist==F)]]

#Two species in spatial dataset not on trait list
#Aythya.collaris
#Ursus.arctos  

#Species on trait list, not in spatial data - check these
#Aix.galericulata
#Allactaga.major
#Allactaga.sibirica
#Allocricetulus.eversmanni
#Anas.platyrhynchos
#Anas.rubripes"            
#Cervus.canadensis #Should seperate this from Cervus elaphys in distribiton data #Also join Alces alces and A. americanus
# Glaucomys.sabrinus   
# Melanitta.americana 
# Melanitta.deglandi  
# Melanitta.stejnegeri
# Micromys.minutus
# Microtus.hyperboreus  
# Napaeozapus.insignis      
# Phasianus.colchicus     
# Phodopus.campbelli
# Phodopus.sungorus
# Pteromys.volans
# Sciurus.carolinensis
# Sciurus.vulgaris
# Spatula.querquedula
# Tamias.striatus


#Not complete!! - check the above



# Species diversity analysis ----------------------------------------------


#Diversity analysis
#Calculate species richness 
sr<-calc(herbivore_dataset,fun=sum,na.rm=T)
#Cells with sr>0
speciesrichness<-sr
speciesrichness[sr==0]<-NA
levelplot(speciesrichness,par.settings=YlOrRdTheme,margin=F)#+
  layer(sp.polygons(noreco_shppp))#Need to get better outline map at somepoint
  

#

# Species based cluster analysis ------------------------------------------

#Community data
herbcomdata<-getValues(herbivore_dataset)
herbcomdata1<-herbcomdata[rowSums(herbcomdata,na.rm=T)>1,]
herbcomdata1[is.na(herbcomdata1)]<-0
#Distance matrix
dist1<-vegdist(herbcomdata1)#Bray curtis
#Hierarchical clustering - Ward.D method
hc1<-hclust(dist1,method='ward.D')
plot(hc1)

#Cut tree (check optimal number of clusters first)
clusters<-cutree(hc1,k=8)

#Make a raster layer for cluster and populate it 
speciesclusts<-speciesrichness
speciesclusts[speciesrichness>1]<-clusters

mycol<-brewer.pal(8,'Dark2')
levelplot(speciesclusts,margin=F,scales=list(draw=FALSE),col.regions=mycol,colorkey=list(at=seq(0.5,8.5,by=1)))+
  layer(sp.polygons(noreco_shppp))#Need to get better outline map at somepoint

#hc1$labels<-rep('llllllllllllllll',times=length(hc1$height)) #Bodge a coloured bar as a label
#plot(as.phylo(hc1),tip.color=mycol[clusters],cex=0.5,no.margin=T)
