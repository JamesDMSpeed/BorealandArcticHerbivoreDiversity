###---------------------------------get rasters from IUCN data-------------------------###

###---------------------------------------Attempt 1------------------------------------###
#get(alcesdstf)
#rastermlayer<-raster(alcesdstf)#create master layer from previous
#getValues(rastermlayer)#check extent of new master layer
#ccanashp<-(readOGR('data_0.shp'))#works if set WD where file is located (ONLY)
#ccanashpp<-spTransform(ccanashp,CRS=crs(polarproj))#Change projection

#attempt to transform new shapefile into raster 
#ccanaraster <- rasterize(ccanashpp, rastermlayer)
#plot(ccanaraster)

#plot(herbivore_dataset[[144]])
#names(herbivore_dataset)
#plot(ccanashpp)
#Get flat outer eges with this
###------------------------------------Attempt 1 FAILED---------------------------------###

###########################################################################################
###--------------------------------Attempt 2 (Fix flat edges)---------------------------###

###---------------------------ONLY HAS TO BE DONE ONCE PER SESSION----------------------###

rml<-calc(herbivore_dataset,fun=sum,na.rm=T)#calculate species richness to get maximum extent map with proper dimensions
plot(rml)#Check that SR map has full extension

 
                                                       
rml2<-rml          #Create an extra modifiable object and transform into basic presence/null absence map
                   #to cut new species into shape
rml2[rml==0]<-NA   #Change all values =0 to null values
rml2[rml>0]<-1     #Change all the values >0 to 1 values 
plot(rml2)         #to check proper p/n coverage

###----------------------------------EVERY TIME----------------------------------------###

#Call and the shape file as an object and change projection
shp<-(readOGR('data_0.shp'))#works if set WD where file is located (ONLY)
shpp<-spTransform(shp,CRS=crs(polarproj))#Change projection

# Transform new shapefile into raster 
raster <- rasterize(shpp, rml2)
plot(raster)

# set the cells associated with the shapfile to the specified value if different value than 1
raster[!is.na(raster)] <- 2

# merge the new raster with the mask raster with mosaic to have 3 different values.
raster2<-mosaic(raster,rml2,fun=mean)
plot(raster2)

#Discriminate the values to maintain only overlapping of the desired species
raster2[raster2==1]<-NA      #Change all the values =1 to NA values
raster2[raster2>1.5]<-NA     #Change all the values >1.5 to NA values
raster2[raster2==1.5]<-1     #Change all the valuse =1.5 to 1 values
plot(raster2)

#Extract the resulting map!

writeRaster(raster2,filename='_Tamias.striatus.tif')


