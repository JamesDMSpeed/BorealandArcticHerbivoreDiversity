#Functional classification
library(FactoMineR)#Factorial analysis of mixed data
require(FactoMineR)
require(vegan) # mantel test
require(ape) # export data to biodiverse
require(RColorBrewer)#Colours to match the results figures
require(picante)#testing phylogenetic conservatism
require(missMDA)#Multivariate imputation

#Load and clean data####

#Preliminary trait table - downloaded on 7th Feb 2019
traittable<-read.csv('Functional classification/TraitTableFeb2019.csv',
                     header=T, na.strings=c("","NA"))
View(traittable)
names(traittable)

tail(traittable)
#Remove empty rows
traittable1<-traittable[traittable$Binomial!="" &!is.na(traittable$Binomial),]
## select variables to dataframe that will be used

#Set missing data to NA

traitvar<-c("Order", "Family", "Genus", "Species", "Binomial", "Elton.Plant_Others",
       "body_mass", "gut_type", "group_size_summer", "group_size_winter", 
       "Litter_clutch_size",  "Population_dynamics", "Habitat_type", "Belowground_feeding",
       "Mobility", "Human_managed",  "Diet_type",   "diet_item_forb", "diet_item_graminoid", 
       "diet_item_shrub",   "diet_item_moss"  , "diet_item_lichen",
       "wintering_strategy","Use_of_vegetation" )

Traits<-traittable1[,which(colnames(traittable1) %in% traitvar)]
summary(Traits)





#Fix some discrepancies in trait categories
levels(Traits$Population_dynamics)[which(levels(Traits$Population_dynamics)=="cyclic_noncyclic ")] = "cyclic_noncyclic" 
levels(Traits$Population_dynamics)[which(levels(Traits$Population_dynamics)=="noncylcic")] = "noncyclic"

levels(Traits$Use_of_vegetation)[which(levels(Traits$Use_of_vegetation)=="groun_vegetation")] = "ground_vegetation"
levels(Traits$Use_of_vegetation)[which(levels(Traits$Use_of_vegetation)=="ground_ vegetation")] = "ground_vegetation" 
levels(Traits$Use_of_vegetation)[which(levels(Traits$Use_of_vegetation)=="gound_vegetation")] = "ground_vegetation"
levels(Traits$Use_of_vegetation)[which(levels(Traits$Use_of_vegetation)=="grpund_vegetation")] = "ground_vegetation" 

levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="facuktative_generalist")] = "facultative_generalist"
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="facultative_generalis")] = "facultative_generalist"
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="fscultative_generalist")] = "facultative_generalist" 
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="fscultative_generalist")] = "facultative_generalist" 

levels(Traits$Habitat_type)[which(levels(Traits$Habitat_type)=="limnic ")] = "limnic"

levels(Traits$gut_type)[which(levels(Traits$gut_type)=="hindgut_fermented")] = "hindgut_fermenter"
levels(Traits$gut_type)[which(levels(Traits$gut_type)=="hingut_fermenter")] = "hindgut_fermenter"

levels(Traits$group_size_summer)[which(levels(Traits$group_size_summer)=="small_groups")] = "small_group"
levels(Traits$group_size_summer)[which(levels(Traits$group_size_summer)=="family_groups")] = "family_group"
levels(Traits$group_size_summer)[which(levels(Traits$group_size_summer)=="soliary")] = "solitary"

levels(Traits$group_size_winter)[which(levels(Traits$group_size_winter)=="large group")] = "large_group" 
levels(Traits$group_size_winter)[which(levels(Traits$group_size_winter)=="large_groups")] = "large_group"

levels(Traits$wintering_strategy)[which(levels(Traits$wintering_strategy)=="active above_ground")] = "active_aboveground"
levels(Traits$wintering_strategy)[which(levels(Traits$wintering_strategy)=="active_above_ground")] = "active_aboveground"
levels(Traits$wintering_strategy)[which(levels(Traits$wintering_strategy)=="active_bellowgroud")] = "active_belowground"
levels(Traits$wintering_strategy)[which(levels(Traits$wintering_strategy)=="hibernation")] = "hibernating" 

Traits$Elton.Plant_Others[Traits$Elton.Plant_Others=="Mentiones in HB of B that it's highly vegetarian"]<-NA
Traits$Elton.Plant_Others[Traits$Elton.Plant_Others=="ND"]<-NA
Traits$Elton.Plant_Others<-as.numeric(as.character(Traits$Elton.Plant_Others))

#Imputing missing data####

#Phylogenetics

#Load phylogeny
# arcborphy<-read.tree('Phylogeny/BestTree_Yet2.newick')
# plot(arcborphy)
# arcborphy$tip.label
# #Change format of tip.labels to match trait data
# arcborphy$tip.label<-gsub('_',' ',arcborphy$tip.label)
# 
# rownames(Traits)<-Traits$Binomial
# 
# matched<-match.phylo.data(arcborphy,Traits)
# 
# phyEstimate(matched$phy,matched$data$Litter_clutch_size[!is.na(matched$data$Litter_clutch_size)])
# t2<-Traits[!is.na(Traits$Litter_clutch_size) & !is.na(Traits$body_mass),c(7,11)]
# phyEstimate(arcborphy,t2)

#Multivariate
#missMDA

Traits$diet_item_forb<-as.factor(Traits$diet_item_forb)
Traits$diet_item_graminoid<-as.factor(Traits$diet_item_graminoid)
Traits$diet_item_shrub<-as.factor(Traits$diet_item_shrub)
Traits$diet_item_moss<-as.factor(Traits$diet_item_moss)
Traits$diet_item_lichen<-as.factor(Traits$diet_item_lichen)
str(Traits)

Traits2<-imputeFAMD(Traits)
Traits2$completeObs
Traits2$completeObs$Litter_clutch_size

Traits<-Traits2$completeObs
################# transform variables   ----####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## body mass to log(body mass); otherwise the large values will dominate enormously
Traits$body_mass_log<-log(Traits$body_mass); hist(Traits$body_mass_log)

## group size to numeric with order; gradient from small to large group size 
Traits$group_size_order<-Traits$group_size_summer
str(Traits)
levels(Traits$group_size_order)<- c(levels(Traits$group_size_order), "1", "2", "3")
Traits$group_size_order[Traits$group_size_order == "solitary"] <- "1"
Traits$group_size_order[Traits$group_size_order == "family_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "small_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "large_group"] <- "3"
Traits$group_size_order<-droplevels(Traits$group_size_order)
Traits$group_size_summer_order<-as.numeric(Traits$group_size_order)

Traits$group_size_order<-Traits$group_size_winter
levels(Traits$group_size_order) <- c(levels(Traits$group_size_order), "1", "2", "3")
Traits$group_size_order[Traits$group_size_order == "solitary"] <- "1"
Traits$group_size_order[Traits$group_size_order == "family_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "small_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "large_group"] <- "3"
Traits$group_size_order<-droplevels(Traits$group_size_order)
Traits$group_size_winter_order<-as.numeric(Traits$group_size_order)

## population dynamics; gradient from little temporal variation ot loads of temporal variation
Traits$population_dynamics_order<-Traits$Population_dynamics
levels(Traits$population_dynamics_order) <- c(levels(Traits$population_dynamics_order), "1", "2", "3")
Traits$population_dynamics_order[Traits$population_dynamics_order == "noncyclic"] <- "1"
Traits$population_dynamics_order[Traits$population_dynamics_order == "cyclic_noncyclic"] <- "2"
Traits$population_dynamics_order[Traits$population_dynamics_order == "cyclic"] <- "3"
Traits$population_dynamics_order<-droplevels(Traits$population_dynamics_order)
Traits$population_dynamics_order<-as.numeric(Traits$population_dynamics_order)

## gut type; gradient from "unefficient" to "efficient"
Traits$gut_type_order<-Traits$gut_type
levels(Traits$gut_type_order) <- c(levels(Traits$gut_type_order), "1", "2", "3")
Traits$gut_type_order[Traits$gut_type_order == "undifferentiated"] <- "1"
Traits$gut_type_order[Traits$gut_type_order == "hindgut_fermenter"] <- "2"
Traits$gut_type_order[Traits$gut_type_order == "ruminant"] <- "3"
Traits$gut_type_order<-droplevels(Traits$gut_type_order)
Traits$gut_type_order<-as.numeric(Traits$gut_type_order)


# ## diet type; gradient from "generalist to specialist " towards increasing specialisation
Traits$diet_type_order<-Traits$Diet_type
levels(Traits$diet_type_order) <- c(levels(Traits$diet_type_order), "1", "2", "3", "4")
Traits$diet_type_order[Traits$diet_type_order == "obligatory_generalist"] <- "1"
Traits$diet_type_order[Traits$diet_type_order == "facultative_generalist"] <- "2"
Traits$diet_type_order[Traits$diet_type_order == "facultative_specialist"] <- "3"
Traits$diet_type_order[Traits$diet_type_order == "obligatory_specialist"] <- "4"
Traits$diet_type_order<-droplevels(Traits$diet_type_order)
Traits$diet_type_order<-as.numeric(Traits$diet_type_order)

# ## winter mode; gradient from "certainly not present " towards increasing presence
Traits$wintering_strategy_order<-Traits$wintering_strategy
levels(Traits$wintering_strategy_order)<- c(levels(Traits$wintering_strategy_order), "1", "2", "3", "4") 
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "active_aboveground"] <- "4"
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "active_above_snow"] <- "4"
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "active_belowground"] <- "3"
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "active_below_snow"] <- "3"
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "hibernating"] <- "2"
Traits$wintering_strategy_order[Traits$wintering_strategy_order == "not_present"] <- "1"
Traits$wintering_strategy_order<-droplevels(Traits$wintering_strategy_order)
Traits$wintering_strategy_order<-as.numeric(Traits$wintering_strategy_order)

#Below ground feeding
levels(Traits$Belowground_feeding)[levels(Traits$Belowground_feeding)=='none']<-'no belowground feeding'

#Mobility
levels(Traits$Mobility)[which(levels(Traits$Mobility)=='no')]<-'immobile'
levels(Traits$Mobility)[which(levels(Traits$Mobility)=='yes')]<-'mobile'

#Use of vegetation layer #Ordered from ground to canopy
Traits$Use_of_vegetation_order<-Traits$Use_of_vegetation
levels(Traits$Use_of_vegetation_order)<-c(levels(Traits$Use_of_vegetation),1,2,3)
Traits$Use_of_vegetation_order[Traits$Use_of_vegetation_order=='ground_vegetation']<-'1'
Traits$Use_of_vegetation_order[Traits$Use_of_vegetation_order=='understory']<-'2'
Traits$Use_of_vegetation_order[Traits$Use_of_vegetation_order=='canopy']<-'3'
Traits$Use_of_vegetation_order<-droplevels(Traits$Use_of_vegetation_order)
Traits$Use_of_vegetation_order<-as.numeric(Traits$Use_of_vegetation_order)

#Convert diet items back to numeric
Traits$diet_item_forb<-as.numeric(Traits$diet_item_forb)
Traits$diet_item_shrub<-as.numeric(Traits$diet_item_shrub)
Traits$diet_item_graminoid<-as.numeric(Traits$diet_item_graminoid)
Traits$diet_item_moss<-as.numeric(Traits$diet_item_moss)
Traits$diet_item_lichen<-as.numeric(Traits$diet_item_lichen)

str(Traits)


#FAMD####
names(Traits)

var<-c( "Binomial", 
        "body_mass_log", "gut_type_order", "group_size_summer_order", "group_size_winter_order", 
        "wintering_strategy_order",  "Litter_clutch_size",  "population_dynamics_order",
        "Habitat_type", "Belowground_feeding",
        "Mobility",  "diet_type_order",   "diet_item_forb", "diet_item_graminoid", 
        "diet_item_shrub",   "diet_item_moss"  , "diet_item_lichen","Use_of_vegetation_order",
        "Elton.Plant_Others")

TOTO<-Traits[,which(colnames(Traits) %in% var)]
#TOTO<-na.omit(TOTO)
head(TOTO)
rownames(TOTO)<- TOTO$Binomial
names(TOTO)
TOTO<-TOTO[,2:19]
Traits_famd<-TOTO ; summary(Traits_famd)

res.famd<-FAMD(Traits_famd,ncp=6) # factorial analysis of mixed data

summary(res.famd)


hc_2<-HCPC(res.famd,method='ward',metric='euclidean',consol=T,nb.clust=-1) # Hierarchical Clustering on Principle Components
plot(hc_2)

borarcherb_functree<-(as.phylo(hc_2$call$t$tree))
x11(5,12)
plot(borarcherb_functree,no.margin=T)
