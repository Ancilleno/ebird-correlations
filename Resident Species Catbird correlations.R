#Dissertation code
#This r-code has been devloped by Ancilleno Davis in partial completion of the 
#requirements for the PhD in Ecology Evolution and Environmental Science at Miami University in Oxford, Ohio.

#The purpose of this Rcode is to import data from eBird for the island of Grand Bahama(GB), The Bahamas
#reorganize and analyze the data to answer the following questions

#1: How do Bahamian Resident species correlate to one another at locations on Grand Bahama Island?

#Set my working directory and load necessary libraries
#When using McNemar's test, narrow the data down to only those that have one species or the other.
#or only localities where either species has been seen

setwd("C:/Users/davisao2/Desktop/open source GIS/eBird r code")
library(cluster)# clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(fpc)
require(ggplot2) #required for poisson regression plots
library(grDevices) #for chull
library(gridExtra)
library(gstat)
require(msm)#required for poisson regression
library(plyr) #this has the count function
library(randomForest)
library(raster)#Make sure you install randomforest before this so random forest does not overwrite.
library(readr)
library(reshape2)#contains dcast for reorganizing data into presence absence etc.
library(rgdal)
library(rgl)
library(Rlab)
require(sandwich)#required for poisson regression
library(sp)
library(spatial)
library(tidyverse)  # data manipulation
library(vegan)#visualization and species accumulation curves
#extract the raster values from the vector points at the locations for each location
#based on the classified raster data set.

require(msm)#required for poisson regression

##Classify species of concern

#Compile list of Bahamian endemic and resident species of importance for comparison
BAHspecieslist<-c("Gray Catbird","Bahama Mockingbird", "Bahama Swallow", "Bahama Warbler", "Bahama Woodstar", "Bahama Yellowthroat","Common Ground Dove", "Northern Mockingbird", "Red-legged Thrush")
BAHGRCAlist<-c("Gray Catbird","Bahama Mockingbird", "Bahama Swallow", "Bahama Warbler", "Bahama Woodstar", "Bahama Yellowthroat","Common Ground Dove", "Northern Mockingbird", "Red-Legged Thrush")

##insert column to classify all bird species as a Bahamian Endemic or not
eBirdGBMay2016$BAHspecies<-eBirdGBMay2016$`COMMON NAME`%in% BAHspecieslist

#Create subset of points for bird species in the Bahamian species list only
eBirdGBMay2016BAH.sub<-subset(eBirdGBMay2016, eBirdGBMay2016$BAHspecies=='TRUE')
eBirdGBMay2016BAH.sub

#Save my subsets to new csv files
write.csv(eBirdGBMay2016BAH.sub, file = "Bahamianspeciesebird.csv")
BahSPSwide<-dcast(eBirdGBMay2016BAH.sub, 
                   LOCALITY  #for every locality create a record
                   ~`COMMON NAME`, #each record has a column for the species common name
                   value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.

summary(BahSPSwide)


#Add presence absence columns for each species
BahSPSwide$BAMO<-
  BahSPSwide$BASW<-
  BahSPSwide$BAWA<-
  BahSPSwide$BAWO<-
  BahSPSwide$BAYE<-
  BahSPSwide$NOMO<-
  BahSPSwide$GRCA<-
  BahSPSwide$RLTH<-0

BahSPSwide
#Create binary presence absence data for Resident Species in the data frame and as vectors
BAMOPA<-  BahSPSwide$BAMO<-BahSPSwide$`Bahama Mockingbird` >0
BASWPA<-  BahSPSwide$BASW<-BahSPSwide$`Bahama Swallow` >0
BAWAPA<-  BahSPSwide$BAWA<-BahSPSwide$`Bahama Warbler` >0
BAWOPA<-  BahSPSwide$BAWO<-BahSPSwide$`Bahama Woodstar` >0
BAYEPA<-  BahSPSwide$BAYE<-BahSPSwide$`Bahama Yellowthroat` >0
GRCAPA<-  BahSPSwide$GRCA<-BahSPSwide$`Gray Catbird` >0
NOMOPA<-  BahSPSwide$NOMO<-BahSPSwide$`Northern Mockingbird` >0
RLTHPA<-  BahSPSwide$RLTH<-BahSPSwide$`Red-legged Thrush` >0

#USING MCNEMAR's TEST'
#McNemar's test specifically checks if one positive changes to a negative or negative to positive.'
#a positive and significant value would mean that the species
BAMOPA
NOMOPA
#McNemars test comparing all species with Gray Catbirds

mcnemar.test(BahSPSwide$GRCA, BahSPSwide$BAMO)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$BASW)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$BAWA)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$BAWO)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$BAYE)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$NOMO)
mcnemar.test(BahSPSwide$GRCA, BahSPSwide$RLTH)
#Mcnemar's test comparing all species with Bahama Mockingbirds
mcnemar.test(BAMOPA, BAMOPA)
mcnemar.test(BAMOPA, BASWPA)
mcnemar.test(BAMOPA, BAWAPA)
mcnemar.test(BAMOPA, BAWOPA)
mcnemar.test(BAMOPA, BAYEPA)
mcnemar.test(BAMOPA, NOMOPA)

#Mcnemar's test comparing all species with Bahama Swallows
mcnemar.test(BASWPA, BAMOPA)
mcnemar.test(BASWPA, BASWPA)
mcnemar.test(BASWPA, BAWAPA)
mcnemar.test(BASWPA, BAWOPA)
mcnemar.test(BASWPA, BAYEPA)
mcnemar.test(BASWPA, NOMOPA)

#Mcnemar's test comparing all species with Bahama Warblers
mcnemar.test(BAWAPA, BAMOPA)
mcnemar.test(BAWAPA, BASWPA)
mcnemar.test(BAWAPA, BAWAPA)
mcnemar.test(BAWAPA, BAWOPA)
mcnemar.test(BAWAPA, BAYEPA)
mcnemar.test(BAWAPA, NOMOPA)

#Mcnemar's test comparing all species with Bahama Woodstar hummingbirds
mcnemar.test(BAWOPA, BAMOPA)
mcnemar.test(BAWOPA, BASWPA)
mcnemar.test(BAWOPA, BAWAPA)
mcnemar.test(BAWOPA, BAWOPA)
mcnemar.test(BAWOPA, BAYEPA)
mcnemar.test(BAWOPA, NOMOPA)

##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(BAMOPA, BAMOPA)
cor(BAMOPA, BASWPA)
cor(BAMOPA, BAWAPA)
cor(BAMOPA, BAWOPA)
cor(BAMOPA, BAYEPA)
cor(BAMOPA, NOMOPA)

cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Swallows

cor(BASWPA, BAMOPA)
cor(BASWPA, BASWPA)
cor(BASWPA, BAWAPA)
cor(BASWPA, BAWOPA)
cor(BASWPA, BAYEPA)
cor(BASWPA, NOMOPA)

cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Warblers

cor(BAWAPA, BAMOPA)
cor(BAWAPA, BASWPA)
cor(BAWAPA, BAWAPA)
cor(BAWAPA, BAWOPA)
cor(BAWAPA, BAYEPA)
cor(BAWAPA, NOMOPA)

cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Woodstar Hummingbirds

cor(BAWOPA, BAMOPA)
cor(BAWOPA, BASWPA)
cor(BAWOPA, BAWAPA)
cor(BAWOPA, BAWOPA)
cor(BAWOPA, BAYEPA)
cor(BAWOPA, NOMOPA)


cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Yellowthroats

cor(BAYEPA, BAMOPA)
cor(BAYEPA, BASWPA)
cor(BAYEPA, BAWAPA)
cor(BAYEPA, BAWOPA)
cor(BAYEPA, BAYEPA)
cor(BAYEPA, NOMOPA)

cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Northern Mockingbird`)

#correlation test comparing all species with Northern Mockingbirds

cor(NOMOPA, BAMOPA)
cor(NOMOPA, BASWPA)
cor(NOMOPA, BAWAPA)
cor(NOMOPA, BAWOPA)
cor(NOMOPA, BAYEPA)
cor(NOMOPA, NOMOPA)


cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Northern Mockingbird`)

##Check these correlations by survey
########################################
########################################

#All NA for species observations during surveys were changed to "0" using Excel
# When all species are recorded, when all species are recorded a lack of a record is a "0" for me
library(readr)
BahSPSwideSurvey<-dcast(eBirdGBMay2016BAH.sub, 
                  LOCALITY +
                    `SAMPLING EVENT IDENTIFIER`#for every locality and Sampling Event create a record
                  ~`COMMON NAME`, #each record has a column for the species common name
                  value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
BahSPSwideSurvey[is.na(BahSPSwideSurvey)]<-0

View(BahSPSwideSurvey)



#Add presence absence columns for each species
BahSPSwideSurvey$BAMO<-
  BahSPSwideSurvey$BASW<-
  BahSPSwideSurvey$BAWA<-
  BahSPSwideSurvey$BAWO<-
  BahSPSwideSurvey$BAYE<-
  BahSPSwideSurvey$RLTH<-
  BahSPSwideSurvey$NOMO<-0

BahSPSwideSurvey
#Create binary presence absence data for Resident Species in the data frame and as vectors
BAMOPASE<-  BahSPSwideSurvey$BAMO<-BahSPSwideSurvey$`Bahama Mockingbird` >0
BASWPASE<-  BahSPSwideSurvey$BASW<-BahSPSwideSurvey$`Bahama Swallow` >0
BAWAPASE<-  BahSPSwideSurvey$BAWA<-BahSPSwideSurvey$`Bahama Warbler` >0
BAWOPASE<-  BahSPSwideSurvey$BAWO<-BahSPSwideSurvey$`Bahama Woodstar` >0
BAYEPASE<-  BahSPSwideSurvey$BAYE<-BahSPSwideSurvey$`Bahama Yellowthroat` >0
GRCAPASE<-  BahSPSwideSurvey$GRCA<-BahSPSwideSurvey$`Gray Catbird` >0
RLTHPASE<-  BahSPSwideSurvey$RLTH<-BahSPSwideSurvey$`Red-legged Thrush` >0
NOMOPASE<-  BahSPSwideSurvey$NOMO<-BahSPSwideSurvey$`Northern Mockingbird` >0

#USING MCNEMAR's TEST'

#McNemars test comparing all species with Northern Mockingbirds
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$BAMO)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$BASW)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$BAWA)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$BAWO)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$BAYE)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$NOMO)
mcnemar.test(BahSPSwideSurvey$GRCA, BahSPSwideSurvey$RLTH)

#Mcnemar's test comparing all species with Bahama Mockingbirds
mcnemar.test(BAMOPASE, BAMOPASE)
mcnemar.test(BAMOPASE, BASWPASE)
mcnemar.test(BAMOPASE, BAWAPASE)
mcnemar.test(BAMOPASE, BAWOPASE)
mcnemar.test(BAMOPASE, BAYEPASE)
mcnemar.test(BAMOPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Swallows
mcnemar.test(BASWPASE, BAMOPASE)
mcnemar.test(BASWPASE, BASWPASE)
mcnemar.test(BASWPASE, BAWAPASE)
mcnemar.test(BASWPASE, BAWOPASE)
mcnemar.test(BASWPASE, BAYEPASE)
mcnemar.test(BASWPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Warblers
mcnemar.test(BAWAPASE, BAMOPASE)
mcnemar.test(BAWAPASE, BASWPASE)
mcnemar.test(BAWAPASE, BAWAPASE)
mcnemar.test(BAWAPASE, BAWOPASE)
mcnemar.test(BAWAPASE, BAYEPASE)
mcnemar.test(BAWAPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Woodstar hummingbirds
mcnemar.test(BAWOPASE, BAMOPASE)
mcnemar.test(BAWOPASE, BASWPASE)
mcnemar.test(BAWOPASE, BAWAPASE)
mcnemar.test(BAWOPASE, BAWOPASE)
mcnemar.test(BAWOPASE, BAYEPASE)
mcnemar.test(BAWOPASE, NOMOPASE)

#McNemars test comparing all species with Red-legged Thrushes
mcnemar.test(RLTHPASE, BAMOPASE)
mcnemar.test(RLTHPASE, BASWPASE)
mcnemar.test(RLTHPASE, BAWAPASE)
mcnemar.test(RLTHPASE, BAWOPASE)
mcnemar.test(RLTHPASE, BAYEPASE)
mcnemar.test(RLTHPASE, NOMOPASE)

##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(BAMOPASE, BAMOPASE)
cor(BAMOPASE, BASWPASE)
cor(BAMOPASE, BAWAPASE)
cor(BAMOPASE, BAWOPASE)
cor(BAMOPASE, BAYEPASE)
cor(BAMOPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Swallows

cor(BASWPASE, BAMOPASE)
cor(BASWPASE, BASWPASE)
cor(BASWPASE, BAWAPASE)
cor(BASWPASE, BAWOPASE)
cor(BASWPASE, BAYEPASE)
cor(BASWPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Warblers

cor(BAWAPASE, BAMOPASE)
cor(BAWAPASE, BASWPASE)
cor(BAWAPASE, BAWAPASE)
cor(BAWAPASE, BAWOPASE)
cor(BAWAPASE, BAYEPASE)
cor(BAWAPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Woodstar Hummingbirds

cor(BAWOPASE, BAMOPASE)
cor(BAWOPASE, BASWPASE)
cor(BAWOPASE, BAWAPASE)
cor(BAWOPASE, BAWOPASE)
cor(BAWOPASE, BAYEPASE)
cor(BAWOPASE, NOMOPASE)


cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Woodstar`, BahSPSwideSurvey$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Yellowthroats

cor(BAYEPASE, BAMOPASE)
cor(BAYEPASE, BASWPASE)
cor(BAYEPASE, BAWAPASE)
cor(BAYEPASE, BAWOPASE)
cor(BAYEPASE, BAYEPASE)
cor(BAYEPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Yellowthroat`, BahSPSwideSurvey$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Mockingbirds

cor(RLTHPASE, BAMOPASE)
cor(RLTHPASE, BASWPASE)
cor(RLTHPASE, BAWAPASE)
cor(RLTHPASE, BAWOPASE)
cor(RLTHPASE, BAYEPASE)
cor(RLTHPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Red-legged Thrush`, BahSPSwideSurvey$`Northern Mockingbird`)

#correlation test comparing all species with Northern Mockingbirds

cor(NOMOPASE, BAMOPASE)
cor(NOMOPASE, BASWPASE)
cor(NOMOPASE, BAWAPASE)
cor(NOMOPASE, BAWOPASE)
cor(NOMOPASE, BAYEPASE)
cor(NOMOPASE, NOMOPASE)


cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Northern Mockingbird`, BahSPSwideSurvey$`Northern Mockingbird`)

#Running a Poisson Regression to see if any of the endemics can be used to predict the Northern Mockingbird
summary(m1 <- glm(BahSPSwideSurvey$`Northern Mockingbird` ~ BAMOPASE + BASWPASE + BAWAPASE + BAWOPASE + BAYEPASE, family="poisson", data=BahSPSwideSurvey))

#Running a Poisson Regression to see if any of the resident species can be used to predict the Bahama Mockingbird
summary(m1 <- glm(BahSPSwideSurvey$`Bahama Mockingbird` ~ NOMOPASE + BASWPASE + BAWAPASE + BAWOPASE + BAYEPASE, family="poisson", data=BahSPSwideSurvey))


eBirdGBMay2016BAH.sub

####By Observer###
BahSPSwideOBS<-dcast(eBirdGBMay2016BAH.sub, 
                      `OBSERVER ID`#for every locality and Sampling Event create a record
                      ~`COMMON NAME`, #each record has a column for the species common name
                      value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
BahSPSwideOBS[is.na(BahSPSwideOBS)]<-0

##I used excel to replace all the NA's with 0

BahSPSwideOBS <- read_csv("C:/Users/davisao2/Desktop/open source GIS/eBird r code/BahSPSwideOBS.csv")
View(BahSPSwideOBS)

#Add presence absence columns for each species
BahSPSwideOBS$BAMO<-
  BahSPSwideOBS$BASW<-
  BahSPSwideOBS$BAWA<-
  BahSPSwideOBS$BAWO<-
  BahSPSwideOBS$BAYE<-
  BahSPSwideOBS$RLTH<-
  BahSPSwideOBS$NOMO<-0

BahSPSwideOBS
#Create binary presence absence data for Resident Species in the data frame and as vectors
BAMOPAOBS<-  BahSPSwideOBS$BAMO<-BahSPSwideOBS$`Bahama Mockingbird` >0
BASWPAOBS<-  BahSPSwideOBS$BASW<-BahSPSwideOBS$`Bahama Swallow` >0
BAWAPAOBS<-  BahSPSwideOBS$BAWA<-BahSPSwideOBS$`Bahama Warbler` >0
BAWOPAOBS<-  BahSPSwideOBS$BAWO<-BahSPSwideOBS$`Bahama Woodstar` >0
BAYEPAOBS<-  BahSPSwideOBS$BAYE<-BahSPSwideOBS$`Bahama Yellowthroat` >0
GRCAPAOBS<-  BahSPSwideOBS$GRCA<-BahSPSwideOBS$`Gray Catbird` >0
RLTHPAOBS<-  BahSPSwideOBS$RLTH<-BahSPSwideOBS$`Red-legged Thrush` >0
NOMOPAOBS<-  BahSPSwideOBS$NOMO<-BahSPSwideOBS$`Northern Mockingbird` >0

#USING MCNEMAR's TEST'

#McNemars test comparing all species with Northern Mockingbirds
mcnemar.test(GRCAPAOBS, BAMOPAOBS)
mcnemar.test(GRCAPAOBS, BASWPAOBS)
mcnemar.test(GRCAPAOBS, BAWAPAOBS)
mcnemar.test(GRCAPAOBS, BAWOPAOBS)
mcnemar.test(GRCAPAOBS, BAYEPAOBS)
mcnemar.test(GRCAPAOBS, NOMOPAOBS)
mcnemar.test(GRCAPAOBS, RLTHPAOBS)
#Mcnemar's test comparing all species with Bahama Mockingbirds
mcnemar.test(BAMOPAOBS, BAMOPAOBS)
mcnemar.test(BAMOPAOBS, BASWPAOBS)
mcnemar.test(BAMOPAOBS, BAWAPAOBS)
mcnemar.test(BAMOPAOBS, BAWOPAOBS)
mcnemar.test(BAMOPAOBS, BAYEPAOBS)
mcnemar.test(BAMOPAOBS, NOMOPAOBS)

#Mcnemar's test comparing all species with Bahama Swallows
mcnemar.test(BASWPAOBS, BAMOPAOBS)
mcnemar.test(BASWPAOBS, BASWPAOBS)
mcnemar.test(BASWPAOBS, BAWAPAOBS)
mcnemar.test(BASWPAOBS, BAWOPAOBS)
mcnemar.test(BASWPAOBS, BAYEPAOBS)
mcnemar.test(BASWPAOBS, NOMOPAOBS)

#Mcnemar's test comparing all species with Bahama Warblers
mcnemar.test(BAWAPAOBS, BAMOPAOBS)
mcnemar.test(BAWAPAOBS, BASWPAOBS)
mcnemar.test(BAWAPAOBS, BAWAPAOBS)
mcnemar.test(BAWAPAOBS, BAWOPAOBS)
mcnemar.test(BAWAPAOBS, BAYEPAOBS)
mcnemar.test(BAWAPAOBS, NOMOPAOBS)

#Mcnemar's test comparing all species with Bahama Woodstar hummingbirds
mcnemar.test(BAWOPAOBS, BAMOPAOBS)
mcnemar.test(BAWOPAOBS, BASWPAOBS)
mcnemar.test(BAWOPAOBS, BAWAPAOBS)
mcnemar.test(BAWOPAOBS, BAWOPAOBS)
mcnemar.test(BAWOPAOBS, BAYEPAOBS)
mcnemar.test(BAWOPAOBS, NOMOPAOBS)

#McNemars test comparing all species with Red-legged Thrushes
mcnemar.test(RLTHPAOBS, BAMOPAOBS)
mcnemar.test(RLTHPAOBS, BASWPAOBS)
mcnemar.test(RLTHPAOBS, BAWAPAOBS)
mcnemar.test(RLTHPAOBS, BAWOPAOBS)
mcnemar.test(RLTHPAOBS, BAYEPAOBS)
mcnemar.test(RLTHPAOBS, NOMOPAOBS)

##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(GRCAPAOBS, BAMOPAOBS)
cor(BAMOPAOBS, BASWPAOBS)
cor(BAMOPAOBS, BAWAPAOBS)
cor(BAMOPAOBS, BAWOPAOBS)
cor(BAMOPAOBS, BAYEPAOBS)
cor(BAMOPAOBS, NOMOPAOBS)

cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Bahama Mockingbird`, BahSPSwideOBS$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Swallows

cor(BASWPAOBS, BAMOPAOBS)
cor(BASWPAOBS, BASWPAOBS)
cor(BASWPAOBS, BAWAPAOBS)
cor(BASWPAOBS, BAWOPAOBS)
cor(BASWPAOBS, BAYEPAOBS)
cor(BASWPAOBS, NOMOPAOBS)

cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Bahama Swallow`, BahSPSwideOBS$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Warblers

cor(BAWAPAOBS, BAMOPAOBS)
cor(BAWAPAOBS, BASWPAOBS)
cor(BAWAPAOBS, BAWAPAOBS)
cor(BAWAPAOBS, BAWOPAOBS)
cor(BAWAPAOBS, BAYEPAOBS)
cor(BAWAPAOBS, NOMOPAOBS)

cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Bahama Warbler`, BahSPSwideOBS$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Woodstar Hummingbirds

cor(BAWOPAOBS, BAMOPAOBS)
cor(BAWOPAOBS, BASWPAOBS)
cor(BAWOPAOBS, BAWAPAOBS)
cor(BAWOPAOBS, BAWOPAOBS)
cor(BAWOPAOBS, BAYEPAOBS)
cor(BAWOPAOBS, NOMOPAOBS)


cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Bahama Woodstar`, BahSPSwideOBS$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Yellowthroats

cor(BAYEPAOBS, BAMOPAOBS)
cor(BAYEPAOBS, BASWPAOBS)
cor(BAYEPAOBS, BAWAPAOBS)
cor(BAYEPAOBS, BAWOPAOBS)
cor(BAYEPAOBS, BAYEPAOBS)
cor(BAYEPAOBS, NOMOPAOBS)

cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Bahama Yellowthroat`, BahSPSwideOBS$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Mockingbirds

cor(RLTHPAOBS, BAMOPAOBS)
cor(RLTHPAOBS, BASWPAOBS)
cor(RLTHPAOBS, BAWAPAOBS)
cor(RLTHPAOBS, BAWOPAOBS)
cor(RLTHPAOBS, BAYEPAOBS)
cor(RLTHPAOBS, NOMOPAOBS)

cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Red-legged Thrush`, BahSPSwideOBS$`Northern Mockingbird`)

#correlation test comparing all species with Northern Mockingbirds

cor(NOMOPAOBS, BAMOPAOBS)
cor(NOMOPAOBS, BASWPAOBS)
cor(NOMOPAOBS, BAWAPAOBS)
cor(NOMOPAOBS, BAWOPAOBS)
cor(NOMOPAOBS, BAYEPAOBS)
cor(NOMOPAOBS, NOMOPAOBS)


cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Bahama Swallow`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Bahama Warbler`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Bahama Woodstar`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Northern Mockingbird`)
cor.test(BahSPSwideSurvey$`Gray Catbird`, BahSPSwideSurvey$`Red-legged Thrush`)


