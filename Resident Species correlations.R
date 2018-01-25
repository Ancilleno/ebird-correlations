#Dissertation code
#This r-code has been devloped by Ancilleno Davis in partial completion of the 
#requirements for the PhD in Ecology Evolution and Environmental Science at Miami University in Oxford, Ohio.

#The purpose of this Rcode is to import data from eBird for the island of Grand Bahama(GB), The Bahamas
#reorganize and analyze the data to answer the following questions

#1: How do Bahamian Resident species correlate to one another at locations on Grand Bahama Island?

#Set my working directory and load necessary libraries

# Working Directory and Libraries -----------------------------------------


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

# Specify Target Species --------------------------------------------------


#Compile list of Bahamian endemic and resident species of importance for comparison
BAHspecieslist<-c("Bahama Mockingbird", "Bahama Swallow", "Bahama Warbler", "Bahama Woodstar", "Bahama Yellowthroat","Common Ground Dove", "Northern Mockingbird", "Red-Legged Thrush")

##insert column to classify all bird species as a Bahamian Endemic or not
eBirdGBMay2016$BAHspecies<-eBirdGBMay2016$`COMMON NAME`%in% BAHspecieslist


# Subset points for Target Species ----------------------------------------

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


# Create Presence Absence data --------------------------------------------


#Add presence absence columns for each species
BahSPSwide$BAMO<-
  BahSPSwide$BASW<-
  BahSPSwide$BAWA<-
  BahSPSwide$BAWO<-
  BahSPSwide$BAYE<-
  BahSPSwide$NOMO<-0

BahSPSwide
#Create binary presence absence data for Resident Species in the data frame and as vectors
BAMOPA<-  BahSPSwide$BAMO<-BahSPSwide$`Bahama Mockingbird` >0
BASWPA<-  BahSPSwide$BASW<-BahSPSwide$`Bahama Swallow` >0
BAWAPA<-  BahSPSwide$BAWA<-BahSPSwide$`Bahama Warbler` >0
BAWOPA<-  BahSPSwide$BAWO<-BahSPSwide$`Bahama Woodstar` >0
BAYEPA<-  BahSPSwide$BAYE<-BahSPSwide$`Bahama Yellowthroat` >0
NOMOPA<-  BahSPSwide$NOMO<-BahSPSwide$`Northern Mockingbird` >0


# Use McNemar's test on presence absence data at localities --------


#USING MCNEMAR's TEST'

#McNemars test comparing all species with Northern Mockingbirds
mcnemar.test(NOMOPA, BAMOPA)
mcnemar.test(NOMOPA, BASWPA)
mcnemar.test(NOMOPA, BAWAPA)
mcnemar.test(NOMOPA, BAWOPA)
mcnemar.test(NOMOPA, BAYEPA)


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
mcnemar.test(BAWAPA, BAWOPA)
mcnemar.test(BAWAPA, BAYEPA)
mcnemar.test(BAWAPA, NOMOPA)

#Mcnemar's test comparing all species with Bahama Woodstar hummingbirds
mcnemar.test(BAWOPA, BAMOPA)
mcnemar.test(BAWOPA, BASWPA)
mcnemar.test(BAWOPA, BAWAPA)
mcnemar.test(BAWOPA, BAYEPA)
mcnemar.test(BAWOPA, NOMOPA)

##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(BAMOPA, BASWPA)
cor(BAMOPA, BAWAPA)
cor(BAMOPA, BAWOPA)
cor(BAMOPA, BAYEPA)
cor(BAMOPA, NOMOPA)

cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Mockingbird`, BahSPSwide$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Swallows

cor(BASWPA, BAMOPA)
cor(BASWPA, BAWAPA)
cor(BASWPA, BAWOPA)
cor(BASWPA, BAYEPA)
cor(BASWPA, NOMOPA)

cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Swallow`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Warblers

cor(BAWAPA, BAMOPA)
cor(BAWAPA, BASWPA)
cor(BAWAPA, BAWOPA)
cor(BAWAPA, BAYEPA)
cor(BAWAPA, NOMOPA)

cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Warbler`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Woodstar Hummingbirds

cor(BAWOPA, BAMOPA)
cor(BAWOPA, BASWPA)
cor(BAWOPA, BAWAPA)
cor(BAWOPA, BAYEPA)
cor(BAWOPA, NOMOPA)


cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Bahama Yellowthroat`)
cor(BahSPSwide$`Bahama Woodstar`, BahSPSwide$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Yellowthroats

cor(BAYEPA, BAMOPA)
cor(BAYEPA, BASWPA)
cor(BAYEPA, BAWAPA)
cor(BAYEPA, BAWOPA)
cor(BAYEPA, NOMOPA)

cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Bahama Yellowthroat`, BahSPSwide$`Northern Mockingbird`)

#correlation test comparing all species with Northern Mockingbirds

cor(NOMOPA, BAMOPA)
cor(NOMOPA, BASWPA)
cor(NOMOPA, BAWAPA)
cor(NOMOPA, BAWOPA)
cor(NOMOPA, BAYEPA)


cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Mockingbird`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Swallow`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Warbler`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Woodstar`)
cor(BahSPSwide$`Northern Mockingbird`, BahSPSwide$`Bahama Yellowthroat`)


# Species correlations by survey ------------------------------------------


#All NA for species observations during surveys were changed to "0" using Excel

# Create records by survey ------------------------------------------------


# When all species are recorded, when all species are recorded a lack of a record is a "0" for me
library(readr)
BahSPSwideSurvey<-dcast(eBirdGBMay2016BAH.sub, 
                  LOCALITY +
                    `SAMPLING EVENT IDENTIFIER`#for every locality and Sampling Event create a record
                  ~`COMMON NAME`, #each record has a column for the species common name
                  value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
write.csv(BahSPSwideSurvey, "BahSPSwideSurvey.csv")

##I used excel to replace all the NA's with 0

BahSPSwideSurvey <- read_csv("C:/Users/davisao2/Desktop/open source GIS/eBird r code/BahSPSwideSurvey.csv")
View(BahSPSwideSurvey)


# Add presence absence for each species by survey -------------------------



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
RLTHPASE<-  BahSPSwideSurvey$RLTH<-BahSPSwideSurvey$`Red-legged Thrush` >0
NOMOPASE<-  BahSPSwideSurvey$NOMO<-BahSPSwideSurvey$`Northern Mockingbird` >0


# Compare Species Presence Absence using McNemar's Test -------------------



#USING MCNEMAR's TEST'

#McNemars test comparing all species with Northern Mockingbirds
mcnemar.test(NOMOPASE, BAMOPASE)
mcnemar.test(NOMOPASE, BASWPASE)
mcnemar.test(NOMOPASE, BAWAPASE)
mcnemar.test(NOMOPASE, BAWOPASE)
mcnemar.test(NOMOPASE, BAYEPASE)

#Mcnemar's test comparing all species with Bahama Mockingbirds
mcnemar.test(BAMOPASE, BASWPASE)
mcnemar.test(BAMOPASE, BAWAPASE)
mcnemar.test(BAMOPASE, BAWOPASE)
mcnemar.test(BAMOPASE, BAYEPASE)
mcnemar.test(BAMOPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Swallows
mcnemar.test(BASWPASE, BAMOPASE)
mcnemar.test(BASWPASE, BAWAPASE)
mcnemar.test(BASWPASE, BAWOPASE)
mcnemar.test(BASWPASE, BAYEPASE)
mcnemar.test(BASWPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Warblers
mcnemar.test(BAWAPASE, BAMOPASE)
mcnemar.test(BAWAPASE, BASWPASE)
mcnemar.test(BAWAPASE, BAWOPASE)
mcnemar.test(BAWAPASE, BAYEPASE)
mcnemar.test(BAWAPASE, NOMOPASE)

#Mcnemar's test comparing all species with Bahama Woodstar hummingbirds
mcnemar.test(BAWOPASE, BAMOPASE)
mcnemar.test(BAWOPASE, BASWPASE)
mcnemar.test(BAWOPASE, BAWAPASE)
mcnemar.test(BAWOPASE, BAYEPASE)
mcnemar.test(BAWOPASE, NOMOPASE)

#McNemars test comparing all species with Red-legged Thrushes
mcnemar.test(RLTHPASE, BAMOPASE)
mcnemar.test(RLTHPASE, BASWPASE)
mcnemar.test(RLTHPASE, BAWAPASE)
mcnemar.test(RLTHPASE, BAWOPASE)
mcnemar.test(RLTHPASE, BAYEPASE)
mcnemar.test(RLTHPASE, NOMOPASE)


# Compare Species Presence Absence using Correlation ----------------------


##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(BAMOPASE, BASWPASE)
cor(BAMOPASE, BAWAPASE)
cor(BAMOPASE, BAWOPASE)
cor(BAMOPASE, BAYEPASE)
cor(BAMOPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Mockingbird`, BahSPSwideSurvey$`Northern Mockingbird`)

#correlation test comparing all species with Bahama Swallows

cor(BASWPASE, BAMOPASE)
cor(BASWPASE, BAWAPASE)
cor(BASWPASE, BAWOPASE)
cor(BASWPASE, BAYEPASE)
cor(BASWPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Warbler`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Swallow`, BahSPSwideSurvey$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Warblers

cor(BAWAPASE, BAMOPASE)
cor(BAWAPASE, BASWPASE)
cor(BAWAPASE, BAWOPASE)
cor(BAWAPASE, BAYEPASE)
cor(BAWAPASE, NOMOPASE)

cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Mockingbird`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Swallow`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Woodstar`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Bahama Yellowthroat`)
cor(BahSPSwideSurvey$`Bahama Warbler`, BahSPSwideSurvey$`Northern Mockingbird`)
#correlation test comparing all species with Bahama Woodstar Hummingbirds

cor(BAWOPASE, BAMOPASE)
cor(BAWOPASE, BASWPASE)
cor(BAWOPASE, BAWAPASE)
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


# Use Poisson Regression to create a model for NOMO using endemics --------


#Running a Poisson Regression to see if any of the endemics can be used to predict the Northern Mockingbird
summary(m1 <- glm(BahSPSwideSurvey$`Northern Mockingbird` ~ BAMOPASE + BASWPASE + BAWAPASE + BAWOPASE + BAYEPASE, family="poisson", data=BahSPSwideSurvey))


# Using a Poisson Regression predict NOMO using all residents -------------


#Running a Poisson Regression to see if any of the resident species can be used to predict the Bahama Mockingbird
summary(m1 <- glm(BahSPSwideSurvey$`Bahama Mockingbird` ~ NOMOPASE + BASWPASE + BAWAPASE + BAWOPASE + BAYEPASE, family="poisson", data=BahSPSwideSurvey))


eBirdGBMay2016BAH.sub


# Reorganize the observations by observer ---------------------------------


####By Observer###
BahSPSwideOBS<-dcast(eBirdGBMay2016BAH.sub, 
                      `OBSERVER ID`#for every locality and Sampling Event create a record
                      ~`COMMON NAME`, #each record has a column for the species common name
                      value.var = "OBSERVATION COUNT") #each Common name column includes the number of surveys in which the species was seen.
write.csv(BahSPSwideOBS, "BahSPSwideOBS.csv")

##I used excel to replace all the NA's with 0

BahSPSwideOBS <- read_csv("C:/Users/davisao2/Desktop/open source GIS/eBird r code/BahSPSwideOBS.csv")
View(BahSPSwideOBS)


# Create presence absence by observer records -----------------------------


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
RLTHPAOBS<-  BahSPSwideOBS$RLTH<-BahSPSwideOBS$`Red-legged Thrush` >0
NOMOPAOBS<-  BahSPSwideOBS$NOMO<-BahSPSwideOBS$`Northern Mockingbird` >0


# Use McNemar's test to compare observation of species by an obser --------


#USING MCNEMAR's TEST'

#McNemars test comparing all species with Northern Mockingbirds
mcnemar.test(NOMOPAOBS, BAMOPAOBS)
mcnemar.test(NOMOPAOBS, BASWPAOBS)
mcnemar.test(NOMOPAOBS, BAWAPAOBS)
mcnemar.test(NOMOPAOBS, BAWOPAOBS)
mcnemar.test(NOMOPAOBS, BAYEPAOBS)
mcnemar.test(NOMOPAOBS, NOMOPAOBS)

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


# Use Pearson's Correlation to compare species observation by obse --------


##using simple Chi squared correlations between Species
#correlation test comparing all species with Bahama Mockingbirds

cor(BAMOPAOBS, BAMOPAOBS)
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


cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Bahama Mockingbird`)
cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Bahama Swallow`)
cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Bahama Warbler`)
cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Bahama Woodstar`)
cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Bahama Yellowthroat`)
cor(BahSPSwideOBS$`Northern Mockingbird`, BahSPSwideOBS$`Northern Mockingbird`)



