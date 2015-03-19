#submitted EEB234 project idea: Using primarily/mostly R, I will be implementing species distribution modeling 
#on Spea hammondii (western spadefoot toad) localities retrieved from GBIF. Using Maxent via the "dismo" package 
#in R I will use 19 bioclimate variables in the present and in past climate scenarios to model past distributions 
#of S. hammondii, as well as its congener Spea multiplicata, to see if their ranges once overlapped in the past to 
#allow for ancient hybridization (currently they do not overlap). Final product will be modeled species 
#distributions of both species for the Last Glacial Maximum (LGM) and Last Interglacial (LIG). with these 
#projections as models for extreme climatic scenarios of the past million years since S. hammondii and 
#S. multiplicata have diverged.
#GBIF downloads come with dozens of columns of information; it is easy to select the desired 
#columns (latitude and longitude only) using R, but I will also present a way to remove unwanted 
#columns using regular expressions.

# SDM with dismo/maxent etc.
# but first manipulate GBIF data with python to get appropriate columns and stuff. Regex etc.
# as guides: http://thebiobucket.blogspot.com/2011/11/retrieve-gbif-species-distribution-data.html
# http://www.molecularecologist.com/2013/04/species-distribution-models-in-r/
# http://cran.r-project.org/web/packages/dismo/dismo.pdf
# maybe useful: http://climate.calcommons.org/lists/datasets
#Use python+regex to edit it to species, longitude, and latitude
#then follow this SDM R procedure: http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf

#do maxent projections with JUST the southern pop of S. hammondii, and then compare with S. multiplicata projection??
#do future climate projections for whole species? for all Spea?

#According to this paper, restricted background point sampling is ineffective at reducing bias: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0097122
#instead, best to continue using random background points but to spatially subsample the occurrence points

#More useful stuff! Formatting and model testing in R. Pretty much does what I did here but better...
# https://github.com/LBAB-Humboldt/parallelMaxent

#Using R to run a Species Distribution Model with climate projections for the western spadefoot toad, _Spea hammondii_
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files")


library(ENMeval)
library(dismo)
library(maptools)
library(maps)
library(mapdata)
library(ggplot2)
library(rJava)
library(rgdal)
library(raster)
library(rworldmap)

speaGBIFraw <- read.table("Spea_hammondii_rawGBIFoutput.txt", header=TRUE, fill = TRUE, sep = "\t", quote="") #raw GBIF data requires these arguments for importing
#contains locality points for Spea hammondii, the western spadefoot toad

#This function, gbifclean(), takes raw output from GBIF or Vertnet and converts
#it into the appropriate format for running in Maxent
#this function should only be used for a single species at a time
gbifclean <- function(gbifraw) {
  latcol <- grep("decimal[Ll]atitude", colnames(gbifraw)) #returns column number for Latitude values
  longcol <- grep("decimal[Ll]ongitude", colnames(gbifraw)) #returns column number for Longitude values
  namecol <- grep("^scientificName$", colnames(gbifraw)) #returns column number for species name(s)
  latlong <- gbifraw[,c(namecol, longcol, latcol)] #subsets the input data to just these 3 columns
  colnames(latlong)[2] <- "Longitude" #changes column name to "Longitude"
  colnames(latlong)[3] <- "Latitude" #changes column name to "Latitude"
  latlong[,1] <- "Spea hammondii" #CHANGE BASED ON INPUT SPECIES! Converts names in column 1 to proper name
  latlong <- latlong[complete.cases(latlong), ] #removes rows with any missing data
  latlong[,2] <- abs(latlong[,2]) * -1 #ensures that longitude values are negative; use as appropriate for region
  latlong[,3] <- abs(latlong[,3]) #ensures that latitude values are positive; use as appropriate for region
  latlong
}
speaGBIFclean <- gbifclean(speaGBIFraw)
write.table(speaGBIFclean, file = "Speahammondii_GBIF_longlat.txt", row.names=FALSE)
head(complete.cases(speaGBIFclean))
head(speaGBIFclean)

#plotgbifclean() then takes the output from gbifclean() (or any dataframe 
#with that speciesname, longitude, latitude format) and plots it on a map
#may have to fix the projection of the "wrld_simpl" map...
plotgbifclean <- function(latlong) {
  data(wrld_simpl) #loads world map
  xlowlim <- min(latlong$Longitude) - 0 #uses values in data to set appropriate limits on plot axes
  xuplim <- max(latlong$Longitude) + 0
  ylowlim <- min(latlong$Latitude) - 0
  yuplim <- max(latlong$Latitude) + 0
  plot(wrld_simpl, xlim = c(xlowlim, xuplim), ylim = c(ylowlim, yuplim), axes=TRUE, ylab="Latitude", xlab="Longitude", cex.axis=0.7)
  points(latlong$Longitude, latlong$Latitude, col='purple', pch=20, cex=0.75) #adds points to map
  box() #puts box around plot
}
plotgbifclean(speaGBIFclean) #plotting the data may reveal aberrant points; these can be fixed later
#points look a bit off; map projection probably different from regular latlong.
#I think WorldClim data isn't projected--just latlong--so it shouldn't affect the Maxent analyses

#library(ggplot2) #I need to play around with this more. Maybe use ggmaps? Will be good to do this to 
#get used to using ggplot and to make more precise and more attractive figures/maps
#mp <- NULL #creates empty thingy
#mapWorld <- borders("world", color="gray50", fill="gray50", xlim=range(spgeo$Longitude), ylim=range(spgeo$Latitude))
#data(spgeo)
#map("usa", xlim=range(spgeo$Longitude), ylim=range(spgeo$Latitude)) + mp
#mp <- ggplot(spgeo, aes(x=Longitude, y=Latitude)) + mapWorld +geom_point(data=spgeo)
#mp
#ggplot(spealoc, aes(x=Longitude, y=Latitude)) + xlim(-125, -115) + ylim(27,39) + geom_point()
#?borders

#Using dismo! Referring here to the dismo vignette
spealoc <- speaGBIFclean[,2:3] #only need lat and long columns for maxent
head(spealoc)
plot(spealoc)
plot(wrld_simpl, add=T, border='blue', lwd=2)

#delete aberrant points
speaLocfix <- subset(spealoc, subset=((Latitude > 28.5) & (Longitude < -120 | Latitude < 39.3) & (Longitude < -118.5 | Latitude < 36.7)))
#removes points outside Spea hammondii's known range; likely other Spea species that look very similar
dim(speaLocfix)
plot(speaLocfix)
plot(wrld_simpl, add=T, border='blue', lwd=2)

#subsampling data - not required; other methods of dealing with bias
#could make this a function...
speaspdf <- SpatialPointsDataFrame( speaLocfix[ c("Longitude", "Latitude") ], data = data.frame(speaLocfix), proj4string = CRS("+proj=longlat +datum=WGS84")) #makes a "spatial points dataframe"...
r <- raster(speaspdf) #converting to raster only works if I project the localities first using line above???
r
?SpatialPointsDataFrame
speaLocfix
head(speaspdf)
res(r) <- 0.5
?res
r <- extend(r, extent(r)+1)
speasel2 <- gridSample(speaspdf, r, n=1) #produces subsample of n=1 point from each raster grid; this will be occurrence data input
p <- rasterToPolygons(r)
plot(p, border='gray')
points(speaspdf)
points(speasel2, cex=1, col='red', pch='x')
#write.table(speasel, file="") #creates a file of this subsampled data
#speasel, the subsampled set, can be used as the ENMeval input...
?gridSample

require(raster)
BClim2_5 = getData("worldclim", var="bio", res=2.5, path="bioclim2.5/") #download 19 bioclim variables, 2.5arcmin resolution
#may be easier/faster to use my pre-cropped files on other computer
ShamRange = extent(-125.5, -90.5, 15.5, 48.5) #crop Bioclim layers to this extent
#do I crop to an extent to compare with other Spea species?
ShamRangenarrow = extent(-125.5, -110, 25, 45)
bclim2.5Shamnarrow = crop(BClim2_5, ShamRangenarrow)
writeRaster(bclim2.5Shamnarrow, filename="bioclim2.5/ShamnarrowBC_2.5.grd", overwrite=T)
bclim2.5Shamnarrow = brick("bioclim2.5/ShamnarrowBC_2.5.grd")
library(ENMeval)
?extent

countrycodes <- getData('ISO3')
countrycodes
#Mexico = MEX, USA = USA

library(raster)
usel <- getData("alt", country="USA")
usel.r <- raster(usel)
mexel <- getData("alt", country="MEX")
elevation <- merge(usel, mexel)

#shamoccENM <- cbind(speaLocfix[,1], speaLocfix[,2]) #have to convert long/lat to matrix in this way before running ENMevaluate
shamnarrowbcENMeval <- ENMevaluate(speasel, bclim2.5Shamnarrow, bg.coords=pseudoabscoords, method="randomkfold", kfolds = 2)
#can use n.bg to set random background points; may be worth pursuing method in molecularecologist.com of 
#getting background points only from areas withon xx km of a presence point
#bg.coords is user-inputted pseudoabsences, which I painstakingly generated...

data(shamnarrowbcENMeval)
shamnarrowbcENMeval@results
plot(shamnarrowbcENMeval@predictions[[which (shamnarrowbcENMeval@results$delta.AICc == 0) ]])
plot(shambcENMeval@predictions[[which (shambcENMeval@results$delta.AICc == 0) ]], xlim=c(-125.5, -110), ylim=c(25, 40))
points(shambcENMeval@occ.pts)
shamnarrowbcENMeval@overlap #? see manual about this
#***I should definitely subset the data before running, such that only 1 or 2 points are in any given cell...

colnames(pseudoabscoords)[1] <- "Longitude"
colnames(pseudoabscoords)[2] <- "Latitude"
?stack

#making a stack of Future layers for climate projection
CMIP2070CN30s <- getData('CMIP5', var='bio', res=2.5, rcp=85, model='CN', year=70, path='bioclim2070_CN_rcp85_') #future; function goes to wrong URL...
#CMIP2070 <- stack(bc2070files) #to easily do stack this way, have to move the files to the working directory...
#cmip2070sham <- crop(CMIP2070, ShamRangenarrow)


#Best future projections for North America: HadCM3, CGCMA, MRI-CGCM, MPI, IPSL.. or GFDL-ESM2G???
#function to crop and rename 19 tif layers for use as projectionlayers
bcGD2070 <- getData('CMIP5', var='bio', res=2.5, rcp=45, model='GD', year=70, path='bioclimGD2070', download=FALSE)
#run this after downloading the files to the appropriate folder
#must rename projectionlayers to "bio1" etc before they will work... ugh
#temporarily change working directory
bcGD2070crop <- crop(bcGD2070, ShamRangenarrow)

#FOr PALEO reconstruction, CNRM-CM5 or IPSL-CM5A-LR appear to be best. Try IPSL first...
#Worldclim website only has a few for LGM though. If you want to use the same model for mid-Holocene and LGM, try MPI-ESM-P
#Really though, you will have to run the model for every paleoclimate model for robustness.
pastfolders <- list.dirs("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past")
pastfolders <- pastfolders[4]
homedir <- "C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files"
setwd(homedir)
getwd()
pastfolders
?list.dirs
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclimGD2070/cmip5/2-5m")
#must manually rename files or rearrange them in this list so that e.g. bio9 is the 9th item in the list; gets screwed up because of dumb sorting
makeprojlayers <- function(directories) { #function to take downloaded future layers and get them to proper format for running maxent...
  for (i in 1:length(directories)) {
    setwd(directories[i])
    projtifs <- list.files() #directory is set before running function with setwd()
    projtifs <- projtifs[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)]
  
    for (i in 1:19) {
      crop(raster(projtifs[i]), ShamRangenarrow, filename=paste("bio",i,".tif", sep=""), overwrite=TRUE)
    }
    projtifsbio <- list.files(pattern="bio")
    projtifsbio <- projtifsbio[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorders appropriately... ugh
  
    for (i in 1:19) {
      writeRaster(raster(projtifsbio[i]), filename=paste("bio",i,".bil",sep=""), format="EHdr")
   }
    if (list.files(pattern="bio1.bil$") == "bio1.bil") {
      cat(".bil rasters successfully created", sep=" ")
    }
  }
  setwd(homedir)
}
list.files("bioclimGD2070", recursive=F)
?raster
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclimGF2070_rcp85")
makeprojlayers(GF_2070_rcp85)
makeprojlayers(pastfolders)

#now turn into a stack...
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past/MIROC-LGM")
projbil <- list.files(pattern=".bil$")
projbil <- projbil[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorder again...
projbilstack <- stack(projbil)
writeRaster(projbilstack, filename="bioclimLGM_MIROC.grd", overwrite=T) #saves stack to a file
MIROCLGMbrick <- brick("bioclimLGM_MIROC.grd")
?stack
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past/MPL-lgm2_5")
projbil <- list.files(pattern=".bil$")
projbil <- projbil[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorder again...
projbilstack <- stack(projbil)
writeRaster(projbilstack, filename="bioclimLGM_MPI.grd", overwrite=T) #saves stack to a file
MPILGMbrick <- brick("bioclimLGM_MPI.grd")

setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past/CCSM4-LGM")
projbil <- list.files(pattern=".bil$")
projbil <- projbil[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorder again...
projbilstack <- stack(projbil)
writeRaster(projbilstack, filename="bioclimLGM_CCSM4.grd", overwrite=T) #saves stack to a file
CCSM4LGMbrick <- brick("bioclimLGM_CCSM4.grd")

setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclimGD2070")
projbil <- list.files(pattern=".bil$")
projbil <- projbil[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorder again...
projbilstack <- stack(projbil)
writeRaster(projbilstack, filename="bioclimGD2070_rcp45.grd", overwrite=T) #saves stack to a file
gd2070brick <- brick("bioclimGD2070_rcp45.grd")
?maxent
###
plot(CCSM4LGMbrick)
plot(projbilstack)
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files")
sham_p = kfold(speasel, 5) #vector of group assignments, splitting speasel into 5 eval groups
sham_a = kfold(pseudoabscoords, 5) #same for the background points
test = 2 #select which of the 5 groups to use as testgroup
train_p = speasel[sham_p!=test, c("Longitude", "Latitude")] #presence points for training model
train_a = pseudoabscoords[sham_a!=test, c("Longitude", "Latitude")] #absence points for training model
test_p = speasel[sham_p==test, c("Longitude", "Latitude")] #presence points for testing model
test_a = pseudoabscoords[sham_a!=test, c("Longitude", "Latitude")] #absence points for testing model

bg <- randomPoints(bclim2.5Shamnarrow, 10000) #pulls random background points from extent of raster
#maxentrun <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxentrun", args=c("-J", "-P"))
maxentpresent_randombg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxentrun_randombg", args=c("-J", "-P"))
#maxentpresent_nopseudo <- maxent(bclim2.5Shamnarrow, speasel, path="maxentrun_nopseudo", args=c("-J", "-P"))
maxentpresent_crossval <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_crossval", args=c("-J", "-P", "replicates=5", "outputgrids=FALSE", "randomtestpoints=20"))
maxentpresent_testargs <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_testargs", randomtestpoints=20, args=c("-J", "-P"))

length(speasel)
?points
#runs maxent with a=backgroundpoints, path=folder to send outputs to, args: -J does jackknife, -P shows response curves
#when run with no pseudoabsence/bg arg, maxent randomly selects 1000 (or 10,000?) by default
#for projecting onto other layers/scenarios: add arg "projectionlayers=layersdirectory"
#maxentfuture <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxentrunfuture", args=c("-J", "-P", "projectionlayers=bioclim2070_CN_rcp85"))
#maxentfuturenopseudo <- maxent(bclim2.5Shamnarrow, speasel, path="maxentrunfuturenopseudo", args=c("-J", "-P", "projectionlayers=bioclim2070_CN_rcp85"))
#maxentfutureGD2070 <- maxent(bclim2.5Shamnarrow, speasel, path="maxentrunfuturenopseudo", args=c("-J", "-P", "projectionlayers=gd2070brick"))
#maxentfutureGD2070 <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxent_GD2070_rcp45", args=c("-J", "-P", "projectionlayers=bioclimGD2070"))
maxentfutureGD2070bg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_GD2070_rcp45_randombg", args=c("-J", "-P", "projectionlayers=bioclimGD2070"))
#have to direct the projectionlayers arg to a folder, not a RasterStack/RasterBrick object...
#maxentGF2070rcp85 <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxent_GF2070_rcp85", args=c("-J", "-P", "projectionlayers=bioclimGF2070_rcp85"))
#maxentGF2070rcp85_nopseudo <- maxent(bclim2.5Shamnarrow, speasel, path="maxent_GF2070_rcp85_nopseudo", args=c("-J", "-P", "projectionlayers=bioclimGF2070_rcp85"))
maxentGF2070rcp85_randombg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_GF2070_rcp85_randombg", args=c("-J", "-P", "projectionlayers=bioclimGF2070_rcp85"))
maxentMPImidHol <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MPI_midHol", args=c("-J", "-P", "projectionlayers=bioclim_past/MPL-midhol2_5"))
maxentMPIlgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MPI_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/MPL-lgm2_5"))
maxentCCSM4lgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/CCSM4_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/CCSM4-LGM"))
maxentMIROClgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MIROC_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/MIROC-LGM"))


#run Maxent GUI to check these?
?maxent
e = evaluate(test_p, test_a, maxentrun, bclim2.5Shamnarrow) #evaluates model; info should also be in maxent.html file
e
par(mfrow=c(2,3))
par(mfrow=c(1,2))
pred_me = predict(maxentrun, bclim2.5Shamnarrow) #final map output showing suitability/probability of occurrence; also in maxent.html file, or should be
plot(pred_me, 1, cex=0.5, legend=T, mar=par("mar"), main="Predicted presence of western spadefoots")
#the run variable matters; the raster stack in there doesn't... I don't think...
?predict
plot(maxentrun)
response(maxentrun)
plot(bg)
evalGD2070 <- evaluate(test_p, test_a, maxentfutureGD2070, gd2070brick)
evalGD2070

#all of these have random background points
predGD2070bg <- predict(maxentfutureGD2070bg, gd2070brick)
plot(predGD2070bg, main="2070 projection, GD2070_rcp45 model")

predGF2070randombg <- predict(maxentGF2070rcp85_randombg, gf2070brick)
plot(predGF2070randombg, main="2070 projection, GF2070_rcp85 model")

predpresrandombg <- predict(maxentpresent_randombg, bclim2.5Shamnarrow)
plot(predpresrandombg, main="Current")
plot(maxentpresent_randombg)
response(predpresrandombg)
points(speasel)



predMPImidHol <- predict(maxentMPImidHol, MPIMidHolbrick) #mid-Hol: warmer summers, cooler winters, overall avg T below present
plot(predMPImidHol, main="mid-Holocene (~6ka) MPI model")

predMPILGM <- predict(maxentMPIlgm, MPILGMbrick)
plot(predMPILGM, main="Last Glacial Maximum (~22ka), MPI model")

predCCSMLGM <- predict(maxentCCSM4lgm, CCSM4LGMbrick)
plot(predCCSMLGM, main="Last Glacial Maximum (~22ka), CCSM4 model")

predMIROCLGM <- predict(maxentMIROClgm, MIROCLGMbrick)
plot(predMIROCLGM, main="Last Glacial Maximum (~22ka), MIROC model")

?getData
?predict
?maxent
?threshold
?response
?plot
predpres_crossval <- predict(maxentpresent_crossval, bclim2.5Shamnarrow)
plot(predpres_crossval, main="Current modeled distribution") #why don't these limits work? ugh
points(speasel, pch=4)
plot(maxentpresent_crossval)
response(maxentpresent_crossval)

predpres_testargs <- predict(maxentpresent_testargs, bclim2.5Shamnarrow)
plot(predpres_testargs)
response(maxentpresent_testargs)
#to plot by category, convert list to factors ... not relevant here but for future ref
dumb = c("Spea", "Bufo", "Rana")
dumbfactor <- as.factor(dumb)
plot(speaLocfix, col=dumbfactor)
legend("topright", legend=dumbfactor, fill=dumbfactor)

### Re-projecting rasters to equal-area
# once you get through these, go back and re-execute subsampling, maxent, etc

#install.packages("spatial.tools")
#library(spatial.tools)
#?spatial_sync_raster

###March 18
  # testing use of args in maxent()
maxentpresent_crossval <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_crossval", args=c("-J", "-P", "replicates=5", "outputgrids=FALSE"))
#maxent() doesn't accept randomtestpoints if doing crossvalidation; only if bootstrapping or subsampling
#outputgrid=FALSE doesn't work? so don't get a summary...
predpres_crossval <- predict(maxentpresent_crossval, bclim2.5Shamnarrow)
plot(predpres_crossval, main="Current modeled distribution") #why don't these limits work? ugh
points(speasel, pch=4)
plot(maxentpresent_crossval)
response(maxentpresent_crossval)

maxentpresent_boots <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_boots", args=c("-J", "-P", "replicates=5", "replicatetype=bootstrap", "randomseed=TRUE", "randomtestpoints=20"))
predpres_boots <- predict(maxentpresent_boots, bclim2.5Shamnarrow)
plot(predpres_boots, main="Current modeled distribution, bootstrapping")
points(speasel, pch=4)
response(maxentpresent_boots)

#start using the "writebackgroundpredictions" option to get the lambdas file to run in... ENMeval? I think?
maxentpresent_boots_1rep <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_boots_1rep", args=c("-J", "-P", "replicatetype=bootstrap", "randomtestpoints=20"))
predpres_boots_1rep <- predict(maxentpresent_boots_1rep, bclim2.5Shamnarrow)
plot(predpres_boots_1rep, main="Current modeled distr., 20% bs testpts")
points(speasel, pch=4, col="red")
response(maxentpresent_boots_1rep)

#For final outputs... do I run with replicates? (no) With test samples? (not sure)

maxentpresent_notest <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_notest", args=c("-J", "-P"))
predpres_notest <- predict(maxentpresent_notest, bclim2.5Shamnarrow)
plot(predpres_notest, main="Current modeled distr., no testpts") 
points(speasel, pch=4, col="red")
response(maxentpresent_boots_1rep)
pairs(bclim2.5Shamnarrow)
