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
  #...and after doing them: it doesn't; they're all in WGS84 and plot fine, though ideally you should use an equal area projection


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


####Function: subsample.occ(), added 3-20-2015
#subsamples occurrence data, writes it to a csv file;
#can also assign function to an object and will get a SpatialPointsDataFrame object with subsampled points
#The subsampled set will be used as occurrence point input for dismo models (e.g. maxent) and also ENMeval
subsample.occ <- function(occurrences, resolution=0.5, crs="+proj=longlat +datum=WGS84", output="spdfSub.csv") {
  #Longitude and Latitude columns in occurrence dataframe must be named as that
  #resolution and projection are given default values but can be changed
  spdf <- SpatialPointsDataFrame( occurrences[ c("Longitude", "Latitude") ], data = data.frame(occurrences), proj4string = CRS(crs)) #makes a "spatial points dataframe" object
  r <- raster(spdf) #converting to raster only works if data has a spatial projection
  res(r) <- resolution
  r <- extend(r, extent(r)+1)
  spdfSub <- gridSample(spdf, r, n=1) #produces subsample of n=1 point from each raster grid; this will be occurrence data input
  write.csv(spdfSub, file=output) #writes csv with name of output arg in function
  p <- rasterToPolygons(r) #allows you to plot the grid (can't plot an empty raster)
  par(mfrow=c(1,2))
  plot(p, border='gray', main='original samples')
  points(spdf) #plots original data
  plot(p, border='gray', main='subsamples')
  points(spdfSub, cex=1, col='red', pch='x') #plots subsampled data
  return(spdfSub)
}
speasel <- subsample.occ(speaspdf, output="Speahammondii_0.5gridSubsample.csv") #will run with the default args and write a csv with default name; 
  #only required arg is dataframe with occurrences

###Getting Bioclim layers
require(raster)
BClim2_5 = getData("worldclim", var="bio", res=2.5, path="bioclim2.5/") #download 19 bioclim variables, 2.5arcmin resolution
ShamRangenarrow = extent(-125.5, -110, 25, 45) #crop Bioclim layers to this extent
bclim2.5Shamnarrow = crop(BClim2_5, ShamRangenarrow)
writeRaster(bclim2.5Shamnarrow, filename="bioclim2.5/ShamnarrowBC_2.5.grd", overwrite=T)
bclim2.5Shamnarrow = brick("bioclim2.5/ShamnarrowBC_2.5.grd")

# ##use following if you want to add elevation, though further steps will be needed for cropping it to align with bioclim layers, etc.
# countrycodes <- getData('ISO3')
# countrycodes
# Mexico = MEX, USA = USA
# 
# library(raster)
# usel <- getData("alt", country="USA")
# usel.r <- raster(usel)
# mexel <- getData("alt", country="MEX")
# mexel.r <- raster(mexel)
# elevation <- merge(usel.r, mexel.r)


# ##Code for doing model evaluation in ENMeval package
# shamoccENM <- cbind(speaLocfix[,1], speaLocfix[,2]) #have to convert long/lat to matrix in this way before running ENMevaluate
# shamnarrowbcENMeval <- ENMevaluate(speasel, bclim2.5Shamnarrow, bg.coords=pseudoabscoords, method="randomkfold", kfolds = 2)
# can use n.bg to set random background points; may be worth pursuing method in molecularecologist.com of 
# getting background points only from areas withon xx km of a presence point
# bg.coords is user-inputted pseudoabsences, which I painstakingly generated...

# data(shamnarrowbcENMeval)
# shamnarrowbcENMeval@results
# plot(shamnarrowbcENMeval@predictions[[which (shamnarrowbcENMeval@results$delta.AICc == 0) ]])
# plot(shambcENMeval@predictions[[which (shambcENMeval@results$delta.AICc == 0) ]], xlim=c(-125.5, -110), ylim=c(25, 40))
# points(shambcENMeval@occ.pts)
# shamnarrowbcENMeval@overlap #? see manual about this
# #***I should definitely subset the data before running, such that only 1 or 2 points are in any given cell...
# 
# colnames(pseudoabscoords)[1] <- "Longitude"
# colnames(pseudoabscoords)[2] <- "Latitude"


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
pastfolders <- list.dirs("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past", 
                         recursive=F)
pastfolders #amend this to a list of just the folders with .tiff climate raster layers you want to run the following function on

homedir <- "C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files"

#I have a lot of folders with a lot of different climate projection scenarios; this function helps automate
#converting them from .tiff to .bil for maxent to use
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

setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclimGF2070_rcp85")
makeprojlayers(GF_2070_rcp85)
makeprojlayers(pastfolders)

#now turn into a stack... could turn this into a function too, maybe even incorporate into the above function?
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files/bioclim_past/MIROC-LGM")
projbil <- list.files(pattern=".bil$")
projbil <- projbil[c(1,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11)] #reorder again...
projbilstack <- stack(projbil)
writeRaster(projbilstack, filename="bioclimLGM_MIROC.grd", overwrite=T) #saves stack to a file
MIROCLGMbrick <- brick("bioclimLGM_MIROC.grd")

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

bg <- randomPoints(bclim2.5Shamnarrow, 10000) #pulls random background points from extent of raster
maxentpresent_randombg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxentrun_randombg", args=c("-J", "-P"))
maxentpresent_crossval <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_crossval", args=c("-J", "-P", "replicates=5", "outputgrids=FALSE", "randomtestpoints=20"))
maxentpresent_testargs <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_testargs", randomtestpoints=20, args=c("-J", "-P"))

#runs maxent with a=backgroundpoints, path=folder to send outputs to, args: -J does jackknife, -P shows response curves
#when run with no pseudoabsence/bg arg, maxent randomly selects 1000 (or 10,000?) by default
#for projecting onto other layers/scenarios: add arg "projectionlayers=layersdirectory"
maxentfutureGD2070bg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_GD2070_rcp45_randombg", args=c("-J", "-P", "projectionlayers=bioclimGD2070"))
#have to direct the projectionlayers arg to a folder, not a RasterStack/RasterBrick object...
maxentGF2070rcp85_randombg <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_GF2070_rcp85_randombg", args=c("-J", "-P", "projectionlayers=bioclimGF2070_rcp85"))
maxentMPImidHol <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MPI_midHol", args=c("-J", "-P", "projectionlayers=bioclim_past/MPL-midhol2_5"))
maxentMPIlgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MPI_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/MPL-lgm2_5"))
maxentCCSM4lgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/CCSM4_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/CCSM4-LGM"))
maxentMIROClgm <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_past/MIROC_LGM", args=c("-J", "-P", "projectionlayers=bioclim_past/MIROC-LGM"))


par(mfrow=c(2,3)) #set plot dimensions
par(mfrow=c(1,2))

#final map output showing suitability/probability of occurrence; also in maxent.html file, or should be
#if did maxent() with more than 1 replicate, this will do it on each replicate individually rather than summarize across them
pred_me = predict(maxentrun, bclim2.5Shamnarrow) 
plot(pred_me, 1, cex=0.5, legend=T, mar=par("mar"), main="Predicted presence of western spadefoots")
#the run variable matters; the raster stack in there doesn't... I don't think...
plot(maxentrun) #plots variable contribution in the model
response(maxentrun) #plots response curves of individual variables
plot(bg)

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

predMPILGM <- predict(maxentMPIlgm, MPILGMbrick) #last glacial maximum, ~21ka
plot(predMPILGM, main="Last Glacial Maximum (~22ka), MPI model")

predCCSMLGM <- predict(maxentCCSM4lgm, CCSM4LGMbrick) #LGM with different climate model
plot(predCCSMLGM, main="Last Glacial Maximum (~22ka), CCSM4 model")

predMIROCLGM <- predict(maxentMIROClgm, MIROCLGMbrick) #LGM with yet another climate model
plot(predMIROCLGM, main="Last Glacial Maximum (~22ka), MIROC model")


#to plot by category, convert list to factors ... not relevant here but for future ref
#dumb = c("Spea", "Bufo", "Rana")
#dumbfactor <- as.factor(dumb)
#plot(speaLocfix, col=dumbfactor)
#legend("topright", legend=dumbfactor, fill=dumbfactor)


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

#testing with bootstrap replicates
maxentpresent_boots <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_boots", args=c("-J", "-P", "replicates=5", "replicatetype=bootstrap", "randomseed=TRUE", "randomtestpoints=20"))
predpres_boots <- predict(maxentpresent_boots, bclim2.5Shamnarrow)
plot(predpres_boots, main="Current modeled distribution, bootstrapping")
points(speasel, pch=4)
response(maxentpresent_boots)

#start using the "writebackgroundpredictions" option to get the lambdas file to run in... ENMeval? I think? maybe skip this
maxentpresent_boots_1rep <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_boots_1rep", args=c("-J", "-P", "replicatetype=bootstrap", "randomtestpoints=20"))
predpres_boots_1rep <- predict(maxentpresent_boots_1rep, bclim2.5Shamnarrow)
plot(predpres_boots_1rep, main="Current modeled distr., 20% bs testpts")
points(speasel, pch=4, col="red")
response(maxentpresent_boots_1rep)

#For final outputs... do I run with replicates? (no) With test samples? (not sure)
#
maxentpresent_notest <- maxent(bclim2.5Shamnarrow, speasel, a=bg, path="maxent_pres_notest", args=c("-J", "-P"))
predpres_notest <- predict(maxentpresent_notest, bclim2.5Shamnarrow)
plot(predpres_notest, main="Current modeled distr., no testpts") 
points(speasel, pch=4, col="red")
response(maxentpresent_boots_1rep)
pairs(bclim2.5Shamnarrow)

