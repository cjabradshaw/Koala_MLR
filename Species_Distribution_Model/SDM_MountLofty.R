## SDM KOALA MOUNT LOFTY RANGES
## APRIL 2024 - Frederik Saltre
# https://biomodhub.github.io/biomod2/index.html


## Clear the workspace
rm(list = ls())

## Load necessary libraries
library(biomod2)
library(terra)
library(ggplot2)
library(rgdal)
library(usdm)
library(maptools)
library(tmap)
library(sf)
library(raster)
library(eechidna)
library(ncdf4)
library(chron)
library(readr)
library(mgcv)
library(GADMTools)
library(raster)
library(matrixStats)

# load our species data
DataSpecies <- read_csv("KoalaData(GKC1&2).csv")
View(DataSpecies)
head(DataSpecies)

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 1: FORMATTING THE PRESENCE ABSENCE DATA
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Remove data points from Kangaroo Island
data.tmp1 <- subset(DataSpecies, latitude > -35.5)
data.tmp2 <- subset(data.tmp1,  longitude > 138)
DataSpecies <- data.tmp2
colnames(DataSpecies) <- c("latitude","longitude","Koala")

## Plot initial data points
gkc1.ll <- data.frame(DataSpecies$latitude, DataSpecies$longitude)
plot(gkc1.ll$DataSpecies.longitude, gkc1.ll$DataSpecies.latitude, col = 'red', pch=4, cex=0.5)
range(gkc1.ll$DataSpecies.latitude)
range(gkc1.ll$DataSpecies.longitude)

## Convert to SpatialPointsDataFrame for geographic operations
gkc1.ll <- data.frame(gkc1.ll$DataSpecies.longitude, gkc1.ll$DataSpecies.latitude)
presences <- data.frame(rep(1,dim(gkc1.ll)[1]))
gkc.sp <- SpatialPointsDataFrame(coords=gkc1.ll, proj4string=CRS("+proj=longlat +datum=WGS84"),data=presences)
gkc.proj <- project(as.matrix(coordinates(gkc.sp)), proj=(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")))

## Save projected data
write.csv(gkc.proj,"KoalaPresenceOnly_proj.csv")
DataSpecies.proj<-gkc.proj
plot(DataSpecies.proj)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 2: FORMATTING THE ENVIRONMENTAL DATA
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Set working directory for environmental data
# load our species data
setwd("~/Desktop/Koala/Model_runs/SDM/Version_2/")

# Load the projected koala presence data from a CSV file into a SpatialPointsDataFrame. 
# This structure is useful for spatial operations within R, enabling easy application of GIS operations.
DataSpecies.proj <- SpatialPointsDataFrame(read_csv("KoalaPresenceOnly_proj.csv"))
DataSpecies.proj <-DataSpecies.proj[,-1]# Remove the first column of the data frame, often an index or ID column not needed for further analysis.
colnames(DataSpecies.proj) <- c("longitude","latitude") # Rename columns to 'longitude' and 'latitude' for clarity and consistency in later spatial operations.
# Add a new column 'presence' and set all its values to 1.
# This step is typical in presence-only species distribution modeling, where only presence points are modeled,
# and absences are either later generated or not modeled directly.
DataSpecies.proj$presence <-1

# Create a new SpatialPointsDataFrame from the longitude and latitude columns.
# This involves re-specifying the coordinate system to ensure that the spatial data is accurately represented.
# The 'proj4string' specifies a projection system, in this case, an Equidistant Cylindrical projection often used for global datasets.
DataSpecies.proj2 <- SpatialPointsDataFrame(coords = DataSpecies.proj[,1:2], 
                                           proj4string=CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")),
                                          data = as.data.frame(DataSpecies.proj$presence))

## Load shapefile for South Australia
# Set the working directory to the folder containing environmental layers specifically for state shapes.
# This directory setup organizes files related to the spatial distribution modeling (SDM) of koalas.
setwd("~/Desktop/Koala/Model_runs/SDM/Version_2/Environmental_layers/States_layers")
# Load a shapefile using the readOGR function from the rgdal package. 'samap' holds the spatial data for the state,
# typically including boundaries that might be used for visualizations or to define study extents.
# 'dsn' specifies the data source name, which is the path to the shapefile without the file extension.
samap <- readOGR(dsn="sa.shp")
plot(samap)
# Transform the coordinate system of the spatial data to an Equidistant Cylindrical projection.
# This projection is often used for global scale mapping as it preserves distances along the equator.
# The transformation makes the spatial data compatible with other spatial datasets that use the same projection,
# facilitating accurate overlays and spatial operations.
sa.proj <- spTransform(samap, CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")))

## Repeat for additional environmental variables: minimum temperature, distance to water, rain, etc.
## Each layer is loaded, converted to a raster, and plotted with species points and SA map
## Load environmental layers and format them
setwd("~/Desktop/Koala/Model_runs/SDM/Version_2/Environmental_layers/")
# 1. distance to road
dist2sroads.asc <- read.asciigrid("Dist2sealRoads.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="D2SROADS")
dist2sroads.rst <- raster(dist2sroads.asc,layer=1, values=TRUE)
plot(dist2sroads.rst)
dist2sroads.xyz<-as.data.frame(dist2sroads.rst, xy=T)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 2. minimum temperature
mintemp.asc <- read.asciigrid("MinTemp.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="MINTEMP")
mintemp.rst <- raster(mintemp.asc)
mintemp.xyz<-as.data.frame(mintemp.rst, xy=T)
plot(mintemp.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 3. distance to water
dist2bwater.asc <- read.asciigrid("Dist2Water.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="D2BWATER")
dist2bwater.rst <- raster(dist2bwater.asc)
dist2bwater.xyz<-as.data.frame(dist2bwater.rst, xy=T)
plot(dist2bwater.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 4. rain
rain.asc <- read.asciigrid("Rain.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="RAIN")
rain.rst <- raster(rain.asc)
rain.xyz<-as.data.frame(rain.rst, xy=T)
plot(rain.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 5. total water index
totwatind.asc <- read.asciigrid("TotalWaterIndex.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="TOTWATIND")
totwatind.rst <- raster(totwatind.asc)
crs(totwatind.rst) <- CRS(as.character("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ptotwatind.rst <- projectRaster(totwatind.rst, crs=as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs"))
ptotwatind2.rst<-resample(ptotwatind.rst,rain.rst,method="ngb")
ptotwatind2.xyz<-as.data.frame(ptotwatind2.rst, xy=T)
plot(ptotwatind2.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 6. fraction of native cover
veg.asc <- read.asciigrid("NativeTreeCover.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="VEG")
veg.rst <- raster(veg.asc)
crs(veg.rst) <- CRS(as.character("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
pveg.rst <- projectRaster(veg.rst, crs=as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs"))
pveg2.rst<-resample(pveg.rst,rain.rst,method="ngb")
pveg2.xyz<-as.data.frame(pveg2.rst, xy=T)
plot(pveg2.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 7. solar exposure
solarexp.asc <- read.asciigrid("SolarExposure.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="SOLAREXP")
solarexp.rst <- raster(solarexp.asc)
solarexp.xyz<-as.data.frame(solarexp.rst, xy=T)
plot(solarexp.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 8. vapour pressure
vappress.asc <- read.asciigrid("VapourPressure.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="VAPPRESS")
vappress.rst <- raster(vappress.asc)
vappress.xyz<-as.data.frame(vappress.rst, xy=T)
plot(vappress.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)

# 9. topography aspect
topoasp.asc <- read.asciigrid("TopographyAspect.asc", as.image = FALSE, plot.image=T, proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="TOPOASP")
topoasp.rst <- raster(topoasp.asc)#the resolution is not the same as the other so need to convert it
crs(topoasp.rst) <- CRS(as.character("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
ptopoasp.rst <- projectRaster(topoasp.rst, crs=as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs"))
ptopoasp2.rst<-resample(ptopoasp.rst,rain.rst,method="ngb")
ptopoasp2.xyz<-as.data.frame(ptopoasp2.rst, xy=T)
plot(ptopoasp2.rst)
points(DataSpecies.proj, pch=3)
plot(sa.proj, add=T)


# Stack all environmental layers into a single raster stack
# This step combines multiple raster layers into one stack object, facilitating collective manipulation and analysis.
# The layers included are: distance to roads, minimum temperature, distance to bodies of water, rainfall,
# total water index, native vegetation cover, solar exposure, vapour pressure, and topography aspect.
env.stack <- stack(dist2sroads.rst, mintemp.rst, dist2bwater.rst, rain.rst, ptotwatind2.rst, pveg2.rst, solarexp.rst, vappress.rst, ptopoasp2.rst)
env.stack
plot(env.stack)

# put all environmental layers into a single dataframe for later plotting
env.xyz<-data.frame(dist2sroads.xyz,mintemp.xyz$MINTEMP,dist2bwater.xyz$D2BWATER,
                    rain.xyz$RAIN,ptotwatind2.xyz$TOTWATIND,pveg2.xyz$VEG,solarexp.xyz$SOLAREXP,
                    vappress.xyz$VAPPRESS,ptopoasp2.xyz$TOPOASP)
colnames(env.xyz)<-c('lon','lat','dist2road','minT','dist2water','rain','totH20','veg','H2Opress','solarexp','topo')                    
write.table(env.xyz,"MLR_EnvironmentalVariables.csv",sep = " ",dec = ".",row.names = FALSE, col.names = TRUE);


## Check for multicollinearity using the variance inflation factor
# High VIF values (generally, VIF > 10) suggest strong multicollinearity that can destabilize model estimations.
# This function calculates the VIF for each variable in the raster stack, helping to identify which variables
# might be redundant or overly correlated with others.
vif(env.stack)

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 3: BIOMOD2 SPECIES DISTRIBUTION MODELLING
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Modeling Process: Utilizes the BIOMOD_Modeling function to fit species distribution models. This function integrates several settings:
# Cross-validation strategy: Ensures that the model is robust and performs well on unseen data by using a part of the data for training and the rest for testing.
# Optimization strategy: 'bigboss' is an exhaustive approach that tries to find the best model settings, balancing model complexity and accuracy.
# Variable importance: Determines the influence of each predictor variable on the model output, which is crucial for understanding the factors driving species distribution.
# Evaluation and Variable Importance: After modeling, the script fetches and prints the evaluation scores and variable importance. These metrics are vital for assessing model quality and interpreting the ecological significance of the predictors.

## Prepare and run species distribution models using the BIOMOD2 framework
setwd("~/Desktop/Koala/Model_runs/SDM/Version_2/")
# Define the species name for the model, which in this case is 'Koala'.
myRespName <- 'Koala'
# Convert the presence/absence data for our species into a numeric format.
# This is required as BIOMOD2 models need binary (0/1) presence/absence data.
myResp <- as.numeric(DataSpecies$Koala)
# Retrieve the coordinates of species data from the preprocessed project.
# These coordinates are used to spatially reference the presence/absence data.
myRespXY <- DataSpecies.proj[,c("longitude","latitude")]

# Format the data for BIOMOD2 modeling.
# This function organizes the response (presence/absence) and explanatory (environmental) variables
# and spatial coordinates into a format suitable for analysis by BIOMOD2.
# - `resp.var`: the response variable (presence/absence)
# - `expl.var`: the explanatory variables (environmental factors)
# - `resp.xy`: coordinates of the data points
# - `resp.name`: name of the response variable
# - `PA.nb.rep`: number of repetitions for pseudo-absence generation
# - `PA.nb.absences`: number of pseudo-absences to be used
# - `PA.strategy`: strategy to select pseudo-absences ('random' in this case)
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,expl.var = env.stack,resp.xy = myRespXY,resp.name = myRespName,PA.nb.rep = 20,PA.nb.absences = 2000,PA.strategy = 'random')

myBiomodData # Display the structured data to ensure correct formatting and completeness.
plot(myBiomodData)


# 2. Defining Models Options using default options.
## Define and run models, perform evaluations, plot results, and check variable importance
setwd("~/Desktop/Koala/Model_runs/SDM/Version_2/Run")
# Perform species distribution modeling using the BIOMOD2 framework.
# `bm.format`: the formatted BIOMOD2 data set prepared earlier
# `modeling.id`: a unique identifier for this set of model runs
# `CV.strategy`: cross-validation strategy, set to 'random' to randomly assign data to folds
# `CV.nb.rep`: number of cross-validation repetitions to perform
# `CV.perc`: percentage of data used for training in each cross-validation fold
# `OPT.strategy`: optimization strategy, 'bigboss' is a comprehensive optimization approach
# `var.import`: set to 100 to compute variable importance for all models
# `metric.eval`: evaluation metrics to use, TSS (True Skill Statistic) and ROC (Receiver Operating Characteristic)
# `nb.cpu`: number of CPU cores to use for parallel processing to speed up computations = here set at 8 but will depend on the machine

# The bigboss set of parameters is available in the dataset OptionsBigboss. This set should give better results than the default set and will be continued to be optimized by the biomod2 Team.
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 100,
                                    metric.eval = c('TSS','ROC'),
                                    # seed.val = 123)
                                    nb.cpu = 8)
myBiomodModelOut
# Retrieve and print evaluation scores for each model run within the BIOMOD2 output.
# This provides insights into the performance of each model based on specified metrics.
BiomodModelEval <- get_evaluations(myBiomodModelOut)
BiomodModelVarImp<-get_variables_importance(myBiomodModelOut)

# Retrieve and print the variable importance scores from the BIOMOD2 output.
# This shows how much each environmental variable contributed to the predictive power of the models.
# This plot helps in quickly assessing model performance across different algorithms.
bm_PlotEvalMean(bm.out = myBiomodModelOut)
# Create boxplots of model evaluation scores, grouped by algorithm.
# This visualization provides insights into the variability and distribution of performance metrics across different modeling algorithms.
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
# Create boxplots of model evaluation scores, grouped by each run and algorithm.
# This helps to see how consistent the model performances are across different runs.
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
# Plot boxplots of variable importance scores, grouped by explanatory variable and algorithm.
# This illustrates which variables are most influential in the models and how their importance varies across different algorithms.
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
# Plot variable importance for each explanatory variable across different runs, grouped by algorithm.
# This detail allows for a finer comparison of variable influence under different modeling conditions.
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
# Plot variable importance for each algorithm, showing how each predictor's importance varies.
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Visualize response curves for selected models, showing the effect of each predictor at its median value.
# This plot helps understand how different environmental conditions influence the predicted presence of the species.
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'median')
# Display response curves for selected models at the minimum value of each predictor.
# This provides insights into how extreme conditions might affect species presence.
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'min')
# Plot bivariate response curves for a specific model, illustrating interactions between two predictors.
# This advanced analysis shows how combinations of environmental factors influence species distribution predictions.
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[3],
                      fixed.var = 'median',
                      do.bivariate = TRUE)

# Project single models
# project a range of models built with the BIOMOD_Modeling function onto new environmental data 
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 4: BIOMOD2 SPECIES ENSEMBLE DISTRIBUTION MODELLING
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Create ensemble models from individual SDMs, evaluate, and plot results
## Additional notes on model operations, confidence intervals, and response curve representations
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,# Input the individual model outputs from BIOMOD_Modeling as the base models.
                                      models.chosen = 'all', # Use all built models for creating the ensemble, which provides a comprehensive overview.
                                      em.by = 'all',   # Ensemble the models by considering all combinations, enhancing robustness.
                                      em.algo = c('EMwmean'), # Use weighted mean as the ensemble method, which considers model performance in weighting.
                                      metric.select = c('TSS'), # Use True Skill Statistic (TSS) to select models for the ensemble.
                                      metric.select.thresh = c(0.8),  # Set a threshold of 0.8 for TSS, only models above this threshold will be included.
                                      metric.eval = c('TSS', 'ROC'), # Evaluate models using TSS and ROC, important for understanding model accuracy and discriminative ability.
                                      var.import = 100,  # Calculate variable importance, useful for ecological insights and model interpretation.
                                      EMci.alpha = 0.05, # Set the alpha for confidence intervals in ensemble predictions, typically 0.05 for 95% confidence.
                                      EMwmean.decay = 'proportional')  # The decay factor for weights in EMwmean, ‘proportional’ balances influence proportionally to performance.
myBiomodEM  # Display the ensemble model object


# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
var_imp<-get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

toto<-myBiomodEM@variables.importance@val
write.table(toto,"Koala_SDM_Variable_Importance_CI.csv",sep = " ",dec = ".",row.names = FALSE, col.names = TRUE);

# we run the model on more time but by merging algorithms and repetitions together and considering pseudo-absence dataset individually. 
# We will later averge them to obtain a confidence interval for the response curves.
# Represent confidence interval around response curves
myBiomodEMCI <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'PA',
                                      em.algo = c('EMwmean'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.8),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 100,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')
myBiomodEMCI

# create and visualize response curves for all the models included in the ensemble model output stored 
# Represent response curves for each Pseudo Absence dataset
RcurvCI.med2<-bm_PlotResponseCurves(bm.out = myBiomodEMCI, 
                                  models.chosen = 'all',
                                  fixed.var = 'mean')
write.table(RcurvCI.med2$tab,"Koala_SDM_ResponseCurve(median)_CI.csv",sep = " ",dec = ".",row.names = FALSE, col.names = TRUE);


# Ensemble Projection 
# project ensemble models built with the BIOMOD_EnsembleModeling function onto new environmental data 
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'ROC',
                                             metric.filter = 'all',
                                             nb.cpu = 8)

BiomodEnsembleProj
myBiomodEnsembleProj <- get_predictions(BiomodEnsembleProj)#we keep only EMwmeanByTSS
plot(myBiomodEnsembleProj$Koala_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData)
plot(sa.proj, add=T)



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  STEP 5: CONVERTING HABITAT SUITABILIOTY IN POPULATION DENSITY
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##calculate population size based on population density
fred<-myBiomodEnsembleProj$Koala_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
fred2<-(fred-fred@data@min)/(fred@data@max-fred@data@min) #we standardize the projections
plot(fred2);plot(sa.proj, add=T)

fred3 <- as.data.frame(fred2, xy = T)
# here we are going to assume that the parc have no protection effect => ther is a detection bias
dens.reg<-106
fred3$dens<-fred3$Koala_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData*dens.reg
Totaldens.reg<-sum(fred3$dens, na.rm = TRUE)
Totaldens.reg #32850.55 total koala in the area

dens.reg.inf<-104
fred3$dens.inf<-fred3$Koala_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData*dens.reg.inf
Totaldens.reg.inf<-sum(fred3$dens.inf, na.rm = TRUE)
Totaldens.reg.inf #32230.73 total minimum koala in the area

dens.reg.sup<-123
fred3$dens.sup<-fred3$Koala_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData*dens.reg.sup
Totaldens.reg.sup<-sum(fred3$dens.sup, na.rm = TRUE)
Totaldens.reg.sup #38119.03 total maximum koala in the area

colnames(fred3)<-c('lon','lat','habitat_suitability','mean_density','min_density','max_density')

write.table(fred3,"Koala_density-SDM_outputs.csv",sep = " ",dec = ".",row.names = FALSE, col.names = TRUE);
