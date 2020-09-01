setwd("~/MEGAsync/Niche models/HSPC/PA-system")

library(raster) ; library(rgdal); library(maptools); library(rgeos)
library(lgcp); library(spatstat); library(coda)


######################
#####Data loading#####
######################

##Load and reproject raster data
raster.data <- stack(c(paste('bio',c(5, 9, 12, 15),'.asc',sep=''),
                       'Maxent-p-alecto.asc'))

proj4string(raster.data) <- CRS('+init=epsg:4326') #Original CRS
raster.data <- projectRaster(raster.data, crs = CRS('+init=epsg:3577')) #Target CRS

#Load and project presence coordinates
Hendra.incidents <- read.csv('Hendra incidents-P-alecto.csv')
coordinates(Hendra.incidents) <- ~Longitude+Latitude

proj4string(Hendra.incidents) <- CRS('+init=epsg:4326')

Hendra.incidents <- spTransform(Hendra.incidents, CRS('+init=epsg:3577'))

#########################################
###Formatting the data to run analyses###
#########################################

##Computing the polygon window

#Points as a spatstat ppp object
sd.p <- readRDS('PA-sd.rds')
##Polygon overlays

cellwidth <- 21000

polyolay2 <- readRDS('PA-polyolay-ext2.rds')

ext <- polyolay2$ext


#################################################
##Poisson offset (correction for popn at risk)###
#################################################

horses <- raster('horses.tif')

horses <- mask(horses, raster.data[[1]])

raster.data <- addLayer(raster.data, horses)

names(raster.data) <- c(names(raster.data)[1:5], 'horses')

##Formatting the environmental data
covar <- data.frame(rasterToPoints(raster.data))
covar <- SpatialPixelsDataFrame(cbind(covar$x, covar$y), data = covar[,3:ncol(covar)])

##Setting up the interpolation routine for every variable type (numerical or categorical)
covar@data=guessinterp(covar@data)
covar@data <- assigninterp(df = covar@data, vars = 'horses',
                           value = "ArealWeightedSum")


##Generate the data from the computational grid and a model formula

Zmat <-getZmat(formula = X ~ bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                     I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                     bio12:bio15 + bio12:I(Maxent.p.alecto^2),
               data = sd.p, pixelcovariates = covar,
               cellwidth = cellwidth, overl = polyolay2, ext = ext)

##Spatial at risk
Zmat_pop <- getZmat(formula = X ~ horses - 1, data = sd.p,
                    pixelcovariates = covar,
                    cellwidth = cellwidth,
                    ext = ext, overl = polyolay2)

Zmat_pop <- Zmat_pop + 1

mm <- length(attr(Zmat_pop, "mcens"))
nn <- length(attr(Zmat_pop, "ncens"))

poisson.offset <- spatialAtRisk(list(X = attr(Zmat_pop, "mcens"),
                                     Y = attr(Zmat_pop, "ncens"),
                                     Zm = matrix(Zmat_pop, mm, nn)))

#################################
###Spatial covariance function###
#################################

cf <- CovFunction(SpikedExponentialCovFct)#RandomFieldsCovFct(model = 'cauchy', additionalparameters = c(2.5)))#(RandomFieldsCovFct(model="bessel",additionalparameters=1))##Best so far fractalB

###########################
###Specifying the priors###
###########################

nbeta <- ncol(Zmat)# nlayers(raster.data) #Number of predictors

priors <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = log(c(2, 5)),
                                                          variance = diag(c(2, 0.5)))),
                    betaprior = PriorSpec(GaussianPrior(mean = rep(0, nbeta),
                                                        variance = diag(rep(10^6, nbeta)))))

#################
#MCMC parameters#
#################
len.burn.thin <- c(500000, 50000, 225)# c(20000000, 2000000, 18000)

################################
######RUNNING MALA SAMPLER######
################################

ft.lgcp = lgcpPredictSpatialPlusPars(formula = X ~ bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                                           I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                                           bio12:bio15 + bio12:I(Maxent.p.alecto^2),
                                     sd = sd.p, Zmat = Zmat,
                                     model.priors = priors, 
                                     model.inits = NULL ,
                                     spatial.covmodel = cf,
                                     cellwidth = cellwidth, 
                                     poisson.offset = poisson.offset,
                                     ext = ext,
                                     output.control = setoutput(gridfunction = dump2dir(paste0('run/n'), forceSave = T),
                                                                gridmeans = MonteCarloAverage(c('mean'), lastonly = T)),
                                     mcmc.control = mcmcpars(mala.length = len.burn.thin[1], burnin = len.burn.thin[2],retain = len.burn.thin[3],
                                                             adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574))
)

saveRDS(ft.lgcp, 'Final-2-PA-system-modelfit.rds')
save.image('Final-2-PA-system.RData')
