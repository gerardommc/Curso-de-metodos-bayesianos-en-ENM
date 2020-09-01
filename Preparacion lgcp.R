setwd("/home/gerardo/MEGAsync/Niche models")

load('Future scenarios Bayesian/Future scenarios bayesian-P-consp.RData')

library(raster) ; library(rgdal); library(maptools); library(dismo); library(rgeos)
library(lgcp); library(spatstat); library(plyr); library(coda)


######################
#####Data loading#####
######################

##Load and reproject raster data
raster.data <- stack(c(paste('BIOCLIM/p-consp/ascii/bio',c(2, 3, 5, 9, 13),'.asc',sep=''),
                       'Bat model results/Maxent/p-consp.asc'))


proj4string(raster.data) <- CRS('+init=epsg:4326') #Original CRS
raster.data <- projectRaster(raster.data, crs = CRS('+init=epsg:3577')) #Target CRS

##Standardise raster data
cel.mean <- cellStats(raster.data, median)
cel.mean[6] <- 0

cel.sd <- cellStats(raster.data, sd)
cel.sd[6] <- 1

raster.data <- (raster.data - cel.mean)/cel.sd

#Load and project presence coordinates
Hendra.incidents <- read.csv('Hendra incidents-P-consp.csv')
coordinates(Hendra.incidents) <- ~Longitude+Latitude

proj4string(Hendra.incidents) <- CRS('+init=epsg:4326')

Hendra.incidents <- spTransform(Hendra.incidents, CRS('+init=epsg:3577'))


#########################################
###Formatting the data to run analyses###
#########################################

##Computing the polygon window
r <- raster.data[[1]]; r[!is.na(r[])] <- 1; r <- buffer(r, 5000)

window <- rasterToPolygons(r, dissolve = T)
spatstat.options(checkpolygons = F)
window <- as(window, 'owin')
window <- simplify.owin(window, dmin = 5000)      
spatstat.options(checkpolygons = T)

#Formatting presence points as a spatstat ppp object
sd.p <- ppp(x = Hendra.incidents@coords[,1], y = Hendra.incidents@coords[,2], 
            window = window) # mask = as.matrix(r))#, 

##Calculating optimal grid size for the point process
scale.analyses <- minimum.contrast(sd.p, model = 'exponential', method = 'K', intens = density(sd.p), transform = log)

chooseCellwidth(sd.p, cwinit = 14908.3)
points(sd.p)

cellwidth <- 14908.3 

env.data <- data.frame(rasterToPoints(raster.data))
env.data <- SpatialPixelsDataFrame(cbind(env.data$x, env.data$y), env.data)

##Creating the samplig grid over whic to calculate point intensity
polyolay <- getpolyol(data = sd.p, pixelcovariates = env.data, cellwidth = cellwidth)
polyolay <- getpolyol(data = sd.p, pixelcovariates = env.data, cellwidth = cellwidth, ext = 4)

save.image('Future scenarios Bayesian/Future scenarios bayesian-P-consp.RData')

##Formatting the environmental data
covar <- data.frame(rasterToPoints(raster.data))
covar <- SpatialPixelsDataFrame(cbind(covar$x, covar$y), data = covar[,3:ncol(covar)])

##Setting up the interpolation routine for every variable type (numerical or categorical)
covar@data=guessinterp(covar@data)
covar@data <- assigninterp(df = covar@data, vars = c('bio2', 'bio3', 'bio5', 'bio9', 'bio13', 'p.consp'),
                           value = "ArealWeightedSum")

##Generate the data from the computational grid and a model formula
Zmat <- getZmat(formula = X ~ bio2 +bio3 + bio5  + bio9 + bio13 + p.consp + I(p.consp^2),
                data = sd.p, pixelcovariates = covar,
                cellwidth = cellwidth, overl = polyolay, ext = 4)

#################################
###Spatial covariance function###
#################################

cf <- CovFunction(SpikedExponentialCovFct)#RandomFieldsCovFct(model = 'matern', additionalparameters = c(1)))#Cauchy, so far the best

###########################
###Specifying the priors###
###########################

nbeta <- 2 + nlayers(raster.data) #Number of predictors

priors <-lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = log(c(2, 500)),
                                                         variance = diag(c(0.15, 2)))),
                   betaprior = PriorSpec(GaussianPrior(mean = rep(0, nbeta),
                                                       variance = diag(rep(10^6, nbeta)))))

#################
#MCMC parameters#
#################
len.burn.thin <- c(2000000, 200000, 1800) #a short chain for the beginning

################################
######RUNNING MALA SAMPLER######
################################
library(doParallel)

registerDoParallel(cores = 3)


ft.lgcp= foreach(i=1:3) %dopar% {lgcpPredictSpatialPlusPars(formula = X ~ bio2 + bio3 + bio5  + bio9 + bio13 + p.consp + I(p.consp^2),
                                   sd = sd.p , Zmat = Zmat,
                                   model.priors = priors, 
                                   model.inits = NULL ,
                                   spatial.covmodel = cf,
                                   ext = 4,
                                   cellwidth = cellwidth, 
                                   poisson.offset = NULL,
                                   output.control = setoutput(gridfunction = dump2dir(paste0('Future scenarios Bayesian/dumped/func/P-consp/',i), forceSave = T),
                                                              gridmeans = MonteCarloAverage(c('mean'), lastonly = T)),
                                   mcmc.control = mcmcpars(mala.length = len.burn.thin[1], burnin = len.burn.thin[2],retain = len.burn.thin[3],
                                                           adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574))
                              )
                        }


save.image('Future scenarios Bayesian/Future scenarios bayesian-P-consp.RData')

###############################
#########Diagnostics###########
###############################

par(mfrow = c(2,5)); priorpost(ft.lgcp[[3]])  #Density plots of parameter values
#
par(mfrow = c(2,5)); traceplots(ft.lgcp[[2]]) #Traces of parameter samples
#
par(mfrow = c(1,1)); postcov(ft.lgcp[[2]]) #The fitted covariance function
#
par(mfrow = c(2,5)); parautocorr(ft.lgcp[[1]])
#


parsummary(ft.lgcp[[1]]) #Summary of parameter estimates


#Autocorrelation
acor <- lgcp::autocorr(ft.lgcp[[1]], lags = c(10, 15, 19), inWindow = window)
par(mfrow = c(1, 3)); plot(acor)

#Trace plots for individual cells
par(mfrow = c(1,1))
tr <- lgcp::extract(ft.lgcp, x = 6, y = 120, t = 1, s = -1)
plot(tr, type = 'l', xlab = 'Iteration', ylab = 'Y')

#Lagged residuals
acf(tr)

par(mfrow = c(1,1))
plot(ft.lgcp)
plot(ft.lgcp, type = 'intensity')
points(sd.p)


###################
####Forecasting####
###################

newdata <- data.frame(rasterToPoints(raster.data))

coeffs <- rbind(ft.lgcp[[1]]$betarec, ft.lgcp[[2]]$betarec, ft.lgcp[[3]]$betarec)

Xmat <- model.matrix( ~ 1 + bio2 + bio3 + bio5 +  bio9 + bio13 + p.consp + I(p.consp^2), data = newdata)

pred <- poisson()$linkinv(coeffs[,1:8]%*%t(Xmat))

pred.intervals <- adply(pred, 2, function(x) {
      data.frame(Mean=mean(x), Median=median(x), HPDinterval(as.mcmc(x)))
})

pred.intervals$x <- newdata$x; pred.intervals$y <- newdata$y

pred.median <- rasterFromXYZ(subset(pred.intervals, select = c('x', 'y', 'Median'))) * 17305600
pred.up <- rasterFromXYZ(subset(pred.intervals, select = c('x', 'y', 'upper'))) * 17305600
pred.low <- rasterFromXYZ(subset(pred.intervals, select = c('x', 'y', 'lower'))) * 17305600

par(mfrow = c(1,1))
plot(pred.median, main = 'Median'); points(Hendra.incidents)
#
plot(pred.up, main = 'Upper (97.5%)'); points(Hendra.incidents)
#
plot(pred.low, main = 'Lower (2.5%')
#

writeRaster(pred.median, 'Future scenarios Bayesian/Model predictions/P consp/Median', 'GTiff', overwrite = T)
writeRaster(pred.up, 'Future scenarios Bayesian/Model predictions/P consp/Upper', 'GTiff', overwrite = T)
writeRaster(pred.low, 'Future scenarios Bayesian/Model predictions/P consp/Lower', 'GTiff', overwrite = T)


#EXceeding probabilities maps

Xmat.locs <- model.matrix(~1 +  bio2 + bio3 + bio5 +  bio9 + bio13 + p.consp + I(p.consp^2), 
                          data = data.frame(raster::extract(raster.data, Hendra.incidents)))

pred.locs <- poisson()$linkinv(coeffs %*% t(Xmat.locs)) 

pred.locs <- adply(pred.locs, 2, function(x) {
      data.frame(Mean=mean(x), Median=median(x), HPDinterval(as.mcmc(x)))
})

#Median
prob.median <- adply(pred, 2, function(x){
      length(which(x > quantile(pred.locs$Median, 0.8)))/length(x)
})

prob.median$x <- newdata$x; prob.median$y <- newdata$y
prob.median.r <- rasterFromXYZ(subset(prob.median, select = c('x', 'y', 'V1')))

plot(prob.median.r); points(Hendra.incidents)

#Upper
prob.95 <- adply(pred, 2, function(x){
      length(which(x > median(pred.locs$upper)))/length(x)
})
prob.95$x <- newdata$x; prob.95$y <- newdata$y
prob.95.r <- rasterFromXYZ(subset(prob.95, select = c('x', 'y', 'V1')))

plot(prob.95.r); points(Hendra.incidents)

#lower
prob.5 <- adply(pred, 2, function(x){
      length(which(x > median(pred.locs$lower)))/length(x)
})
prob.5$x <- newdata$x; prob.5$y <- newdata$y
prob.5.r <- rasterFromXYZ(subset(prob.5, select = c('x', 'y', 'V1')))

plot(prob.5.r)

save.image('Future scenarios Bayesian/Future scenarios bayesian-P-consp.RData')




