setwd("~/MEGAsync/Niche models")

load('Future scenarios Bayesian/Future scenarios bayesian-P-alecto-OFFSET.RData')

library(raster) ; library(rgdal); library(maptools); library(rgeos)
library(lgcp); library(spatstat); library(plyr); library(coda); library(ggplot2)


######################
#####Data loading#####
######################

##Load and reproject raster data
raster.data <- stack(c(paste('BIOCLIM/Hendra/bio',c(5, 9, 12, 15),'.asc',sep=''),
                       'BIOCLIM/Hendra/Maxent-p-alecto.asc'))

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
r <- raster.data[[1]]; r[!is.na(r[])] <- 1; r <- buffer(r, 500)

window <- rasterToPolygons(r, dissolve = T)
spatstat.options(checkpolygons = F)
window <- as(window, 'owin')
window <- simplify.owin(window, dmin = 5000)      
spatstat.options(checkpolygons = T)

#Formatting presence points as a spatstat ppp object
sd.p <- ppp(x = Hendra.incidents@coords[,1], y = Hendra.incidents@coords[,2], 
            window = window) # mask = as.matrix(r))#, 

saveRDS(sd.p, 'objects/PA-sd.rds')

##Calculating optimal grid size for the point process
scale.analyses <- minimum.contrast(sd.p, model = 'exponential', method = 'g', intens = density(sd.p), transform = log)

chooseCellwidth(sd.p, cwinit = 21000)

cellwidth <- 21000

env.data <- data.frame(rasterToPoints(raster.data))
env.data <- SpatialPixelsDataFrame(cbind(env.data$x, env.data$y), env.data)

##Creating the samplig grid over whic to calculate point intensity
polyolay2 <- getpolyol(data = sd.p, pixelcovariates = env.data, cellwidth = cellwidth, ext = 4)

saveRDS(polyolay2, file = 'objects/PA-polyolay.rds')

#save(polyolay, 'Future scenarios Bayesian/poly-p-alecto.rda')
#save(polyolay2, 'Future scenarios Bayesian/poly2-p-alecto.rda')

save.image('Future scenarios Bayesian/Future scenarios bayesian-P-alecto-OFFSET.RData')

#################################################
##Poisson offset (correction for popn at risk)###
#################################################

##Loading horse data
horses <- read.csv('Complete-horse-properties.csv')
horses.x <- rep(horses$coords_x1.N.19.9, horses$number_hor.N.19.12)
horses.y <- rep(horses$coords_x2.N.19.9, horses$number_hor.N.19.12)

horses.x <- c(horses.x, coordinates(Hendra.incidents)[,1])
horses.y <- c(horses.y, coordinates(Hendra.incidents)[,2])

horses <- data.frame(x = horses.x, y = horses.y); 
#horses <- rbind(horses, data.frame(x = coordinates(Hendra.incidents)[,1], y = coordinates(Hendra.incidents)[,2]))

horses.l <- list()
for(i in 1:100){
      x = jitter(horses$x, amount = 2000); y = jitter(horses$y, amount = 2000)
      horses.l[[i]] <- data.frame(x = x, y = y)
}

r <- raster.data[[1]]#Correct the extent of the area

horses.r <- stack(lapply(horses.l, function(x){
      x <- rasterize(SpatialPoints(x), r, fun = 'count');
      x[is.na(x[])] <- 0; return(x)}))

horses.r.m <- mean(horses.r)

horses.r.m[horses.r.m[]==0] <- 0.1

writeRaster(horses.r.m, 'BIOCLIM/horses', 'GTiff', overwrite = T)

horses.r.m <- mask(horses.r.m, raster.data[[1]])

raster.data <- addLayer(raster.data, horses.r.m)

names(raster.data) <- c(names(raster.data)[1:5], 'horses')

##Formatting the environmental data
covar <- data.frame(rasterToPoints(raster.data))
covar <- SpatialPixelsDataFrame(cbind(covar$x, covar$y), data = covar[,3:ncol(covar)])

##Setting up the interpolation routine for every variable type (numerical or categorical)
covar@data=guessinterp(covar@data)
covar@data <- assigninterp(df = covar@data, vars = 'horses',
                           value = "ArealWeightedSum")


##Generate the data from the computational grid and a model formula

Zmat2 <- getZmat(formula = X ~ bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                       I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                       bio12:bio15 + bio12:I(Maxent.p.alecto^2),
                 data = sd.p, pixelcovariates = covar,
                 cellwidth = cellwidth, overl = polyolay2, ext = 4)

##Spatial at risk
Zmat_pop <- getZmat(formula = X ~ horses - 1, data = sd.p,
                    pixelcovariates = covar,
                    cellwidth = cellwidth,
                    ext = 4, overl = polyolay)

mm <- length(attr(Zmat_pop, "mcens"))
nn <- length(attr(Zmat_pop, "ncens"))

poisson.offset <- spatialAtRisk(list(X = attr(Zmat_pop, "mcens"),
                                     Y = attr(Zmat_pop, "ncens"),
                                     Zm = matrix(Zmat_pop, mm, nn)))

#################################
###Spatial covariance function###
#################################

cf <- CovFunction(RandomFieldsCovFct(model = 'cauchy', additionalparameters = c(2.5)))#(RandomFieldsCovFct(model="bessel",additionalparameters=1))##Best so far fractalB

###########################
###Specifying the priors###
###########################

nbeta <- ncol(Zmat2)# nlayers(raster.data) #Number of predictors

priors <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = log(c(2, 200)), #Better than previous
                                                          variance = diag(c(2, 5)))),
                    betaprior = PriorSpec(GaussianPrior(mean = rep(0, nbeta),
                                                        variance = diag(rep(10^6, nbeta)))))

#mean = log(c(3, 50)),
#variance = diag(c(0.15, 2))))

#mean = log(c(1.5, 200)),
#variance = diag(c(2, 2)))),

#mean = log(c(2, 200)), #Better than previous
#variance = diag(c(2, 5))))
 
#mean = log(c(2, 200)), #Mean was good, variance not so, traces of parameters were better
#variance = diag(c(3, 5)))),

#mean = log(c(2, 500)), #Slightly better traces but autocorr increased for all parameters
#variance = diag(c(3, 5))))

#mean = log(c(3, 500)), #Not good at all, lots of moves rejected
#variance = diag(c(1, 5))))

#mean = log(c(2, 500)), #Rejected moves
#variance = diag(c(2, 2))

#mean = log(c(2, 500)), #Good but autocorr of mean is very high again
#variance = diag(c(3, 10))

#mean = log(c(2, 50)), #Pretty good, seems
#variance = diag(c(0.3, 10))))

#################
#MCMC parameters#
#################
len.burn.thin <- c(5000, 500, 15) #c(2000000, 200000, 1800)

################################
######RUNNING MALA SAMPLER######
################################

library(doParallel)

registerDoParallel(cores = 3)

# foreach(i=1:3) %dopar% {

ft.lgcp = lgcpPredictSpatialPlusPars(formula = X ~ bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                                           I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                                           bio12:bio15 + bio12:I(Maxent.p.alecto^2),
                                      sd = sd.p , Zmat = Zmat2,
                                      model.priors = priors, 
                                      model.inits = NULL ,
                                      spatial.covmodel = cf,
                                      cellwidth = cellwidth, 
                                      poisson.offset = poisson.offset,
                                      ext = 4,
                                      output.control = setoutput(gridfunction = dump2dir(paste0('Future scenarios Bayesian/dumped/func/P-alecto/offset/prueba'), forceSave = T),
                                                                 gridmeans = MonteCarloAverage(c('mean'), lastonly = T)),
                                      mcmc.control = mcmcpars(mala.length = len.burn.thin[1], burnin = len.burn.thin[2],retain = len.burn.thin[3],
                                                              adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574))
)

save.image('Future scenarios Bayesian/Future scenarios bayesian-P-alecto-OFFSET.RData')

###############################
#########Diagnostics###########
###############################

par(mfrow = c(2,7));priorpost(ft.lgcp)  #Density plots of parameter values
#
par(mfrow = c(2,7));traceplots(ft.lgcp) #Traces of parameter samples
#
par(mfrow = c(2,7));parautocorr(ft.lgcp) #lagged residuals
#
parsummary(ft.lgcp) #Summary of parameter estimates

par(mfrow = c(1,1));postcov(ft.lgcp) #The fitted covariance function
#

#Autocorrelation
acor <- lgcp::autocorr(ft.lgcp, lags = c(10, 15, 19), inWindow = window)
par(mfrow = c(1, 3)); plot(acor)

#Trace plots for individual cells
par(mfrow = c(1,1))
tr <- lgcp::extract(ft.lgcp[[1]], x = 6, y = 120, t = 1, s = -1)
plot(tr, type = 'l', xlab = 'Iteration', ylab = 'Y')

#Lagged residuals
acf(tr)


plot(ft.lgcp[[1]]); plot(ft.lgcp[[2]]); plot(ft.lgcp[[3]])
plot(ft.lgcp, type = 'intensity')
points(sd.p)


###################
####Forecasting####
###################

newdata <- data.frame(rasterToPoints(raster.data))

coeffs <- ft.lgcp$betarec#rbind(ft.lgcp[[1]]$betarec, ft.lgcp[[2]]$betarec, ft.lgcp[[3]]$betarec)

colnames(coeffs) <- colnames(Zmat2)

write.csv(coeffs, 'Models coefficients/P. alecto/coeffs.csv')

Xmat <- model.matrix( ~ bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                            I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                            bio12:bio15 + bio12:I(Maxent.p.alecto^2), 
                      data = newdata)

pred <- poisson()$linkinv(coeffs[,1:ncol(Xmat)]%*%t(Xmat))

pred.interval <- adply(pred, 2, function(x){
      data.frame(Mean=mean(x), Median=median(x), HPDinterval(as.mcmc(x)))
})

pred.interval$x <- newdata$x; pred.interval$y <- newdata$y

pred.median <- rasterFromXYZ(subset(pred.interval, select = c('x', 'y', 'Median'))) * 17305600
pred.up <- rasterFromXYZ(subset(pred.interval, select = c('x', 'y', 'upper')))  * 17305600
pred.low <- rasterFromXYZ(subset(pred.interval, select = c('x', 'y', 'lower'))) * 17305600

Thr <- HPDinterval(as.mcmc(raster::extract(pred.median, Hendra.incidents)), 0.05)

par(mfrow = c(1,1))
plot(pred.median > Thr[2] , main = 'Median')
#
plot(pred.up > Thr[2], main='Upper (97.5%)')
#
plot(pred.low > Thr[2], main = 'Lower (2.5%)')
#

plot(pred.median, main = 'Median')
#
plot(pred.up, main='Upper (97.5%)')
#
plot(pred.low, main = 'Lower (2.5%)')
#

#EXceeding probabilities maps

Xmat.locs <- model.matrix(~bio5 + bio9 + bio12 + bio15 + Maxent.p.alecto + 
                                I(Maxent.p.alecto^2) + bio5:Maxent.p.alecto + bio12:Maxent.p.alecto + 
                                bio12:bio15 + bio12:I(Maxent.p.alecto^2), 
                          data = data.frame(raster::extract(raster.data, Hendra.incidents)))

pred.locs <- poisson()$linkinv(coeffs %*% t(Xmat.locs)) 

pred.locs <- adply(pred.locs, 2, function(x) {
      data.frame(Mean=mean(x), Median=median(x), HPDinterval(as.mcmc(x)))
})

#Median
prob.median <- adply(pred, 2, function(x){
      length(which(x > quantile(pred.locs$Median, 0.2)))/length(x)
})

prob.median$x <- newdata$x; prob.median$y <- newdata$y
prob.median.r <- rasterFromXYZ(subset(prob.median, select = c('x', 'y', 'V1')))

par(mfrow = c(1,1))
plot(prob.median.r)#; points(Hendra.incidents)
#

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






#All running now, but with a few errors. Now we have to see how diagnostics a done with this package
#or if other covariance structures are more appropriate than exponential



