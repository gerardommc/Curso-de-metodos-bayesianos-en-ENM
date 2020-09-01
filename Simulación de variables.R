#Gerardo Martín

#Este script genera las variables que se van a utilizar en el taller de modelación
#de nichos con métodos bayesianos. El script está basado en el tutorial para modelar 
#la autocorrelación espacial de Petr Keil (https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/)

#Comenzamos por cargar los paquetes

library(mvtnorm)   # para simular distribuciones normales multivariadas
library(raster)    # para trabajar con objetos de clase ráster
library(foreach)   # para hacer loops de manera mas avanzada
library(dismo)     # para generar las presencias

#Creamos la función que va a calcular la distancia desde cada píxel a los demán

#Esta función sólo requiere como argumento la longitud de los lados. Por el momento
#Sólo trabajaremos con un espacio cuadrado

matriz.dists <- function(long.lado)
{
      coords.filas <- rep(1:long.lado, times=long.lado)
      coords.column <- rep(1:long.lado, each=long.lado)
      filas.column <- data.frame(coords.filas, coords.column)
      D <- dist(filas.column, method="euclidean", diag=TRUE, upper=TRUE)
      D <- as.matrix(D)
      return(D)
}

#Establecemos la función que va a simular las variables (Sólamente traducida al español de la original de Petr Keil)

#

sup.correl <- function(long.lado, media.global, lambda){
      
      require(mvtnorm)
      
      D <- matriz.dists(long.lado)
      # Aquí se escala la matríz por la función exponencial (puede ser otra función!)
      #Este es el efecto del espacio, las celdas vecinas se parecen más entre si
      SIGMA <- exp(-lambda*D)
      mu <- rep(media.global, times=long.lado*long.lado)
      # Aquí generamos el contenido de la matriz que después será un objeto ráster
      # El contenido será generado con una distribución normal multivariada
      M <- matrix(nrow=long.lado, ncol=long.lado)
      M[] <- rmvnorm(1, mu, SIGMA)
      return(M)
}

matriz.a.raster <- function(matriz, extension)
{
      require(raster)
      
      rast <- raster(matriz)
      rast@extent@xmax <- extension
      rast@extent@ymax <- extension
      return(rast)
}


#######Simulando unas superficies....

#Para obtener resultados similares vamos a determinar una semilla, que es el número
# A partir del cuál se generará la serie de números aleatorios

set.seed(14392)

#Vamos a trabajar con 5 variables que representarán diversos aspectos del clima. El efecto de la
#autocorrelación espacial en cada variable será diferente y generado de manera aleatoria
#con una distribución uniforme

lambdas <- runif(5, 0.01, 0.7) # efectos de autocorrelación
long.lado <- 50 #Tamaño de las capas ráster
medias.globales <- runif(5, -20, 20) #Medias de cada capa
dists <- matriz.dists(long.lado = long.lado)

capas <- foreach(i = 1:5) %do% {
      mat <- sup.correl(long.lado = long.lado, media.global = medias.globales[i], lambda = lambdas[i])
      r <- matriz.a.raster(mat, extension = 50)
      return(r)
}

capas <- stack(capas)

saveRDS(capas, 'capas-experimento.rds')

####Ahora vamos a generar los puntos de presencia con un método del buen lichoso

##Comenzamos extraer los datos contenidos en "capas" 

capas.df <- data.frame(rasterToPoints(capas))

#Calculamos la distancia mahalanobis al centroide del nicho contenido en las cinco capas

cov = cov(capas.df[,paste0('layer.',1:5)])# Primero necesitamos la matriz de covarianza
# que controla el efecto de la correlación entre variables en la distancia

#Ahora sí calculamos las distancias al centroide dado por las medias que utilizamos
#para simular las capas
dist.mahal <- mahalanobis(capas.df[,paste0('layer.',1:5)], center = medias.globales, cov = cov) 

#Convertimos en ráster las distancias que calculamos
dist.r <- rasterFromXYZ(data.frame(capas.df[,c('x', 'y')], mahal = dist.mahal))

#Hacemos una transformación logística de las distancias, y como queremos que los bichis
#estén cerca del centroide trabajamos con la distancia * (-1), 0 es el centroide!
prob.pres <- exp(dist.r*(-1))/(1 + exp(dist.r*(-1)))

#Ahora que ya tenemos la superficie que nos establece la probabilidad de presencia
#Generamos los puntos, vamos a trabajar con 200

presencias <- randomPoints(mask = prob.pres, n = 200, prob = T)

#Guardamos las presencias simuladas
saveRDS(presencias, 'presencias.rds')


