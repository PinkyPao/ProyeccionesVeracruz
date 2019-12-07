### Proyecto final  Analisis demografico III
### Autores: Lopez Velazquez Omar Luciano y Vazquez Castillo Paola
### Proyección de población del Estado de Veracruz de Ignacio de la Llave
### Fecha de creacion: 3 de diciembre 2019 
### Ultima modificacion: 6 de diciembre 2019 

#########################################################################
### Metodo Lee-Carter completo para hacer la ecuacion demografica basica
#########################################################################

rm(list=ls())
# gc libera memoria RAM
gc()

require(forecast)
require(ggplot2)
require(mvtnorm)

setwd("~/Mirror/Maestría/Semestre III/Analisis demografico III/Final/proyecto final")
load("Veracruz.RData")

# Tasas de mortalidad
mx<-mx2
mx<-mx[,21:46]

dev.off()
par(mfrow=c(1,2))
matplot(log(mx2[1:110,]),type="l",
        main="Hombres",
        xlab = "Edad",
        ylab="log(mx)")
matplot(log(mx2[111:220,]),type="l", 
        xlab = "Edad",
        main="Mujeres",
        ylab="log(mx)")

# Tasas de fecundidad
fx<-fx1[1:35,]
matplot(c(15:49),fx,type="l",xlab = "Edad",ylab = "Tasas específicas")
matplot(c(1970:2015),colSums(fx),type = "l",xlab = "Año",ylab = "Tasa global de fecundidad")
fx<-fx[,21:46]

# Migracion internacional
ixIH<-as.data.frame(Ix_IH_T.fit)
ixIH<-ixIH[,21:46]
ixIM<-as.data.frame(Ix_IM_T.fit)
ixIM<-ixIM[,21:46]
exIH<-as.data.frame(Ex_IH_T.fit)
exIH<-exIH[,21:46]
exIM<-as.data.frame(Ex_IM_T.fit)
exIM<-exIM[,21:46]
# Migracion nacional
ixNH<-as.data.frame(Ix_NH_T.fit)
ixNM<-as.data.frame(Ix_NM_T.fit)
exNH<-as.data.frame(Ex_NH_T.fit)
exNM<-as.data.frame(Ex_NM_T.fit)

# Definir objetos de edades y tiempos base
edades<-dim(mx)[1]
edades.fec<-dim(fx)[1]
tiempo.mort<-dim(mx)[2]
yrini.mort<-1990
yrini.fec<-1990
tiempo.fec<-dim(fx)[2]
yrbase<-2015
horizonte<-25
yrfin<-yrbase+horizonte
tiempo.tot<-tiempo.mort+horizonte
edades.mig<-dim(ixIM)[1]
tiempo.mig<-dim(ixIM)[2]
yrini.mig<-1990

lc.svd<-function(m, edades, tiempo1, tiempo2, ln){
  
  if(ln == TRUE){
    lm <- log(m)
  } else {
    lm <- m
  }
  
  ax <- rowMeans(lm[, tiempo1: tiempo2])
  
  lm_a <- lm-ax
  
  d <- matrix(0, nr = min(edades, tiempo2), nc = min(edades, tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  # Descomponiendo. Me qeudo con todas las componentes, para despues decidir con cuantas me voy a quedar
  kt <- (d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
  lc.svd <- list(ax = ax, bx = bx, kt = kt, D = d)
  
}

tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,] <- (lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,] <- Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,] <- Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(NA,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,] <- Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,] <- Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}


# Aplicando Lee-Carter a cada tasa
# Tomando el calendario de edades desde 1990
#matplot((mx),type="l")
lc.mort<- lc.svd(mx, edades, 
                 tiempo1 = 1, 
                 tiempo2 = tiempo.mort, ln=TRUE )

# Tomando el calendario de edades desde 1990
#matplot(colSums(fx),type="l")
lc.fec<- lc.svd(fx, edades.fec,
                tiempo1=1,
                tiempo2=tiempo.fec, ln=TRUE)
# Migracion Internacional
# Tomando el calendario de edades desde 1990 para homogeneizar
lc.inmF<- lc.svd(ixIM, edades.mig,
                 tiempo1=1,
                 tiempo2=tiempo.mig, ln=TRUE)
lc.inmM<- lc.svd(ixIH, edades.mig,
                 tiempo1=1,
                 tiempo2=tiempo.mig, ln=TRUE)
lc.emigF<- lc.svd(exIM, edades.mig,
                  tiempo1=1,
                  tiempo2=tiempo.mig, ln=TRUE)
lc.emigM<- lc.svd(exIH, edades.mig,
                  tiempo1=1,
                  tiempo2=tiempo.mig, ln=TRUE)

# Como la migracion nacional la vamos a mantener constante no es necesario incorporarla a los calculos de proyeccion Lee-Carter 

# Ajustando las componentes
# Mortalidad
kt1.fit<-auto.arima(lc.mort$kt[1,],trace = T,d=1)
# Fecundidad
ft1.fit<-auto.arima(lc.fec$kt[1,],trace = T,d=1,allowdrift = T)
# Migracion
# Inmigracion
it1F.fit<-auto.arima(lc.inmF$kt[1,], trace = T, allowdrift = F)
it1M.fit<-auto.arima(lc.inmM$kt[1,],trace = T, allowdrift = F)
# emigracion
et1F.fit<-auto.arima(lc.emigF$kt[1,], trace = T, allowdrift = F)
et1M.fit<-auto.arima(lc.emigM$kt[1,],trace = T, allowdrift = F)

# Proyecciones
h<-25
kt.for <- forecast(kt1.fit,h,c(95))
ft.for <- forecast(ft1.fit,h,c(95))
itF.for <- forecast(it1F.fit,h,c(95))
itM.for <- forecast(it1M.fit,h,c(95))
etF.for <- forecast(et1F.fit,h,c(95))
etM.for <- forecast(et1M.fit,h,c(95))

# Modelo de cada componente demografica
mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))
fx.for <- exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))
ixF.for <- rbind(exp(lc.inmF$ax + lc.inmF$bx[,1]%*%t(itF.for$mean)), matrix(0,20,25))
ixM.for <- rbind(exp(lc.inmM$ax + lc.inmM$bx[,1]%*%t(itM.for$mean)), matrix(0,20,25))
exF.for <- rbind(exp(lc.emigF$ax + lc.emigF$bx[,1]%*%t(etF.for$mean)), matrix(0,20,25))
exM.for <- rbind(exp(lc.emigM$ax + lc.emigM$bx[,1]%*%t(etM.for$mean)), matrix(0,20,25))


# Revisamos que si salga bien- Todo COOL :) 
matplot(log(mx.for),type = "l")
matplot((fx.for),type = "l")
matplot((ixF.for),type = "l")
#matplot((Ex_IM_T.fit),type = "l")
matplot((ixM.for),type = "l")
matplot((exF.for),type = "l")
matplot((exM.for),type = "l")

matplot((ixNH),type = "l")
#matplot((Ex_IM_T.fit),type = "l")
matplot((ixNM),type = "l")
matplot((exNH),type = "l")
matplot((exNM),type = "l")

# Hago las funciones de sobrevivencia
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex=1)$Sx
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex=2)$Sx


# Ya en este punto tengo las tasas toca proyectar a la poblacion 
# Método de Cohorte-Componente
# Primero la poblacion femenina

PxF <- as.matrix(Px[Px$Sexo=="Mujeres", -c(1,2)])
PxM <- as.matrix(Px[Px$Sexo=="Hombres", -c(1,2)])

# Matriz de la poblacion que voy a proyectar
# Poblacion a inicio de year
PxF.for <- matrix(0,110,26)
PxM.for <- matrix(0,110,26)
PxF.for[,1] <- PxF[, "2016"] 
PxM.for[,1] <-  PxM[, "2016"]
# Poblacion a mitad de periodo
NxF <- as.matrix(Nx[Nx$Sexo=="Mujeres", -c(1,2)])
NxM <- as.matrix(Nx[Nx$Sexo=="Hombres", -c(1,2)])
NxF.for <- matrix(0,110,26)
NxM.for <- matrix(0,110,26)
NxF.for[,1] <- NxF[, "2015"] 
NxM.for[,1] <-  NxM[, "2015"]
# Nacimientos
Bx <- matrix(0,35,25)
BF<- vector(length = 25)
BM <- vector(length = 25)

matplot(colSums(fx),type="l")

####################
# Poblacion de mujeres
####################

for (i in 2:26) {
  
  # Edades intermedias  
  PxF.for[2:109,i] <- (PxF.for[1:108,i-1] + NxF.for[1:108,i-1]*0.5*ixF.for[1:108,i-1]) * SxF.for[1:108,i-1] +
    NxF.for[2:109,i-1]*0.5*ixF.for[2:109,i-1] -
    NxF.for[1:108,i-1]*exF.for[1:108,i-1]
  
  ################################
  # Grupo de edades abiertas
  ################################
  
  PxF.for[110,i] <- (PxF.for[109,i-1]+
                       0.5*NxF.for[109,i-1]*ixF.for[109,i-1])* SxF.for[109,i-1] +
    NxF.for[110,i-1]*0.5*ixF.for[110,i-1]-
    NxF.for[109,i-1]*0.5*ixF.for[109,i-1]+
    (PxF.for[110,i-1] + 
       0.5*NxF.for[110,i-1]*ixF.for[110,i-1])* SxF.for[110,i-1]-NxF.for[110,i-1]*exF.for[110,i-1]
  
  ################################
  # Calcular los nacimientos
  ###############################
  
  Bx[,i-1] <- fx.for[,i-1]*(PxF.for[16:50,i-1]+
                              0.5*NxF.for[16:50,i-1]*ixF.for[16:50,i-1]+
                              PxF.for[16:50,i])*0.5
  
  # Nacimientos de mujeres
  BF[i-1]<-(1/2.05)*sum(Bx[,i-1])
  
  # Ahora poblacion femenina del primer grupo de edad
  PxF.for[1,i] <- BF[1]*SxF.for[1,i-1]+ NxF.for[1,i-1]*0.5*ixF.for[1,i-1]+ NxF.for[1,i-1]*exF.for[1,i-1]
  
  # Ya tengo el vector de mi priera poblacion proyectada completo
  
  # Calculo la nueva poblacion a mitad de year
  
  NxF.for[,i]<-0.5*(PxF.for[,i-1]+PxF.for[,i])
  
}

matplot(NxF.for,type="l")
plot(BF[1:25],type="l")

####################
# Poblacion de hombres
####################

for (i in 2:26) {
  
  ################################
  # Edades intermedias  
  ################################
  
  PxM.for[2:109,i] <- (PxM.for[1:108,i-1] + NxM.for[1:108,i-1]*0.5*ixM.for[1:108,i-1]) * SxM.for[1:108,i-1] +
    NxM.for[2:109,i-1]*0.5*ixM.for[2:109,i-1] -
    NxM.for[1:108,i-1]*exM.for[1:108,i-1]
  
  ################################
  # Grupo de edades abiertas
  ################################
  
  PxM.for[110,i] <- (PxM.for[109,i-1]+
                       0.5*NxM.for[109,i-1]*ixM.for[109,i-1])* SxM.for[109,i-1] +
    NxM.for[110,i-1]*0.5*ixM.for[110,i-1]-
    NxM.for[109,i-1]*0.5*ixM.for[109,i-1]+
    (PxM.for[110,i-1] + 
       0.5*NxM.for[110,i-1]*ixM.for[110,i-1])* SxM.for[110,i-1]-NxM.for[110,i-1]*exM.for[110,i-1]
  
  ################################
  # Calcular los nacimientos
  ###############################
  
  
  # Nacimientos de hombres
  BM[i-1]<-(1.05/2.05)*sum(Bx[,i-1])
  
  # Ahora poblacion masculina del primer grupo de edad
  PxM.for[1,i] <- BM[1]*SxM.for[1,i-1]+ NxM.for[1,i-1]*0.5*ixM.for[1,i-1]+ NxM.for[1,i-1]*exM.for[1,i-1]
  
  # Calculo la nueva poblacion a mitad de year
  
  NxM.for[,i]<-0.5*(PxM.for[,i-1]+PxM.for[,i])
  
}

PxF.for<-round(PxF.for,0)
PxM.for<-round(PxM.for,0)
NxF.for<-round(NxF.for,0)
NxM.for<-round(NxM.for,0)

#matplot((NxM.for),type="l")
colSums(NxM.for)
colSums(NxF.for)

# Estimaciones de CONAPO a 2040
# Mujeres: 4679401
# Hombres: 4294846 

