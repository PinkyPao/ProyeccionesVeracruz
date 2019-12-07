### Proyecto final  Analisis demografico III
### Autores: Lopez Velazquez Omar Luciano y Vazquez Castillo Paola
### Limpieza de bases de datos
### Fecha de creacion: 25 de noviembre 2019 
### Ultima modificacion: 6 de diciembre 2019 

# Limpieza de memorias
rm(list=ls())
gc()

# Carpeta donde tenemos los archivos virgenes (descargas del dropbox) 
setwd("~/Mirror/Maestría/Semestre III/Analisis demografico III/Final/proyecto final")

# Cargamos librerias a usar
library(tidyverse)
library(readxl)
library(mosaic)
library(nls.multstart)
library(pracma)

###############
### Poblacion a mitad de periodo
###############

Nx<-read_xlsx("poblac_1950_2016.xlsx",sheet = 4)
Nx<- Nx[c(6560:6783),]
Nx<- Nx[-c(1,112,113,224),]
names(Nx)[c(2:47)]<-c(1970:2015)
Nx$Sexo<-c(rep("Hombres",110),rep("Mujeres",110))
names(Nx)[1]<-c("age")
Nx<-Nx[,c(48,1,2:47)]

###############
### Poblacion a inicio del year
###############

Px<-read_xlsx("poblac_1950_2016.xlsx",sheet = 3)
Px<- Px[c(6560:6783),]
Px<- Px[-c(1,112,113,224),]
names(Px)[c(2:48)]<-c(1970:2016)
Px$Sexo<-c(rep("Hombres",110),rep("Mujeres",110))
names(Px)[1]<-c("age")
Px<-Px[,c(49,1,2:48)]

##############
### MORTALIDAD
##############

# Hombres
mxH<-read_xlsx("tablas_mortalidad_1950_2015.xlsx",sheet="Veracruz")
mxH<-mxH[-c(1:3),3]
mxH$Sexo<-c(rep("Hombres",5152))
mxH$year<-c(rep(1970:2015, each=112))
mxH<-mxH[complete.cases(mxH[,1]),]
names(mxH)[1]<-"mx"
mxH$age<-c(rep(0:109, 46))
mxH<-mxH[,c("Sexo","age","year","mx")]
mxH<-mxH %>% spread(year,mx)

# Mujeres
mxM<-read_xlsx("tablas_mortalidad_1950_2015.xlsx",sheet="Veracruz")
mxM<-mxM[-c(1:3),11]
mxM$Sexo<-c(rep("Mujeres",5152))
mxM$year<-c(rep(1970:2015, each=112))
mxM<-mxM[complete.cases(mxM[,1]),]
names(mxM)[1]<-"mx"
mxM$age<-c(rep(0:109, 46))
mxM<-mxM[,c("Sexo","age","year","mx")]
mxM<-mxM %>% spread(year,mx)
# Tasas de mortalidad juntas
mx<-bind_rows(mxH,mxM)
# matplot(mxH[,-c(1,2)],type="l")
# matplot(mxM[,-c(1,2)],type="l")
# matplot(mx[,-c(1,2)],type="l")
mx2<-as.data.frame(lapply(mx[,-c(1,2)],as.numeric))
matplot(log(mx2),type="l")

###############
### FECUNDIDAD
###############

fec<-read.csv("fecundidad_1950_2015.csv", header = T)
fec<-fec[c(296:303),]
names(fec)<-c("age",1970:2015)

# Ahora para proyectar necesitamos tener las tasas especificas de fecundidad 
# Vamos a calcularlas con el metodo de Heather-Booth

# Primero hay que poner la base en forma tidy para que sea compatible con la funcion de ASFR que ya tenemos de clase

fx5<- fec[-8,]
fx5<-fx5 %>% gather("year","fx",-age)
fx5<-fx5[,c("year","age","fx")]

# Ahora si- corro la funcion de descoposicion de ASFR (age-specific-fertility-rates)
asfr <-function(fx5,year){
  fx0<-fx5[fx5$year==year,"fx"]
  Fx<-5*cumsum(fx0)
  TGF<-Fx[7]
  #Proporciones de las fecundidades acumuladas con respecto a la TGF
  FxF<-Fx/TGF
  # Edaes quinquenales (marcas de clase de los intervalos) 
  x5<-seq(17.5,47.5,5)
  ### Calculo la funcion de Booth
  Yx<-log(-log(FxF))
  
  ### Ahora defino la regresion. 
  # Regresion de los valores que tengo de la F(x) de Booth y mis marcas de clase en las edades
  Yx.lm<-lm(Yx[-7]~x5[-7])
  a<-as.numeric(Yx.lm$coefficients[1])
  b<- as.numeric(Yx.lm$coefficients[2])
  A <- as.numeric(exp(-exp(a)))
  # A es la proporcion de la fecundidad alcanzada a la edad media a la fecundidad
  B <- as.numeric(exp(b))
  # B es la varianza de A
  # Genero los puntos donde quiero estimar la fec por edades simples
  x1<-c(15:50)
  
  Fx.estim<-as.numeric(TGF*(A^(B^(x1))))
  
  fx<- Fx.estim[2:36]-Fx.estim[1:35]
  fx<- c(Fx.estim[1],fx)
  return(fx)
}

# Creo una matriz para guardar las edades especificas
fx1<-data.frame(matrix(0,36,46))
row.names(fx1)<-c(15:50)
names(fx1)<-c(1970:2015)
# Llenando la matriz
for (i in 1:46) {
  fx1[,i]<-asfr(fx5,1969+i)
}
# Tasas especificas de fecundidad ASFR 
matplot(fx1,type = "l")

# Tasas globales de fecundidad
# TGF.fx1<-colSums(fx1)
# plot(TGF.fx1,type="l")
# Ya tengo las tasas especificas de fecundidad 

############## 
### MIGRACION  
##############

# Migracion internacional 
# Inmigrantes 
# Hombres
Ix_IH<-read_xlsx("mig_internacion_1970_2015.xlsx",sheet = 3)
Ix_IH<-Ix_IH[1227:1244,]
# Numero
Ix_IH_N<-Ix_IH[,1:10]
# Tasa
Ix_IH_T<-Ix_IH[,c(1,12:20)]
Ix_IH_T2<-as.data.frame(lapply(Ix_IH_T[,-1],as.numeric))/1000
#matplot(Ix_IH_T2,type="l")

Ix_IM<-read_xlsx("mig_internacion_1970_2015.xlsx",sheet = 3)
Ix_IM<-Ix_IM[1247:1264,]
# Numero
Ix_IM_N<-Ix_IM[,1:10]
# Tasa
Ix_IM_T<-Ix_IM[,c(1,12:20)]
Ix_IM_T2<-as.data.frame(lapply(Ix_IM_T[,-1],as.numeric))/1000
# Asi se ve la inmigracion internacional de mujeres en 2010-2015
#matplot(Ix_IM_T2,type="l")

# Emigrantes
# Hombres
Ex_IH<-read_xlsx("mig_internacion_1970_2015.xlsx",sheet = 4)
Ex_IH<-Ex_IH[1227:1244,]
Ex_IH_N<-Ex_IH[,1:10]
Ex_IH_T<-Ex_IH[,c(1,12:20)]
Ex_IH_T2<-as.data.frame(lapply(Ex_IH_T[,-1],as.numeric))/1000
#matplot(Ex_IH_T2,type = "l")

#Mujeres
Ex_IM<-read_xlsx("mig_internacion_1970_2015.xlsx",sheet = 4)
Ex_IM<-Ex_IM[1247:1264,]
Ex_IM_N<-Ex_IM[,1:10]
Ex_IM_T<-Ex_IM[,c(1,12:20)]
Ex_IM_T2<-as.data.frame(lapply(Ex_IM_T[,-1],as.numeric))/1000
#matplot(Ex_IM_T2,type = "l")

# Migracion interna

Mig_mat_70<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 2)
Mig_mat_75<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 3)
Mig_mat_80<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 4)
Mig_mat_85<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 5)
Mig_mat_90<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 6)
Mig_mat_95<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 7)
Mig_mat_00<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 8)
Mig_mat_05<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 9)
Mig_mat_10<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 10)

# Inmigracion Nacional
# Veracruz es la columna 31 y las filas 1228-1245 y 1248-1266
#1970
Ix_NH_70<-as.data.frame(lapply(Mig_mat_70[1228:1245,-1],as.numeric))
Ix_NH_70<-rowSums(Ix_NH_70)
Ix_NM_70<-as.data.frame(lapply(Mig_mat_70[1248:1265,-1],as.numeric))
Ix_NM_70<-rowSums(Ix_NM_70)

#1975
Ix_NH_75<-as.data.frame(lapply(Mig_mat_75[1228:1245,-1],as.numeric))
Ix_NH_75<-rowSums(Ix_NH_75)
Ix_NM_75<-as.data.frame(lapply(Mig_mat_75[1248:1265,-1],as.numeric))
Ix_NM_75<-rowSums(Ix_NM_75)

#1980
Ix_NH_80<-as.data.frame(lapply(Mig_mat_80[1228:1245,-1],as.numeric))
Ix_NH_80<-rowSums(Ix_NH_80)
Ix_NM_80<-as.data.frame(lapply(Mig_mat_80[1248:1265,-1],as.numeric))
Ix_NM_80<-rowSums(Ix_NM_80)

#1985
Ix_NH_85<-as.data.frame(lapply(Mig_mat_85[1228:1245,-1],as.numeric))
Ix_NH_85<-rowSums(Ix_NH_85)
Ix_NM_85<-as.data.frame(lapply(Mig_mat_85[1248:1265,-1],as.numeric))
Ix_NM_85<-rowSums(Ix_NM_85)

#1990
Ix_NH_90<-as.data.frame(lapply(Mig_mat_90[1228:1245,-1],as.numeric))
Ix_NH_90<-rowSums(Ix_NH_90)
Ix_NM_90<-as.data.frame(lapply(Mig_mat_90[1248:1265,-1],as.numeric))
Ix_NM_90<-rowSums(Ix_NM_90)

#1995
Ix_NH_95<-as.data.frame(lapply(Mig_mat_95[1228:1245,-1],as.numeric))
Ix_NH_95<-rowSums(Ix_NH_95)
Ix_NM_95<-as.data.frame(lapply(Mig_mat_95[1248:1265,-1],as.numeric))
Ix_NM_95<-rowSums(Ix_NM_95)

#2000
Ix_NH_00<-as.data.frame(lapply(Mig_mat_00[1228:1245,-1],as.numeric))
Ix_NH_00<-rowSums(Ix_NH_00)
Ix_NM_00<-as.data.frame(lapply(Mig_mat_00[1248:1265,-1],as.numeric))
Ix_NM_00<-rowSums(Ix_NM_00)

#2005
Ix_NH_05<-as.data.frame(lapply(Mig_mat_05[1228:1245,-1],as.numeric))
Ix_NH_05<-rowSums(Ix_NH_05)
Ix_NM_05<-as.data.frame(lapply(Mig_mat_05[1248:1265,-1],as.numeric))
Ix_NM_05<-rowSums(Ix_NM_05)

#2010
Ix_NH_10<-as.data.frame(lapply(Mig_mat_10[1228:1245,-1],as.numeric))
Ix_NH_10<-rowSums(Ix_NH_10)
Ix_NM_10<-as.data.frame(lapply(Mig_mat_10[1248:1265,-1],as.numeric))
Ix_NM_10<-rowSums(Ix_NM_10)

# Matriz Inmigracion nacional

Ix_NH_N<-cbind(Mig_mat_70[1228:1245,1],
             Ix_NH_70,Ix_NH_75,Ix_NH_80,Ix_NH_85,Ix_NH_90,Ix_NH_95,Ix_NH_00,Ix_NH_05,Ix_NH_10)
Ix_NM_N<-cbind(Mig_mat_70[1228:1245,1],
               Ix_NM_70,Ix_NM_75,Ix_NM_80,Ix_NM_85,Ix_NM_90,Ix_NM_95,Ix_NM_00,Ix_NM_05,Ix_NM_10)

# Para estudiar el modelo de Rogers-Castro se necesitan tasas asi que
APV_migN<-read_xlsx("mig_internos_1970_2015.xlsx",sheet = 1)
APV_migNH<-as.data.frame(lapply(APV_migN[1228:1245,2:10],as.numeric))
APV_migNM<-as.data.frame(lapply(APV_migN[1248:1265,2:10],as.numeric))

# calcular las tasas- Son del orden correcto asi que estan bien :)
Ix_NH_T<-Ix_NH_N[,-1]/APV_migNH
Ix_NM_T<-Ix_NM_N[,-1]/APV_migNM
# matplot(Ix_NH_T,type="l")
# matplot(Ix_NM_T,type="l")

# Emigracion nacional

# Extraigo la informacion de cada year

Ex_NH_70a<-as.data.frame(lapply(Mig_mat_70[,31],as.numeric))
Ex_NH_75a<-as.data.frame(lapply(Mig_mat_75[,31],as.numeric))
Ex_NH_80a<-as.data.frame(lapply(Mig_mat_80[,31],as.numeric))
Ex_NH_85a<-as.data.frame(lapply(Mig_mat_85[,31],as.numeric))
Ex_NH_90a<-as.data.frame(lapply(Mig_mat_90[,31],as.numeric))
Ex_NH_95a<-as.data.frame(lapply(Mig_mat_95[,31],as.numeric))
Ex_NH_00a<-as.data.frame(lapply(Mig_mat_00[,31],as.numeric))
Ex_NH_05a<-as.data.frame(lapply(Mig_mat_05[,31],as.numeric))
Ex_NH_10a<-as.data.frame(lapply(Mig_mat_10[,31],as.numeric))

# Matriz que voy a llenar con las filas
Ex_NH_70<-matrix(0,18,32)
Ex_NH_75<-matrix(0,18,32)
Ex_NH_80<-matrix(0,18,32)
Ex_NH_85<-matrix(0,18,32)
Ex_NH_90<-matrix(0,18,32)
Ex_NH_95<-matrix(0,18,32)
Ex_NH_00<-matrix(0,18,32)
Ex_NH_05<-matrix(0,18,32)
Ex_NH_10<-matrix(0,18,32)

Ex_NM_70<-matrix(0,18,32)
Ex_NM_75<-matrix(0,18,32)
Ex_NM_80<-matrix(0,18,32)
Ex_NM_85<-matrix(0,18,32)
Ex_NM_90<-matrix(0,18,32)
Ex_NM_95<-matrix(0,18,32)
Ex_NM_00<-matrix(0,18,32)
Ex_NM_05<-matrix(0,18,32)
Ex_NM_10<-matrix(0,18,32)

# Lleno las matrices de todos los estados year x year
for (i in 1:32) {
  Ex_NH_70[,i]<-  Ex_NH_70a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_70[,i]<-  Ex_NH_70a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_75[,i]<-  Ex_NH_75a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_75[,i]<-  Ex_NH_75a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_80[,i]<-  Ex_NH_80a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_80[,i]<-  Ex_NH_80a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_85[,i]<-  Ex_NH_85a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_85[,i]<-  Ex_NH_85a[(30+((i-1)*42)):(47+((i-1)*42)),]
 
  Ex_NH_90[,i]<-  Ex_NH_90a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_90[,i]<-  Ex_NH_90a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_95[,i]<-  Ex_NH_95a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_95[,i]<-  Ex_NH_95a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_00[,i]<-  Ex_NH_00a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_00[,i]<-  Ex_NH_00a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_05[,i]<-  Ex_NH_05a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_05[,i]<-  Ex_NH_05a[(30+((i-1)*42)):(47+((i-1)*42)),]
  
  Ex_NH_10[,i]<-  Ex_NH_10a[(10+((i-1)*42)):(27+((i-1)*42)),]
  Ex_NM_10[,i]<-  Ex_NH_10a[(30+((i-1)*42)):(47+((i-1)*42)),]
   
}

# Ahora si escribo la matriz final de la emigracion
Ex_NH_N<-matrix(0,18,9)
Ex_NM_N<-matrix(0,18,9)

# Lleno la matriz

Ex_NH_N[,1]<-rowSums(Ex_NH_70)
Ex_NH_N[,2]<-rowSums(Ex_NH_75)
Ex_NH_N[,3]<-rowSums(Ex_NH_80)
Ex_NH_N[,4]<-rowSums(Ex_NH_85)
Ex_NH_N[,5]<-rowSums(Ex_NH_90)
Ex_NH_N[,6]<-rowSums(Ex_NH_95)
Ex_NH_N[,7]<-rowSums(Ex_NH_00)
Ex_NH_N[,8]<-rowSums(Ex_NH_05)
Ex_NH_N[,9]<-rowSums(Ex_NH_10)

Ex_NM_N[,1]<-rowSums(Ex_NM_70)
Ex_NM_N[,2]<-rowSums(Ex_NM_75)
Ex_NM_N[,3]<-rowSums(Ex_NM_80)
Ex_NM_N[,4]<-rowSums(Ex_NM_85)
Ex_NM_N[,5]<-rowSums(Ex_NM_90)
Ex_NM_N[,6]<-rowSums(Ex_NM_95)
Ex_NM_N[,7]<-rowSums(Ex_NM_00)
Ex_NM_N[,8]<-rowSums(Ex_NM_05)
Ex_NM_N[,9]<-rowSums(Ex_NM_10)

# calcular las tasas- Son del orden correcto asi que estan bien :)
Ex_NH_T<-Ex_NH_N/APV_migNH
Ex_NM_T<-Ex_NM_N/APV_migNM
#matplot(Ex_NM_T,type = "l")
#matplot(Ex_NH_T,type = "l")

##################################
### Modelo de Rogers y Castro
##################################

# ahora necesito desagregar las tasas por edades 
# tengo Ex_IH_T Ex_IM_T- Emigracion internacional hombres y mujeres
# Ex_NH_T Ex_NM_T- Emigracion nacinal hombres y mujeres
# Ix_IH_T2 Ix_IM_T2- Inmigracion internacional hombres y mujeres
# Ix_NH_T2 Ix_NM_T2- Inmigracion nacinal hombres y mujeres

# Modelo de Rogers-Castro

# colnames(Ix_IM_T2)<-c(seq(1970,2010,5))
# # Internacional
# matplot((Ix_IM_T2)*1000,type = "l")
# legend(60, 2.1, legend = colnames(Ix_IM_T2), col = 1:9, lty = 1:9, bty="n")

# matplot(Ix_IH_T2,type = "l")
# matplot(Ex_IM_T2,type = "l")
# matplot(Ex_IH_T2,type = "l")
# 
# # Nacional (como la migracion interna se va a considerar constante consideraremos la del ultimo año)
# matplot(Ix_NM_T[,9],type = "l")
# matplot(Ix_NH_T[,9],type = "l")
# matplot(Ex_NM_T[,9],type = "l")
# matplot(Ex_NH_T[,9],type = "l")

# En Veracruz no se observan grandes curvas del retiro, por lo que modelaremos con el modelo de 7 parametros

# Fitting de la primera parte de la curva (0-14)
# Ejemplo de la emigracion nacional de mujeres

# Primera parte del modelo: migracion infantil
# solo me estoy tomando las jovenes - marcas de clase
#  a<- seq(2.5,12.5,5)
# # # selecciono los valores de y en las edades jovenes
#  y<-Ex_NM_T[1:3,9]
# # # Hago el modelo
# # # Modelo exponencial 
# #y.modelo <- lm(log(y)~(a))
# # #summary(y.modelo)
# # # aplico el modelo 
# # #y.ajustadas<-exp(y.modelo$coefficients[1]-y.modelo$coefficients[2]*b)
# # #plot(Ex_NM_T[1:3,9],type = "l")
# # #plot(y.ajustadas)
# # # De este modo las y.ajustadas son los valores de las tasas de los 0-14
# # 
# # # Creo un vector de todas las edaes en el grupo
#  b<-c(0:14)
#  q<-fitModel(y~A*exp(-B*a))
#  y.ajustadas<-coef(q)[1]*exp(-coef(q)[2]*b)
# # 
# # # Segunda parte del modelo. Curva unimodal sesgada a la izquierda
# # # Para la PEA
# # #Marcas de clase
#  a2<- seq(17.5,62.5,5)
# # # Tasas para la regresion
#  y2<-Ex_IM_T2[4:13,9]
#  plot(y2,type = "l")
# 
# nls_multstart(y2~(A)*exp(-B*(a2-C)-exp(-D*(a2-C))),
#              start_lower = c(A=0.02,B=15,C=0.05,D=0.08),
#              start_upper = c(A=0.19,B=38,C=0.33,D=1.49),
#              convergence_count = 100,
#              na.action = na.omit,
#              iter = 1)


 
# fitModel(y2~P*exp(-Q*(a2))+
#             (A)*exp(-B*(a2-C)-exp(-D*(a2-C))),start = list(P=0.02,Q=0.10,
#               A=0.07,B=20,C=0.15,D=0.42))
# # 
 # nls(y2~P*exp(-Q*(a2))+
 #       (A)*exp(-B*(a2-C)-exp(-D*(a2-C))),start = list(P=0.02,Q=0.10,
 #                                                      A=0.07,B=20,C=0.15,D=0.42))

# Mujeres
# nls_multstart(y2~P*exp(-Q*(a2))+(A)*exp(-B*(a2-C)-exp(-D*(a2-C))),
#               start_lower = c(P=0.005,Q=0.015,A=0.02,B=15,C=0.05,D=0.08),
#               start_upper = c(P=0.04,Q=0.41,A=0.19,B=38,C=0.33,D=1.49),
#               convergence_count = 100,
#               na.action = na.omit,
#               iter = 1)



#########
# plan B- Spline cubico para las tasas- De acuerdo a la IUSSP es valido hacer esto
#########

matplot(Ix_IM_T2,type = "l")
matplot(Ix_IH_T2,type = "l")
matplot(Ex_IM_T2,type = "l")
matplot(Ex_IH_T2,type = "l")

# Nacional (como la migracion interna se va a considerar constante consideraremos la del ultimo año)
matplot(Ix_NM_T[,9],type = "l")
matplot(Ix_NH_T[,9],type = "l")
matplot(Ex_NM_T[,9],type = "l")
matplot(Ex_NH_T[,9],type = "l")

# Para tener las tasas por edad Nacionales (constantes)
# marcas de clase
a2<-seq(2.5,87.5,5)
Ix_NM_T.fit<-cubicspline(a2,Ix_NM_T[,9],c(0:90))
Ix_NH_T.fit<-cubicspline(a2,Ix_NH_T[,9],c(0:90))
Ex_NM_T.fit<-cubicspline(a2,Ex_NM_T[,9],c(0:90))
Ex_NH_T.fit<-cubicspline(a2,Ex_NH_T[,9],c(0:90))

dev.off()
par(mfrow=c(2,2))
matplot(Ix_NM_T[,9],type = "l",ylab = "Tasas (Mujeres)",xlab = "Grupos de edad")
matplot(c(0:90),Ix_NM_T.fit,type = "l",ylab = "",xlab = "Edades simples")
matplot(Ix_NH_T[,9],type = "l",ylab = "Tasas (Hombres)",xlab = "Grupos de edad")
matplot(c(0:90),Ix_NH_T.fit,type = "l",ylab = "",xlab = "Edades simples")

dev.off()
par(mfrow=c(2,2))
matplot(Ex_NM_T[,9],type = "l",ylab = "Tasas (Mujeres)",xlab = "Grupos de edad")
matplot(c(0:90),Ex_NM_T.fit,type = "l",ylab = "",xlab = "Edades simples")
matplot(Ex_NH_T[,9],type = "l",ylab = "Tasas (Hombres)",xlab = "Grupos de edad")
matplot(c(0:90),Ex_NH_T.fit,type = "l",ylab = "",xlab = "Edades simples")

# Internacionales

Ix_IM_T.fit<-matrix(0,91,46)
Ix_IH_T.fit<-matrix(0,91,46)
Ex_IM_T.fit<-matrix(0,91,46)
Ex_IH_T.fit<-matrix(0,91,46)

for (i in 1:9) {
  Ix_IM_T.fit[,((5*(i-1)+1):(5*i))]<-cubicspline(a2,Ix_IM_T2[,i],c(0:90))
  Ix_IH_T.fit[,((5*(i-1)+1):(5*i))]<-cubicspline(a2,Ix_IH_T2[,i],c(0:90))
  Ex_IM_T.fit[,((5*(i-1)+1):(5*i))]<-cubicspline(a2,Ex_IM_T2[,i],c(0:90))
  Ex_IH_T.fit[,((5*(i-1)+1):(5*i))]<-cubicspline(a2,Ex_IH_T2[,i],c(0:90))
  Ix_IM_T.fit[,46]<-Ix_IM_T.fit[,45]
  Ix_IH_T.fit[,46]<-Ix_IH_T.fit[,45]
  Ex_IM_T.fit[,46]<-Ex_IM_T.fit[,45]
  Ex_IH_T.fit[,46]<-Ex_IH_T.fit[,45]
  
}

Ix_IH_T.fit[89:91,41:46]<-0.0000000000000000001

dev.off()
par(mfrow=c(2,2))
matplot(Ix_IM_T2,type = "l",ylab = "Tasas (Mujeres)",xlab = "Grupos de edad")
matplot(c(0:90),Ix_IM_T.fit,type = "l",ylab = "",xlab = "Edades simples")
matplot(Ix_IH_T2,type = "l",ylab = "Tasas (Hombres)",xlab = "Grupos de edad")
matplot(c(0:90),Ix_IH_T.fit,type = "l",ylab = "",xlab = "Edades simples")

dev.off()
par(mfrow=c(2,2))
matplot(Ex_IM_T2,type = "l",ylab = "Tasas (Mujeres)",xlab = "Grupos de edad")
matplot(c(0:90),Ex_IM_T.fit,type = "l",ylab = "",xlab = "Edades simples")
matplot(Ex_IH_T2,type = "l",ylab = "Tasas (Hombres)",xlab = "Grupos de edad")
matplot(c(0:90),Ex_IH_T.fit,type = "l",ylab = "",xlab = "Edades simples")

save(Ix_IM_T.fit,Ix_IH_T.fit,Ex_IM_T.fit,Ex_IH_T.fit,Ix_NM_T.fit,Ix_NH_T.fit,Ex_NM_T.fit,Ex_NH_T.fit,Px,Nx,mx2,fx1, file="Veracruz.RData")


