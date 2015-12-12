source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
#biocLite("made4")
library(Biobase)
library(GEOquery)
library(genefilter)
#install.packages("ISLR")
#install.packages("pgirmess")
library(ISLR)
library(made4)
library(ade4)
#library(pgirmess)

cargaDatos<- function(archivo) { load(archivo); as.list(environment()) }

#Datos previo procesamiento de made4 <- Affy
dataLibreriaISLR<-assayData(NCI60)
#nci.labs<-NCI60$labs
#nci.data<-NCI60$data
#Metodos asociados
#assayData, featureNames, featureNames, notes, sampleNames, storageMode,genefinder 

#Datos parseados con made4
data(NCI60)
#print (names(NCI60))
dataISLRProcesada.affy<-NCI60$Affy
dataISLRProcesada.ross<-NCI60$Ross
EtiquetasGenes<-NCI60$Annot

#Datos completos obtenidos desde
#http://www.bioinf.ucd.ie/people/aedin/R/full_datasets/
dataComplete<-assayData(cargaDatos("NCI60.rda"))

#Tabla que indica la pertenencia de la clase por cada probe
#dataComplete.clases<-dataComplete$NCI60$classes
#Tabla de expresion de genes, etiquetada segun la identificacion de Affymetrix
dataComplete.affy<-dataComplete$NCI60$Affy
#NSCLC Corresponde al cancer de pulmon
cantidadProbes<-dim(dataComplete.affy)[2]
cantidadGenes<-dim(dataComplete.affy)[1]
nombresClases<-c()
for (nombreProbe in colnames(dataComplete.affy)){
  nombresClases<-c(nombresClases,strsplit(nombreProbe,"_")[[1]][1])
}
#Se ha renonmbrado las clases, quedando todos los grupos con nombres comunes
colnames(dataComplete.affy) <- nombresClases
#Se obtiene una lista sin repeticiones del nombre de las clases
nombresClases<-unique(nombresClases)
matrizUbicacionesClases=c()
for (nombreClase in nombresClases){
  vectorTemporal <- which(colnames(dataComplete.affy)==nombreClase)
  #print(vectorTemporal)
  matrizUbicacionesClases<-c(matrizUbicacionesClases,vectorTemporal)
}
#Proceso de diferenciacion multiclase usando Kruskall-Wallis
pvalue<-0

for(a in 1:cantidadGenes) {
  
  mat_exp_t<-c()
  for (posClases in  matrizUbicacionesClases){
    if (a==1){
    #print (dataComplete.affy[posClases])
    }
    c(mat_exp_t, dataComplete.affy[a,posClases])
    
  }
  tabla<-data.frame(cbind(mat_exp_t,nombresClases))
  pvalue[a]<-kruskal.test(tabla[,1]~nombresClases,data=tabla)$p.value
}

dataComplete.affyOrdenada<-dataComplete.affy[order(pvalue),]
#dataComplete.affyOrdenada<-dataComplete.affy[order(pvalue)]
#pvalue<-pvalue[order(pvalue)]

#print(length(which(pvalue<=0.05)))