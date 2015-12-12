source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
#biocLite("made4")
library(Biobase)
library(GEOquery)
library(genefilter)
#install.packages("ISLR")
library(ISLR)
library(made4)
library(ade4)

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
dataComplete<-cargaDatos("NCI60.rda")

#Tabla que indica la pertenencia de la clase por cada probe
dataComplete.clases<-dataComplete$NCI60$classes
#Tabla de expresion de genes, etiquetada segun la identificacion de Affymetrix
dataComplete.affy<-dataComplete$NCI60$Affy