source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
#biocLite("made4")
library(Biobase)
library(GEOquery)
library(genefilter)
#install.packages("ISLR")
#install.packages("rgdal")
#install.packages("pgirmess")
library(ISLR)
library(made4)
library(ade4)
library(pgirmess)

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
nombresClasesSinRepeticion<-unique(nombresClases)
matrizUbicacionesClases<- c()
contador<-1
for (nombreClase in nombresClasesSinRepeticion){
  vectorTemporal <- which(colnames(dataComplete.affy)==nombreClase)
  #print(vectorTemporal)
  #matrizUbicacionesClases<-c(matrizUbicacionesClases,vectorTemporal)
  matrizUbicacionesClases[[contador]]<-vectorTemporal
  contador<-contador+1
}
#print (matrizUbicacionesClases)
#Proceso de ordenamiento multiclase usando Kruskall.test
pvalue<-0
cantidadOK<-0
for(filaGen in 1:cantidadGenes) {
  
  #mat_exp_t<-0
  tabla<-list()
  contador<-1
  for (clase in nombresClasesSinRepeticion){
    tabla[[nombresClasesSinRepeticion[contador]]]<-dataComplete.affy[filaGen,][matrizUbicacionesClases[[contador]]]
    #tabla<-setNames(c(dataComplete.affy[filaGen,][matrizUbicacionesClases[[contador]]]),nombresClases[contador])
    #mat_exp_t[[contador]]<-dataComplete.affy[filaGen,][matrizUbicacionesClases[[contador]]]
    contador<-contador+1
  }

  pvalue[filaGen]<-kruskal.test(tabla)$p.value
  if (pvalue[filaGen]<0.05){
    #print (pvalue[filaGen])
    cantidadOK<-cantidadOK+1
  }
}
print (paste("Cantidad de datos con pvalue < 0.05: ",cantidadOK))

#Proceso de ordenamiento por pValue
probes<-dataComplete.affy[order(pvalue),]

#Proceso de diferenciacion multiclase usando Kruskall-Wallis

cantidadGenesARevisar<-10
genesDiferenciadores<-c()

#Para cada gen de la matriz probes (que corresponde a la matriz de genes ordenada por Pvalue)
for(filaGen in 1:cantidadGenesARevisar)
  {
  #Modo Lista
  tablaK<-list()
  #Modo Vector
  #tablaK<-c()
  contador<-1
  #Se crea una lista con la fila actual, para agrupar las distintas clases
  for (clase in nombresClasesSinRepeticion){
    tablaK[[nombresClasesSinRepeticion[contador]]]<-probes[filaGen,][matrizUbicacionesClases[[contador]]]
    #tablaK=c(tablaK,probes[filaGen,][matrizUbicacionesClases[[contador]]])
    contador<-contador+1
  }
  cantidadDiferenciasEncontradas<-0
  #Se aplica Kruskall-Wallis al vector agrupado por clase
  datosKruskalWallis<-kruskalmc(nombresClases,categ=nombresClases,data=tablaK)
  for (numeroFila in  1:dim(datosKruskalWallis$dif.com)[1]){
    #Si la fila corresponde a la comparacion de CNS y comprueba la diferencia significativa mediante TRUE
    if (length(grep("CNS",rownames(datosKruskalWallis$dif.com)[numeroFila]))>0&datosKruskalWallis$dif.com$difference[numeroFila]==TRUE) {
      #Se indica mediante un contador
      cantidadDiferenciasEncontradas<-cantidadDiferenciasEncontradas+1
    }
  }
  #Si la cantidad de diferencias significativas son igual o mayores a dos
  if (cantidadDiferenciasEncontradas>=2){
    #Pues, el gen representa en parte
    genesDiferenciadores<-c(genesDiferenciadores,rownames(probes)[filaGen])
  }
}