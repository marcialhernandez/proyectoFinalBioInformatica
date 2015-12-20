source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
#biocLite("made4")
#install.packages( pkgs= "gplots")
library(Biobase)
library(GEOquery)
library(genefilter)
#install.packages("ISLR")
#install.packages("rgdal")
#install.packages("pgirmess")
#install.packages("agricolae")
#install.packages("nortest")
library(ISLR)
library(made4)
library(ade4)
library(pgirmess)
library(agricolae)
require (nortest)
library(gplots)

############ Funciones ############ 

asignaColor <- function(valorClase) { 
  if (valorClase==1){
    return ("#FF0000") 
  }
  
  else {
    return ("#0000FF")
  }
}

cargaDatos<- function(archivo) { load(archivo); as.list(environment()) }

#Proceso de demostracion de normalidad

pruebaNormalidad<-function(filaNumerica){
  if (shapiro.test(filaNumerica)$p.value>0.1 & ad.test(filaNumerica)$p.value>0.1 & lillie.test(filaNumerica)$p.value>0.1){
    return (TRUE)
  }
  
  else{
    return (FALSE)
  }
}

pruebaHomocedasticidad<-function(filaLista){
  if (var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$BREAST))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$CNS))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$COLON))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$LEUK))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$MELAN))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$OVAR))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$PROSTATE))$p.value>0.1 & var.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$RENAL))$p.value>0.1){
    return (TRUE)
  }
  
  else{
    return (FALSE)
  }
}

#Pruebas t-student de a pares entre el cancer de pulmon y todos los otros tipos de cancer
#Para filas que hayan pasado las pruebas de normalidad y homocedasticidad
#El criterio a seguir es que si existen diferencias significativas con por lo menos
#dos tipos de cancer diferentes, se considera como un gen diferenciador del cancer de pulmon
pruebaDifSignificativasCancer<-function(filaLista,nombreGenActual){
  cantidadDiferenciasSignificativas<-0
  NSCLCBREAST<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$BREAST),var.equal = TRUE)$p.value
  NSCLCCNS<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$CNS),var.equal = TRUE)$p.value
  NSCLCCOLON<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$COLON),var.equal = TRUE)$p.value
  NSCLCLEUK<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$LEUK),var.equal = TRUE)$p.value
  NSCLCMELAN<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$MELAN),var.equal = TRUE)$p.value
  NSCLCOVAR<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$OVAR),var.equal = TRUE)$p.value
  NSCLCPROSTATE<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$PROSTATE),var.equal = TRUE)$p.value
  NSCLCRENAL<-t.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$RENAL),var.equal = TRUE)$p.value
  if (NSCLCBREAST<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCCNS<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCCOLON<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCLEUK<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCMELAN<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCOVAR<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCPROSTATE<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (NSCLCRENAL<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (cantidadDiferenciasSignificativas>=3){
    datosGenesParametricos<-paste("BREAST - NSCLC ",NSCLCBREAST,"\n","CNS - NSCLC ",NSCLCCNS,"\n","COLON - NSCLC ",NSCLCCOLON,"\n","LEUK - NSCLC ",NSCLCLEUK,"\n","MELAN - NSCLC ",NSCLCMELAN,"\n","NSCLC - OVAR ",NSCLCOVAR,"\n","NSCLC - PROSTATE ",NSCLCPROSTATE,"\n","NSCLC - RENAL ",NSCLCRENAL)
    logFile <- file(paste("genesDiferenciadoresParametricos/",nombreGenActual,".txt"))
    write(datosGenesParametricos,file=logFile,append=FALSE)
    close(logFile)
    return (TRUE)
  }
  else{
    return (FALSE)
  }
}

#Pruebas Wilcoxon Signed-Rank Test de a pares entre el cancer de pulmon y todos los otros tipos de cancer
#Para filas que no hayan pasado las pruebas de normalidad y homocedasticidad
#El criterio a seguir es que si existen diferencias significativas con por lo menos
#dos tipos de cancer diferentes, se considera como un gen diferenciador del cancer de pulmon
#Con atributo paired = FALSE, es equivalente a Mann-Whitney test
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
pruebaDifSignificativasCancerNonNormal<-function(filaLista){
  cantidadDiferenciasSignificativas<-0
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$BREAST),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$CNS),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$COLON),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$LEUK),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$MELAN),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$OVAR),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$PROSTATE),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (wilcox.test(as.numeric(filaLista$NSCLC),as.numeric(filaLista$RENAL),paired=FALSE)$p.value<0.01){
    cantidadDiferenciasSignificativas<-cantidadDiferenciasSignificativas+1
  }
  if (cantidadDiferenciasSignificativas>=3){
    return (TRUE)
  }
  else{
    return (FALSE)
  }
}

############ Main ############ 

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

# pValuesPorProbe<-c()
# cuentaAnormalidades<-0
# demostracionesHipotesisNormalidad<-list()
# for (probe in 1:dim(probes)[1]){
#   #pValuesPorProbe<-c(shapiro.test(as.numeric(probes[probe,]))$p.value,ad.test(as.numeric(probes[probe,]))$p.value,lillie.test(as.numeric(probes[probe,]))$p.value)
#   demostracionesHipotesisNormalidad[probe]<-pruebaNormalidad(as.numeric(probes[probe,]))
#   #for (valorHipotesis in 1:3){
#   #  if (pValuesPorProbe[valorHipotesis]<0.01){
#   #    demostracionesHipotesisNormalidad[probe]<-FALSE
#   #  }
#   #}
#   if (demostracionesHipotesisNormalidad[probe]==FALSE){
#     cuentaAnormalidades<-cuentaAnormalidades+1
#   }
# }
# if (cuentaAnormalidades>dim(probes)[2]/2){
#   print ("Los datos tienden a tener una distribucion anormal")
# }

#Proceso de diferenciacion multiclase

cantidadGenesARevisar<-cantidadOK
#cantidadGenesARevisar<-30
genesDiferenciadoresDistAnormal<-c()
cantidadGenesDistNormal<-0
cantidadGenesDistANormal<-0
posGenesDiferenciadores<-c()
genesDiferenciadoresDistNormal<-c()
demostracionesHipotesisNormalidad<-list()
#Para cada gen de la matriz probes (que corresponde a la matriz de genes ordenada por Pvalue)
for(filaGen in 1:cantidadGenesARevisar)
  {
  #Modo Lista
  tablaK<-list()
  #Modo Vector
  tablaVector<-c()
  contador<-1
  #Se crea una lista con la fila actual, para agrupar las distintas clases
  for (clase in nombresClasesSinRepeticion){
    tablaK[[nombresClasesSinRepeticion[contador]]]<-probes[filaGen,][matrizUbicacionesClases[[contador]]]
    tablaVector=c(tablaVector,probes[filaGen,][matrizUbicacionesClases[[contador]]])
    contador<-contador+1
  }
  #Se aplica Kruskall-Wallis al vector agrupado por clase
  #datosKruskalWallis<-kruskalmc(nombresClases,categ=nombresClases,data=tablaK,probs=0.1)
  
  if (pruebaNormalidad(as.numeric(probes[filaGen,]))==TRUE & pruebaHomocedasticidad(tablaK)==TRUE){
    print (paste("fila ",filaGen," - gen: ",rownames(probes)[filaGen]," cumple condiciones One Way ANOVA"))
    demostracionesHipotesisNormalidad[filaGen]<-TRUE
    #Proceso de diferenciacion multiclase usando t-students multiples
    cantidadGenesDistNormal<-cantidadGenesDistNormal+1
    if (pruebaDifSignificativasCancer(tablaK,rownames(probes)[filaGen])==TRUE){
      posGenesDiferenciadores<-c(posGenesDiferenciadores,filaGen)
      #Implica datos parametricos
      #Proceso de diferenciacion multiclase usando multiples t-student (Cancer de pulmon contra todos)
      genesDiferenciadoresDistNormal<-c(genesDiferenciadoresDistNormal,rownames(probes)[filaGen])
    }
  }
  
  else{ 
    #Implica datos no parametricos
    #Proceso de diferenciacion multiclase usando Kruskall-Wallis
    print (paste("fila ",filaGen," - gen: ",rownames(probes)[filaGen]," cumple condiciones Kruskall-Wallis"))
    demostracionesHipotesisNormalidad[filaGen]<-FALSE
    cantidadGenesDistANormal<-cantidadGenesDistANormal+1
    cantidadDiferenciasEncontradas<-0
    #Vector que contendra las filas de las comparaciones con el tipo de cancer de pulmon, en caso que sea significativo, se guardara
    #informeTemporal<-c(paste("gen: ",rownames(probes)[filaGen],"\n"))
    informeTemporal<-c()
    #datosKruskalWallis<-kruskalmc(nombresClases,categ=nombresClases,data=probes[filaGen,],probs=0.1) <<<<-- Metodo deshechado, pues no se puede cambiar el ajuste bonferroni a los pvalues
    datosKruskalWallis<-kruskal(as.numeric(tablaVector),nombresClases,alpha=0.01,p.adj="BH",group=FALSE)
    for (numeroFila in  1:dim(datosKruskalWallis$comparison)[1]){
      #Si la fila corresponde a la comparacion de NSCLC y ademas tiene un pValue < 0.01
      if (length(grep("NSCLC",rownames(datosKruskalWallis$comparison)[numeroFila]))>0) {
        #print (datosKruskalWallis$comparison[numeroFila,])
        #Se indica mediante un contador
        informeTemporal<-c(informeTemporal,paste(rownames(datosKruskalWallis$comparison)[numeroFila],datosKruskalWallis$comparison[numeroFila,]$pvalue),"\n")
        if (datosKruskalWallis$comparison[numeroFila,]$pvalue<0.01){
        cantidadDiferenciasEncontradas<-cantidadDiferenciasEncontradas+1
        }
      }
    }
    #print (cantidadDiferenciasEncontradas)
    #Si la cantidad de diferencias significativas son igual o mayores a dos
    if (cantidadDiferenciasEncontradas>=3){
      #Pues, el gen representa en parte
      genesDiferenciadoresDistAnormal<-c(genesDiferenciadoresDistAnormal,rownames(probes)[filaGen])
      posGenesDiferenciadores<-c(posGenesDiferenciadores,filaGen)
      logFile <- file(paste("genesDiferenciadoresNoParametricos/",rownames(probes)[filaGen],".txt"))
      write(paste(informeTemporal),file=logFile,append=FALSE)
      close(logFile)
    }
  }
  print ("__________________")
  
}

#Reasignacion de grupos, generacion de los grupos "cancer de pulmon", y cualquier otro

nombresClasesSimple<-c()
for(filaClase in 1:length(nombresClases)){
  if (nombresClases[filaClase]=="NSCLC"){
    nombresClasesSimple<-c(nombresClasesSimple,1)
  }
  else{
    nombresClasesSimple<-c(nombresClasesSimple,2)
  }
}

mapeadoColor<-unlist(lapply(nombresClasesSimple, asignaColor))

print ("GenesComportamientoParametricoDiferenciadores: ")
print (paste(genesDiferenciadoresDistNormal))
print ("GenesComportamientoNoParametricoDiferenciadores: ")
print (paste(genesDiferenciadoresDistAnormal))

print ("Creando ClusterizacionKMeans.pdf")
pdf("ClusterizacionKMeansInvestigacionPrincipal.pdf", paper="a4", width=8, height=8)
probes2<-probes[posGenesDiferenciadores,]
probes3<-t(probes2)[,nrow(probes2):1]
#probes2<-probes[1:15,]
clusterizacionKmedias <- kmeans (probes3,2)
grafico<-plot(probes3, col = clusterizacionKmedias$cluster, type="n")
text(probes3,labels=nombresClasesSimple,col=clusterizacionKmedias$cluster)
points(clusterizacionKmedias$centers, col = 1:2, pch = 16, cex = 2)
title(main="Grafico de K-medias en la investigacion principal")
dev.off()

print ("Creando EuclidianHeatMapInvestPrincipal")
pdf("EuclidianHeatMapInvestPrincipal.pdf", paper="a4", width=8, height=8)
#png("EuclidianHeatMapInvest1ConTratamiento.png")
euclidianHeatMap<-heatmap.2(data.matrix(probes2),col=greenred(50),xlab='',ylab='',
                            labRow=rownames(probes2),ColSideColors=mapeadoColor, scale="row", key=TRUE,
                            symkey=FALSE, density.info="none", trace="none", cexRow=0.5,
                            distfun = function(x) dist(x,method = 'euclidean'),labCol=colnames(probes2),
                            hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()