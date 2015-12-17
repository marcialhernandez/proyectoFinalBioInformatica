source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("qvalue")
#biocLite("genefilter")
#install.packages( pkgs= "gplots")
#biocLite("pgirmess")
print ("Precargando Librerias")
library(Biobase)
library(GEOquery)
library(genefilter)
library(gplots)
library(gtools)
library(qvalue)
library(pgirmess)
require (nortest)

#Investigacion tomada de:
#http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2489

#Titulo: Large airway epithelium response to cigarette smoking (HuGeneFL)
#Sumario: Analysis of large airway epithelial cells of phenotypically normal smokers. 
#Lung cancer tumors exhibit neuroendocrine properties, and chronic smokers have 
#increased numbers of neuroendocrine cells. Results provide insight into the effect of 
#cigarette smoking on neuroendocrine cells.

options('download.file.method.GEOquery'='auto')

asignaColor <- function(valorClase) { 
  if (valorClase=="2"){
    return ("#FF0000") 
  }
  else {
    return ("#0000FF")
  }
}

pruebaNormalidad<-function(filaNumerica){
  if (shapiro.test(filaNumerica)$p.value>0.1 & ad.test(filaNumerica)$p.value>0.1 & lillie.test(filaNumerica)$p.value>0.1){
    return (TRUE)
  }
  
  else{
    return (FALSE)
  }
}

pruebaHomocedasticidad<-function(filaGrupo1,filaGrupo2){
  if (var.test(filaGrupo1,filaGrupo2)$p.value>0.1 ){
    return (TRUE)
  }
  
  else{
    return (FALSE)
  }
}

print ("Carga de datos Soft a variable global")
datos <- getGEO(filename='GDS2489.soft.gz')

#Impresion de caracteristicas de los datos por salida estandar
print (Meta(datos)$description[1])
print (Meta(datos)$feature_count)
print (Meta(datos)$sample_type)

cantidadGenes<-strtoi(Meta(datos)$feature_count)
cantidadMuestras<-strtoi(Meta(datos)$sample_count)

print ("Conversion de datos a objeto tipo ExpressionSet")
eset <- GDS2eSet(datos, do.log2=TRUE,GPL=NULL,AnnotGPL=FALSE)

mat_exp<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)

for(a in 1:cantidadGenes)
{
  for(b in 1:cantidadMuestras)
  {
    mat_exp[a,b]<-exprs(eset)[a,b]
  }
}
pvalue<-0
featureDataEset<-featureData(eset)

colnames(mat_exp)<-sampleNames(phenoData(eset))
rownames(mat_exp)<-sampleNames(featureDataEset)

print ("[Validacion] Imprime una porcion de la matriz de expresion")
print(mat_exp[1:50,1:5])

#2 <- control
#1 <- fumadores
clases<-c(pData(phenoData(eset))$stress)
cantidadControl<-18
cantidadFumadores<-26
c1 <- which(clases==1) #Lista de ubicaciones de probes de fumadores 
c2 <- which(clases==2) #Lista de ubicaciones de probes de control 

mat_exp2<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)
colnames(mat_exp2)<-sampleNames(phenoData(eset))
rownames(mat_exp2)<-sampleNames(featureDataEset)

for(i in 1:cantidadGenes){
  mat_exp2[i,]=c(mat_exp[i,c2],mat_exp[i,c1])
  
  if (pruebaNormalidad(mat_exp2[i,])==TRUE & pruebaHomocedasticidad(mat_exp[i,c2],mat_exp[i,c1])==TRUE){
    
    #En caso que sea parametrica
    pvalue[i]<-t.test(mat_exp[i,c2], mat_exp[i,c1])$p.value
  }
  
  #En caso que sea no-parametrica
  else{
    #Con atributo paired = FALSE, es equivalente a Mann-Whitney test
    #https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
    pvalue[i]<-wilcox.test(mat_exp[i,c2], mat_exp[i,c1],paired=FALSE)$p.value
  }
}

mapeadoColor<-unlist(lapply(clases, asignaColor))

#Ordenar matriz
mat_exp2<-mat_exp[order(pvalue),]
pvalue<-pvalue[order(pvalue)]


ExpresionInvestigacion2 <- data.frame(rownames(mat_exp2),pvalue)
colnames(ExpresionInvestigacion2)<-c("Gen","pValue")

print ("[Validacion] Los 15 genes mas expresados son: ")
print(ExpresionInvestigacion2[1:15,1:2])

logFile3 <- file("InvestigacionRelacionada2-GDS2489/GenesExpresados.txt")
write.table(ExpresionInvestigacion2[1:200,1:2],file=logFile3,append=FALSE,sep="\t")
close(logFile3)

print ("Creando EuclidianHeatMapInvest2.pdf")
pdf("EuclidianHeatMapInvest2.pdf", paper="a4", width=8, height=8)
euclidianHeatMap<-heatmap.2(mat_exp2[1:30,],col=greenred(50),xlab='',ylab='',
                                         labRow=rownames(mat_exp2[1:30,]),ColSideColors=mapeadoColor, scale="row", key=TRUE,
                                         symkey=FALSE, density.info="none", trace="none", cexRow=0.5,
                                         distfun = function(x) dist(x,method = 'euclidean'),labCol=colnames(mat_exp2[1:30,]),
                                         hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()