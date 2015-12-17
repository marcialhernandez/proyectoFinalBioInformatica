source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("qvalue")
#biocLite("genefilter")
#install.packages( pkgs= "gplots")
#biocLite("pgirmess")
#install.packages("stringr", dependencies=TRUE)
print ("Precargando Librerias")
library(Biobase)
library(GEOquery)
library(genefilter)
library(gplots)
library(gtools)
library(qvalue)
library(pgirmess)
require (nortest)
library(stringr)

#Investigacion tomada de:
#http://www.ncbi.nlm.nih.gov/pubmed/15867382

#Titulo: Lung cancer cell line response to motexafin gadolinium: time course

#Sumario: Expression profiling of lung cancer cell line A549 following treatment 
#with 50 uM of the metal cation-containing chemotherapeutic drug motexafin gadolinium (MGd). 
#Cells examined at 4, 12, and 24 hours following treatment. Results provide insight into 
#the mechanism of action of MGd. 

options('download.file.method.GEOquery'='auto')

asignaColorSinTratamiento <- function(valorClase) { 
  if (valorClase=="control4h"){
    return ("#FF0000") 
  }
  
  else if (valorClase=="control12h"){
    return ("#FFA500")
  }
  else {
    return ("#0000FF")
  }
}

asignaColorConTratamiento <- function(valorClase) { 
  if (valorClase=="motexafingadolinium4h"){
    return ("#FF0000") 
  }
  
  else if (valorClase=="motexafingadolinium12h"){
    return ("#FFA500")
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

pruebaHomocedasticidad<-function(filaGrupo1,filaGrupo2,filaGrupo3){
  if (var.test(filaGrupo1,filaGrupo2)$p.value>0.1 &var.test(filaGrupo1,filaGrupo3)$p.value>0.1&var.test(filaGrupo2,filaGrupo3)$p.value>0.1){
    return (TRUE)
  }
  
  else{
    return (FALSE)
  }
}

#Se aplica un ajuste parametrico usando metodo parametrico
testParametricoMultiplesGrupos<-function(grupo1,grupo2,grupo3){
  return (p.adjust(as.numeric(c(t.test(grupo1,grupo2)$p.value,t.test(grupo1,grupo3)$p.value,t.test(grupo2,grupo3)$p.value)),"bonferroni"))
}

#Se aplica un ajute a los pvalues no parametrico usando metodo no parametrico
#Con atributo paired = FALSE, es equivalente a Mann-Whitney test
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
testNoParametricoMultiplesGrupos<-function(grupo1,grupo2,grupo3){
  return (p.adjust(as.numeric(c(wilcox.test(grupo1, grupo2 ,paired=FALSE)$p.value,wilcox.test(grupo1, grupo3 ,paired=FALSE)$p.value,wilcox.test(grupo2, grupo3 ,paired=FALSE)$p.value)),"BH"))
}


print ("Carga de datos Soft a variable global")
datos <- getGEO(filename='GDS1204.soft.gz')

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

featureDataEset<-featureData(eset)

colnames(mat_exp)<-sampleNames(phenoData(eset))
rownames(mat_exp)<-sampleNames(featureDataEset)

print ("[Validacion] Imprime una porcion de la matriz de expresion")
print(mat_exp[1:15,1:5])

#Se identificaron 2 clases, aquellos con tratamiento y aquellos sin tratamiento
#Pero a su vez, estan separados por horas, por 4-12 y 24 horas con y sin tratamiento
#Por ello se combinan ambos atributos, obteniendo 6 clases
clases<-str_replace_all(paste0(pData(phenoData(eset))$agent,pData(phenoData(eset))$time), fixed(" "), "")
#Y con esto se obtienen valores numericos que representan a cada clase

c1 <- which(clases=="control4h")
c2 <- which(clases=="control12h")
c3 <- which(clases=="control24h")
c4 <- which(clases=="motexafingadolinium4h")
c5 <- which(clases=="motexafingadolinium12h")
c6 <- which(clases=="motexafingadolinium24h")

pvalueSinTratamiento<-0
pvalueConTratamiento<-0

mat_exp2<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)
mat_exp3<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)

#Analisis de grupos
for(i in 1:cantidadGenes){
  mat_exp2[i,]=c(mat_exp[i,c1],mat_exp[i,c2],mat_exp[i,c3]) #Expresion de genes del grupo sin tratamiento
  mat_exp3[i,]=c(mat_exp[i,c4],mat_exp[i,c5],mat_exp[i,c6]) #Expresion de genes del grupo con tratamiento
  
  ###################################### Analisis de expresion sin tratamiento  ###################################### 
  if (pruebaNormalidad(mat_exp2[i,])==TRUE & pruebaHomocedasticidad(mat_exp[i,c1],mat_exp[i,c2],mat_exp[i,c3])==TRUE){
    
    #En caso que sea parametrica
    pvalueSinTratamiento[i]<-testParametricoMultiplesGrupos(mat_exp[i,c1],mat_exp[i,c2],mat_exp[i,c3])
  }
  
  #En caso que sea no-parametrica
  else{

    pvalueSinTratamiento[i]<-testNoParametricoMultiplesGrupos(mat_exp[i,c1],mat_exp[i,c2],mat_exp[i,c3])
  }

  ###################################### Analisis de expresion con tratamiento  ###################################### 
  
  if (pruebaNormalidad(mat_exp3[i,])==TRUE & pruebaHomocedasticidad(mat_exp[i,c4],mat_exp[i,c5],mat_exp[i,c6])==TRUE){
    
    #En caso que sea parametrica
    pvalueConTratamiento[i]<-testParametricoMultiplesGrupos(mat_exp[i,c4],mat_exp[i,c5],mat_exp[i,c6])
  }
  
  #En caso que sea no-parametrica
  else{
    
    pvalueConTratamiento[i]<-testNoParametricoMultiplesGrupos(mat_exp[i,c4],mat_exp[i,c5],mat_exp[i,c6])
  }
}

#Ordenar matriz
mat_exp2<-mat_exp[order(pvalueSinTratamiento),]
mat_exp3<-mat_exp[order(pvalueConTratamiento),]
pvalueSinTratamiento<-pvalueSinTratamiento[order(pvalueSinTratamiento)]
pvalueConTratamiento<-pvalueConTratamiento[order(pvalueConTratamiento)]

ExpresionInvestigacion1SinTratamiento <- data.frame(rownames(mat_exp2),pvalueSinTratamiento)
ExpresionInvestigacion1ConTratamiento <- data.frame(rownames(mat_exp3),pvalueConTratamiento)

print ("[Validacion] Los 15 genes mas expresados del conjunto de datos sin tratamiento son: ")
print(ExpresionInvestigacion1SinTratamiento[1:15,1:2])

print ("[Validacion] Los 15 genes mas expresados del conjunto de datos con tratamiento son: ")
print(ExpresionInvestigacion1ConTratamiento[1:15,1:2])

logFile1 <- file("InvestigacionRelacionada1-GDS1204/GenesExpresadosSinTratamiento.txt")
write.table(ExpresionInvestigacion1SinTratamiento[1:200,1:2],file=logFile1,append=FALSE,sep="\t")
#close(logFile1)

logFile2 <- file("InvestigacionRelacionada1-GDS1204/GenesExpresadosConTratamiento.txt")
write.table(ExpresionInvestigacion1ConTratamiento[1:200,1:2],file=logFile2,append=FALSE,sep="\t")
#close(logFile2)