library(GEOquery)
library(oligo)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(limma)
library(affy)
library(dplyr)
library(knitr)
library(RColorBrewer)
library(hugene20sttranscriptcluster.db)

#getting the processed data
gse56481 <- getGEO('GSE56481', GSEMatrix=TRUE)
show(gse56481)
class(gse56481)
class(gse56481[[1]])
length(gse56481)
gse56481_expres<-gse56481[[1]]

#getting and reading the CEL files with comparison
getGEOSuppFiles('GSE56481')
untar("GSE56481/GSE56481_RAW.tar", exdir = 'Data/')
celpath = "C:\\Users\\neera\\Desktop\\R programming\\f\\Data"
list = list.files(celpath,full.names=TRUE)
gse56481_expres$supplementary_file
pd <- pData(gse56481_expres)
pd
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse56481_celdata<-read.celfiles(paste0('Data/',pd$cel_file),phenoData=phenoData(gse56481_expres))
show(gse56481_celdata)
image(gse56481_celdata[,1])
pData(gse56481_celdata)[,c("geo_accession","diagnosis:ch1","facs:ch1")]

#RMA and plotting CEL file densities 
hist(gse56481_celdata,transfo=log2,which=c("all"),main="CEL file densities before normalisation and background correction")

oligo_backgroundc<- backgroundCorrect(gse56481_celdata)
hist(oligo_backgroundc,transfo=log2,lwd=2,xlab='log intensity', which=c("all"),
     main="CEL file densities after background correction")

oligo_normalised <- normalize(oligo_backgroundc,method='quantile',which='all')
hist(oligo_normalised,transfo=log2,lwd=2,xlab='log intensity', which=c("all"),
     main="CEL file densities after quantile normalisation")

oligo_summarised <- oligo::rma(oligo_normalised,background=FALSE,normalize=FALSE)
hist(oligo_summarised,xlab='log intensity', main="CEL file densities after RMA without background corrrection and normalization")

gse56481_eset <- oligo::rma(gse56481_celdata)
hist(gse56481_eset,xlab='log intensity', main="CEL file densities after RMA")

#selecting the variables grouping
varLabels(gse56481_eset)
gse56481_eset$`diagnosis:ch1`
gse56481_eset$`facs:ch1`
pd <- pData(gse56481_eset)
pd <- rename(pd,diagnosis="diagnosis:ch1",facs="facs:ch1")
pd$diagnosis <- as.factor(pd$diagnosis)
levels(pd$diagnosis) <- c("GPA","Healthy.Control")
pd$group <- as.factor(paste(pd$facs,pd$diagnosis))
levels(pd$group) <- c("CD4.GPA","CD4.Healthy.Control","CD4.CD8.GPA","CD4.CD8.Healthy.Control","CD8.GPA","CD8.Healthy.Control")

#Designing a model
design <- model.matrix(~ 0 + pd$group)
colnames(design) <- levels(pd$group)
design

#making a contrast matrix
Table_data<-pd[c("diagnosis","facs","group")]
contrasts_matrix <- makeContrasts(CD4.D8 = CD4.CD8.GPA - CD4.CD8.Healthy.Control,
                                  CD4 = CD4.GPA - CD4.Healthy.Control,
                                  CD8 = CD8.GPA - CD8.Healthy.Control,
                                  Control = CD4.CD8.Healthy.Control - CD4.Healthy.Control- CD8.Healthy.Control,
                                  GPA = CD4.CD8.GPA - CD4.GPA - CD8.GPA,
                                  interaction=(CD4.CD8.Healthy.Control - CD4.Healthy.Control- CD8.Healthy.Control) - (CD4.CD8.GPA - CD4.GPA - CD8.GPA),
                                  levels=design)
contrasts_matrix
kable(contrasts_matrix)

#Empirical Bayes correction in limma
gse56481_fit <- lmFit(gse56481_eset,design)
gse56481_fit2 <- contrasts.fit(gse56481_fit,contrasts=contrasts_matrix)
gse56481_fit2 <- eBayes(gse56481_fit2)
tT<-topTable(gse56481_fit2)
summary(decideTests(gse56481_fit2,lfc=1))

#downstream analysis
volcanoplot(gse56481_fit2,coef=1,main=sprintf("Volcano plot"))
Genes_int <- topTable(gse56481_fit2,number=Inf,p.value = 0.05,lfc=10)
volcanoplot(gse56481_fit2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(Genes_int)), highlight=10)
points(Genes_int[['logFC']],-log10(Genes_int[['P.Value']]),col='red')
eset_of_interest <- gse56481_eset[rownames(Genes_int),]
heatmap(exprs(eset_of_interest))
heatmap(exprs(eset_of_interest),
        labCol=gse56481_eset[['Diagnosis']],
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
pheatmap(eset_of_interest[100:150])
pheatmap(eset_of_interest[250:300])

#annotation of genomic data the probe IDs
columns(hugene20sttranscriptcluster.db)
keytypes(hugene20sttranscriptcluster.db)
ps2 <- topTable(gse56481_fit2,number=Inf,p.value = 0.05,lfc=10)
ps2_up <- rownames(ps2)
df<-AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")                   
dplyr::mutate(df,GENENAME=stringr::str_trunc(GENENAME,30))
write.csv(df,file = "./Annotation.csv")