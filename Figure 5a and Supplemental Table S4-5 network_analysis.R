
rm(list = ls())
##Root
library(edgeR)
library(phyloseq)
library(Hmisc)
DiffZotuTable = read.csv("edgeR_result_defoliation.csv")
DiffZotuTable = subset(DiffZotuTable,PValue<0.001)
DiffZotuTable = subset(DiffZotuTable,abs(logFC)>log2(1.5))
DiffZotuTableRe = subset(DiffZotuTable,Compartment == "Root")

IntData = readRDS("IntegratedBacteriaData.rds")
IntData = subset_samples(IntData, Treatment %in% c("Intact"))
IntData = subset_samples(IntData, Compartment %in% c("Root"))

SampleData<-sample_data(IntData)
OtuTable<-data.frame(otu_table(IntData))
Treat<- factor(SampleData$Treatment)
Time<- factor(substring(SampleData$Timepoint,4,4))
y <- DGEList(counts=OtuTable, group=Treat)
y <- calcNormFactors(y)
OtuTableCpm<- cpm(y)

IntDataSel = transform_sample_counts(IntData, function(x) x / sum(x) )
IntDataSel = filter_taxa(IntDataSel, function(x) mean(x) > 0.0001, TRUE)
IntDataSel = filter_taxa(IntDataSel,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
OtuTableCpm<-OtuTableCpm[rownames(otu_table(IntDataSel)),]

NetworkData<-data.frame(t(data.frame(OtuTableCpm)))
occor<-rcorr(as.matrix(NetworkData),type="spearman")
r.cor = occor$r # 取相关性矩阵R值
p.cor = occor$P # 取相关性矩阵p值
p.adj <- p.adjust(p.cor, method="BH")
r.matrix<-r.cor
p.matrix<-p.cor

for(i in 1:nrow(r.matrix)){
    r.matrix[i,i:ncol(r.matrix)]<-rep(0,length( r.matrix[i,i:ncol(r.matrix)]))
    }
	
EdgeTable<-data.frame()
for(i in 1:nrow(r.matrix)){
    if(length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])>0){
       tab<-data.frame(Source=rep(rownames(r.matrix)[i],length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])),Target=colnames(r.matrix)[which(abs(r.matrix[i,])>0)],r=r.matrix[i,][which(abs(r.matrix[i,])>0)],p=p.matrix[i,][which(abs(r.matrix[i,])>0)])
       rownames(tab)<-paste(tab$Source,tab$Target,sep="_")
       EdgeTable<-rbind(EdgeTable,tab)
       }
    }
	
EdgeTable<-EdgeTable[which(EdgeTable$r != 1),]
EdgeTable$adjust.p<-p.adjust(EdgeTable$p, method="BH")
EdgeTable<-subset(EdgeTable,abs(r)>0.8)
EdgeTable<-subset(EdgeTable,adjust.p<0.001)

name<-unique(c(as.character(EdgeTable$Source),as.character(EdgeTable$Target)))
NodesTable<-data.frame(id=name,label=name)
NodesTable$labelDiff<-"other"
NodesTable[which(NodesTable$label %in% DiffZotuTable$Row.names),]$labelDiff<-rep("diff",length(NodesTable[which(NodesTable$label %in% DiffZotuTable$Row.names),]$labelDiff))

TaxTable<-data.frame(tax_table(IntDataSel))
TaxTable$label<-rownames(TaxTable)
NodesTableAddTax<- merge(NodesTable,TaxTable, by.x = "label")
OtuTableSel<-data.frame(otu_table(IntDataSel))
OtuAbunMatrix<-data.frame(name=rownames(OtuTableSel),abun=rowMeans(OtuTableSel))
OtuAbunMatrix$label<-rownames(OtuAbunMatrix)
NodesTableAddTaxAbun<- merge(NodesTableAddTax,OtuAbunMatrix, by.x = "label")

write.csv(EdgeTable,paste("OutputDir/Edges_Table_Endo_network",".csv",sep=""),row.names=FALSE)
write.csv(NodesTableAddTaxAbun,paste("OutputDir/Nodes_Table_Endo_network",".csv",sep=""),row.names=FALSE)





        
##Rhizosphere
library(edgeR)
library(phyloseq)
library(Hmisc)
DiffZotuTable = read.csv("edgeR_result_defoliation.csv")
DiffZotuTable = subset(DiffZotuTable,PValue<0.001)
DiffZotuTable = subset(DiffZotuTable,abs(logFC)>log2(1.5))
DiffZotuTableRe = subset(DiffZotuTable,Compartment == "Rhizosphere")

IntData = readRDS("IntegratedBacteriaData.rds")
IntData = subset_samples(IntData, Treatment %in% c("Intact"))
IntData = subset_samples(IntData, Compartment %in% c("Rhizosphere"))

SampleData<-sample_data(IntData)
OtuTable<-data.frame(otu_table(IntData))
Treat<- factor(SampleData$Treatment)
Time<- factor(substring(SampleData$Timepoint,4,4))
y <- DGEList(counts=OtuTable, group=Treat)
y <- calcNormFactors(y)
OtuTableCpm<- cpm(y)

IntDataSel = transform_sample_counts(IntData, function(x) x / sum(x) )
IntDataSel = filter_taxa(IntDataSel, function(x) mean(x) > 0.0001, TRUE)
IntDataSel = filter_taxa(IntDataSel,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
OtuTableCpm<-OtuTableCpm[rownames(otu_table(IntDataSel)),]

NetworkData<-data.frame(t(data.frame(OtuTableCpm)))
occor<-rcorr(as.matrix(NetworkData),type="spearman")
r.cor = occor$r # 取相关性矩阵R值
p.cor = occor$P # 取相关性矩阵p值
p.adj <- p.adjust(p.cor, method="BH")
r.matrix<-r.cor
p.matrix<-p.cor

for(i in 1:nrow(r.matrix)){
r.matrix[i,i:ncol(r.matrix)]<-rep(0,length( r.matrix[i,i:ncol(r.matrix)]))
}

EdgeTable<-data.frame()
for(i in 1:nrow(r.matrix)){
    if(length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])>0){
       tab<-data.frame(Source=rep(rownames(r.matrix)[i],length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])),Target=colnames(r.matrix)[which(abs(r.matrix[i,])>0)],r=r.matrix[i,][which(abs(r.matrix[i,])>0)],p=p.matrix[i,][which(abs(r.matrix[i,])>0)])
       rownames(tab)<-paste(tab$Source,tab$Target,sep="_")
       EdgeTable<-rbind(EdgeTable,tab)
       }
    }
	
EdgeTable<-EdgeTable[which(EdgeTable$r != 1),]
EdgeTable$adjust.p<-p.adjust(EdgeTable$p, method="BH")
EdgeTable<-subset(EdgeTable,abs(r)>0.8)
EdgeTable<-subset(EdgeTable,adjust.p<0.001)

name<-unique(c(as.character(EdgeTable$Source),as.character(EdgeTable$Target)))
NodesTable<-data.frame(id=name,label=name)
NodesTable$labelDiff<-"other"
NodesTable[which(NodesTable$label %in% DiffZotuTable$Row.names),]$labelDiff<-rep("diff",length(NodesTable[which(NodesTable$label %in% DiffZotuTable$Row.names),]$labelDiff))

TaxTable<-data.frame(tax_table(IntDataSel))
TaxTable$label<-rownames(TaxTable)
NodesTableAddTax<- merge(NodesTable,TaxTable, by.x = "label")
OtuTableSel<-data.frame(otu_table(IntDataSel))
OtuAbunMatrix<-data.frame(name=rownames(OtuTableSel),abun=rowMeans(OtuTableSel))
OtuAbunMatrix$label<-rownames(OtuAbunMatrix)
NodesTableAddTaxAbun<- merge(NodesTableAddTax,OtuAbunMatrix, by.x = "label")

write.csv(EdgeTable,paste("OutputDir/Edges_Table_Rhizo_network",".csv",sep=""),row.names=FALSE)
write.csv(NodesTableAddTaxAbun,paste("OutputDir/Nodes_Table_Rhizo_network",".csv",sep=""),row.names=FALSE)
        
	