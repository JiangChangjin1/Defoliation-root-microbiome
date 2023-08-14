
rm(list = ls())
library(phyloseq)
library(Hmisc)
library(edgeR)
IntData = readRDS("IntegratedBacteriaData.rds")     
for(Tre in c("Intact","Defoliation")){
    for(Var in c("MH63","NIP")){
       for(Com in c("Root","Rhizosphere")){
           IntDataT =subset_samples(IntData, Treatment %in% Tre)
		   IntDataTV =subset_samples(IntDataT, Variety %in% Var)
           IntDataTVC<- subset_samples(IntDataTV, Compartment %in% Com)
           IntDataTVC_Rel <- transform_sample_counts(IntDataTVC, function(x) x / sum(x) )
           IntDataTVC_RelHigh = filter_taxa(IntDataTVC_Rel, function(x) mean(x) > 0.0001, TRUE)
           IntDataTVC_RelHigh = filter_taxa(IntDataTVC_RelHigh,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
           SampleData<-sample_data(IntDataTVC)
           OtuTable<-data.frame(otu_table(IntDataTVC))
           y <- DGEList(counts=OtuTable)
           y <- calcNormFactors(y)
           OtuTableCPM<- cpm(y)
           OtuTableCPM_Select<-OtuTableCPM[rownames(otu_table(IntDataTVC_RelHigh)),]
           comm.data<-data.frame(t(data.frame(OtuTableCPM_Select)))
           occor<-rcorr(as.matrix(comm.data),type="spearman")
           r.cor = occor$r 
           p.cor = occor$P 
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
           TaxTable<-data.frame(tax_table(IntDataTVC))
           TaxTable$label<-rownames(TaxTable)
           NodesTableAddTax<- merge(NodesTable,TaxTable, by.x = "label")
           OtuTableSel<-data.frame(otu_table(IntDataTVC))
           OtuAbunMatrix<-data.frame(name=rownames(OtuTableSel),abun=rowMeans(OtuTableSel))
		   
           OtuAbunMatrix$label<-rownames(OtuAbunMatrix)
           NodesTableAddTaxAbun<- merge(NodesTableAddTax,OtuAbunMatrix, by.x = "label")
           write.csv(EdgeTable,paste("OutputDir/",Tre,Var,Com,"network_edge_test",".csv",sep=""),row.names=FALSE)
           write.csv(NodesTableAddTaxAbun,paste("OutputDir/Supplemental Figure S9/",Tre,Var,Com,"network_nodes",".csv",sep=""),row.names=FALSE)      
           }
        }
    }
