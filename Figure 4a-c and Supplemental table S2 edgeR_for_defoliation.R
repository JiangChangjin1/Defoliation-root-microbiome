

rm(list = ls())
library(phyloseq)
library(edgeR)
IntData = readRDS("IntegratedBacteriaData.rds")
ResultAll<-c()
for(Com in c("Root","Rhizosphere")){
  for(Exp in c("Exp1","Exp2")){
      for(Var in c("NIP","MH63")){
      IntDataC<- subset_samples(IntData, Compartment %in% Com)
      IntDataCE<- subset_samples(IntDataC, Experiment %in% Exp)
      IntDataCEV<- subset_samples(IntDataCE, Variety %in% Var )
	  
      IntDataCEV_Sel = filter_taxa(IntDataCEV,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
      SampleData<-data.frame(sample_data(IntDataCEV_Sel))
      SampleData$Treatment <- factor(SampleData$Treatment,levels=c("Intact","Defoliation"))
      SampleData$Timepoint <- factor(SampleData$Timepoint,unique(SampleData$Timepoint))
      OtuTable<-data.frame(otu_table(IntDataCEV_Sel))
	  
      Treat<- factor(SampleData$Treatment)
      Time<- factor(substring(SampleData$Timepoint,4,4))
      y <- DGEList(counts=OtuTable, group=Treat)
      y <- calcNormFactors(y)
      design <- model.matrix(~Time+Treat)#, data=targets,robust=TRUE)
      y <- estimateDisp(y, design)#, robust=TRUE)
      fit <- glmQLFit(y, design)#, robust=TRUE)
      qlf <- glmQLFTest(fit)
      Result=qlf$table
      Result$fdr<-p.adjust(Result$PValue,method="BH")
	  
      TaxTable<-data.frame(tax_table(IntDataCEV))
      keep<-which(rownames(TaxTable) %in% rownames(Result))
      TaxTable<-TaxTable[keep,]
      ResultAddTax = merge(Result,TaxTable, by = "row.names")
      rownames(ResultAddTax)<-ResultAddTax$Row.names
      ResultAddTax$Experiment<-rep(Exp,nrow(TaxTable))
      ResultAddTax$Compartment<-rep(Com,nrow(TaxTable))
      ResultAddTax$Variety<-rep(Var,nrow(TaxTable))
      IntDataCEV_SelAbun <- transform_sample_counts(IntDataCEV_Sel, function(x) x / sum(x) )
      OtuTableAbun<-data.frame(otu_table(IntDataCEV_SelAbun))
	  
      if(Exp == "Exp1"){
        ResultAddTax$rel.TP1<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day1")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP2<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP3<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP4<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day7")$SampleID)][rownames(ResultAddTax),])
      }else{
        ResultAddTax$rel.TP1<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day2")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP2<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP3<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day4")$SampleID)][rownames(ResultAddTax),])
        ResultAddTax$rel.TP4<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(ResultAddTax),]) 
      }
      ResultAddTax$rel.Intact<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Treatment == "Intact")$SampleID)][rownames(ResultAddTax),])
      ResultAddTax$rel.Defoliation<-rowMeans(OtuTableAbun[,which(colnames(OtuTableAbun) %in% subset(SampleData,Treatment == "Defoliation")$SampleID)][rownames(ResultAddTax),])
      ResultAll<-rbind(ResultAll,ResultAddTax)
    }
  }
}


write.csv(ResultAll,"OutputDir/edgeR_result_defoliation.csv")


