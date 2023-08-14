rm(list = ls())
library(phyloseq)
library(edgeR)
FunTable<-read.csv("zotus_tab_tax_final_normal.faprotax.csv")
rownames(FunTable)<-FunTable[,1]
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
	  
	  FunTableSelect<-data.frame(FunTable[,which(colnames(FunTable) %in% rownames(SampleData))])
      FunTableSelect =  round(data.frame(FunTableSelect*100000))
	  Treat<- factor(SampleData$Treatment)
      Time<- factor(substring(SampleData$Timepoint,4,4))
      y <- DGEList(counts=FunTableSelect, group=Treat)
      y <- calcNormFactors(y)
	  design <- model.matrix(~Time+Treat)#, data=targets,robust=TRUE)
      y <- estimateDisp(y, design)#, robust=TRUE)
      fit <- glmQLFit(y, design)#, robust=TRUE)
      qlf <- glmQLFTest(fit)
	  Result=qlf$table
      Result$fdr<-p.adjust(Result$PValue,method="BH")
      Result$Experiment<-rep(Exp,nrow(Result))
      Result$Compartment<-rep(Com,nrow(Result))
      Result$Variety<-rep(Var,nrow(Result))

      if(Exp == "Exp1"){
        Result$rel.TP1<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day1")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day7")$SampleID)][rownames(Result),])
      }else{
        Result$rel.TP1<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day2")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day4")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),]) 
      }
	  
      Result$rel.Intact<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Treatment == "Intact")$SampleID)][rownames(Result),])
      Result$rel.Defoliation<-rowMeans(FunTableSelect[,which(colnames(FunTableSelect) %in% subset(SampleData,Treatment == "Defoliation")$SampleID)][rownames(Result),])
	  Result$function_name = rownames(Result)
      ResultAll<-rbind(ResultAll,Result)
    }
  }
}

	  

write.csv(ResultAll,"OutputDir/edgeR_function_result_defoliation.csv")
	  
	  
	  
	  
	  
	  
	  
	   
	   