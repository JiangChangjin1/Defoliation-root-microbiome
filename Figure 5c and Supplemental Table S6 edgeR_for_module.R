
rm(list = ls())
set.seed(5)
library(phyloseq)
library(edgeR)
IntData = readRDS("IntegratedBacteriaData.rds")
ResultAll<-c()
for(Com in c("Rhizosphere")){
  for(Exp in c("Exp1","Exp2")){
      for(Var in c("NIP","MH63")){
      IntDataC<- subset_samples(IntData, Compartment %in% Com)
      IntDataCE<- subset_samples(IntDataC, Experiment %in% Exp)
      IntDataCEV<- subset_samples(IntDataCE, Variety %in% Var )
	  
      IntDataCEV_Sel = filter_taxa(IntDataCEV,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
      SampleData<-data.frame(sample_data(IntDataCEV_Sel))
      SampleData$Treatment <- factor(SampleData$Treatment,levels=c("Intact","Defoliation"))
      SampleData$Timepoint <- factor(SampleData$Timepoint,unique(SampleData$Timepoint))
      Otutable<-data.frame(otu_table(IntDataCEV_Sel))
	  
	  Nodes<-read.csv("rs_int_0.8_0.001_new_diff_nodes2.csv")
	  Moduletable<-c()
	  for(m in unique(Nodes$modularity_class)){
      NodesSel<-subset(Nodes,modularity_class == m)
	  OtutableSel<-Otutable[which(rownames(Otutable) %in% NodesSel$Id),]
	  Moduletable<-rbind(Moduletable,colSums(OtutableSel))
	  }
	  colnames(Moduletable)<-colnames(Otutable)
	  rownames(Moduletable)<-paste("module",unique(Nodes$modularity_class),sep="")
	  ModuleOther<-colSums(Otutable)-colSums(Moduletable)
	  ModuletableAll<-rbind(ModuleOther,Moduletable)
	  rownames(ModuletableAll)[1]<-"module_other"
	  
      Treat<- factor(SampleData$Treatment)
      Time<- factor(substring(SampleData$Timepoint,4,4))
      y <- DGEList(counts=ModuletableAll, group=Treat)
      y <- calcNormFactors(y)
      design <- model.matrix(~Time+Treat)
      y <- estimateDisp(y, design)
      fit <- glmQLFit(y, design)
      qlf <- glmQLFTest(fit)
      Result=qlf$table
      Result$fdr<-p.adjust(Result$PValue,method="BH")
      Result$Experiment<-rep(Exp,nrow(Result))
      Result$Compartment<-rep(Com,nrow(Result))
      Result$Variety<-rep(Var,nrow(Result))
	  
      if(Exp == "Exp1"){
        Result$rel.TP1<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day1")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day7")$SampleID)][rownames(Result),])
      }else{
        Result$rel.TP1<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day2")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day4")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),]) 
      }
	  
      Result$rel.Intact<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Treatment == "Intact")$SampleID)][rownames(Result),])
      Result$rel.Defoliation<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Treatment == "Defoliation")$SampleID)][rownames(Result),])
	  Result$module_name = rownames(Result)
      ResultAll<-rbind(ResultAll,Result)
    }
  }
}

	
write.csv(ResultAll,"OutputDir/edgeR_module_result_defoliation_rhizo.csv")




		
	
rm(list = ls())
set.seed(5)
library(phyloseq)
library(edgeR)
IntData = readRDS("IntegratedBacteriaData.rds")
ResultAll<-c()
for(Com in c("Root")){
  for(Exp in c("Exp1","Exp2")){
      for(Var in c("NIP","MH63")){
      IntDataC<- subset_samples(IntData, Compartment %in% Com)
      IntDataCE<- subset_samples(IntDataC, Experiment %in% Exp)
      IntDataCEV<- subset_samples(IntDataCE, Variety %in% Var )
	  
      IntDataCEV_Sel = filter_taxa(IntDataCEV,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
      SampleData<-data.frame(sample_data(IntDataCEV_Sel))
      SampleData$Treatment <- factor(SampleData$Treatment,levels=c("Intact","Defoliation"))
      SampleData$Timepoint <- factor(SampleData$Timepoint,unique(SampleData$Timepoint))
      Otutable<-data.frame(otu_table(IntDataCEV_Sel))
	  
	  Nodes<-read.csv("re_int_0.8_0.001_new_diff_nodes2.csv")
	  Moduletable<-c()
	  for(m in unique(Nodes$modularity_class)){
      NodesSel<-subset(Nodes,modularity_class == m)
	  OtutableSel<-Otutable[which(rownames(Otutable) %in% NodesSel$Id),]
	  Moduletable<-rbind(Moduletable,colSums(OtutableSel))
	  }
	  colnames(Moduletable)<-colnames(Otutable)
	  rownames(Moduletable)<-paste("module",unique(Nodes$modularity_class),sep="")
	  ModuleOther<-colSums(Otutable)-colSums(Moduletable)
	  ModuletableAll<-rbind(ModuleOther,Moduletable)
	  rownames(ModuletableAll)[1]<-"module_other"

      Treat<- factor(SampleData$Treatment)
      Time<- factor(substring(SampleData$Timepoint,4,4))
      y <- DGEList(counts=ModuletableAll, group=Treat)
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
        Result$rel.TP1<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day1")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day7")$SampleID)][rownames(Result),])
      }else{
        Result$rel.TP1<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day2")$SampleID)][rownames(Result),])
        Result$rel.TP2<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day3")$SampleID)][rownames(Result),])
        Result$rel.TP3<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day4")$SampleID)][rownames(Result),])
        Result$rel.TP4<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Timepoint == "day5")$SampleID)][rownames(Result),]) 
      }
	  
      Result$rel.Intact<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Treatment == "Intact")$SampleID)][rownames(Result),])
      Result$rel.Defoliation<-rowMeans(ModuletableAll[,which(colnames(ModuletableAll) %in% subset(SampleData,Treatment == "Defoliation")$SampleID)][rownames(Result),])
	  Result$module_name = rownames(Result)
      ResultAll<-rbind(ResultAll,Result)
    }
  }
}
write.csv(ResultAll,"OutputDir/edgeR_module_result_defoliation_endo.csv")


			




