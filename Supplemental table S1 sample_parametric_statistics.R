
rm(list = ls())
library(phyloseq)
OriginalOtuTable<-import_biom("zotus_tab_tax_final.biom")
AbundanceTable<-otu_table(OriginalOtuTable)
AnnotationInformation<-tax_table(OriginalOtuTable)
SampleInformation<-import_qiime_sample_data("sample_data.txt")
IntegratedData<-phyloseq(AbundanceTable,AnnotationInformation,SampleInformation)
colnames(tax_table(IntegratedData))<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
saveRDS(object = IntegratedData,file = "IntegratedData.rds")

SampleInformation$MergeExpTime<-paste(SampleInformation$Experiment,SampleInformation$Timepoint,sep="_")
sample_data(IntegratedData)<-SampleInformation
IntegratedData<- subset_samples(IntegratedData, !MergeExpTime %in% c("Exp2_day1"))

TotalSequenceNumberOfSamples<-colSums(data.frame(otu_table(IntegratedData)))
MitochondrialData <- subset_taxa(IntegratedData, Family %in% "Mitochondria")
TotalNumberOfMitochondria<-colSums(data.frame(otu_table(MitochondrialData)))
ArchaeaData <- subset_taxa(IntegratedData, Kingdom %in% "Archaea")
TotalNumberOfArchaea<-colSums(data.frame(otu_table(ArchaeaData)))
UnassignedData<- subset_taxa(IntegratedData, Kingdom %in% "Unassigned")
TotalNumberOfUnassigned<-colSums(data.frame(otu_table(UnassignedData)))

IntegratedDataRemoveMito <- subset_taxa(IntegratedData, !Family %in% "Mitochondria")
IntegratedDataRemoveMitoAU  <- subset_taxa(IntegratedDataRemoveMito, !Kingdom %in% c("Archaea","Unassigned"))
IntegratedBacteriaData <- subset_taxa(IntegratedDataRemoveMitoAU, !Phylum %in% "Cyanobacteria")
saveRDS(object = IntegratedBacteriaData,file = "IntegratedBacteriaData.rds")

TotalSequenceNumberOfBacteria<-colSums(data.frame(otu_table(IntegratedBacteriaData)))
OtuTable<-data.frame(otu_table(IntegratedBacteriaData))
AllDataTable<-data.frame(sample_id=colnames(OtuTable),all_reads=TotalSequenceNumberOfSamples,all_reads_nonbacteria=TotalSequenceNumberOfSamples-TotalSequenceNumberOfBacteria,all_reads_bacteria=TotalSequenceNumberOfBacteria,all_reads_bacteria_per=paste(round(100*(TotalSequenceNumberOfBacteria/TotalSequenceNumberOfSamples), 4), "%", sep=""))
sample_data<-data.frame(sample_data(IntegratedBacteriaData))
sample_data$sample_id<-rownames(sample_data)
AllDataTableMerge<-merge(AllDataTable,sample_data, by.x = "sample_id")
write.csv(AllDataTableMerge,"OutputDir/AllDataTableMerge.csv")


