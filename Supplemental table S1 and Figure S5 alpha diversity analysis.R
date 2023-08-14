
set.seed(5)
IntData = readRDS("IntegratedBacteriaData.rds")

erie<-IntData
min_lib <- min(sample_sums(erie))

nsamp=nsamples(erie)

trials=100

richness<-matrix(nrow=nsamp,ncol=trials)

row.names(richness)<-sample_names(erie)

evenness<-matrix(nrow=nsamp,ncol=trials)

row.names(evenness)<-sample_names(erie)

Shannon<-matrix(nrow=nsamp,ncol=trials)

row.names(Shannon)<-sample_names(erie)

ACE<-matrix(nrow=nsamp,ncol=trials)

row.names(ACE)<-sample_names(erie)

Fisher<-matrix(nrow=nsamp,ncol=trials)

row.names(Fisher)<-sample_names(erie)

Chao1<-matrix(nrow=nsamp,ncol=trials)

row.names(Chao1)<-sample_names(erie)

Simpson<-matrix(nrow=nsamp,ncol=trials)

row.names(Simpson)<-sample_names(erie)








set.seed(4)
for(i in 1:100){
#Subsample
r<-rarefy_even_depth(erie,
sample.size=min_lib,
verbose=FALSE,
replace=TRUE)

rich<-as.numeric(as.matrix(estimate_richness(r,measures='Observed')))

richness[,i]<-rich

even<-as.numeric(as.matrix(estimate_richness(r,measures='InvSimpson')))

evenness[,i]<-even

Shan<-as.numeric(as.matrix(estimate_richness(r,measures='Shannon')))

Shannon[,i]<-Shan

Cha<-as.numeric(as.matrix(estimate_richness(r,measures='Chao1')[,1]))

Chao1[,i]<-Cha

Sim<-as.numeric(as.matrix(estimate_richness(r,measures='Simpson')))

Simpson[,i]<-Sim


Fis<-as.numeric(as.matrix(estimate_richness(r,measures='Fisher')))

Fisher[,i]<-Fis

A<-as.numeric(as.matrix(estimate_richness(r,measures='ACE')[,1]))

ACE[,i]<-A

}
SampleID<-row.names(richness)
mean<-apply(richness,1,mean)
sd<-apply(richness,1,sd)
measure<-rep('Richness',nsamp)
rich_stats<-data.frame(SampleID,mean,sd,measure)


SampleID<-row.names(evenness)
mean<-apply(evenness,1,mean)
sd<-apply(evenness,1,sd)
measure<-rep('Inverse Simpson',nsamp)
even_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Shannon)
mean<-apply(Shannon,1,mean)
sd<-apply(Shannon,1,sd)
measure<-rep('Shannon',nsamp)
Shannon_stats<-data.frame(SampleID,mean,sd,measure)


SampleID<-row.names(ACE)
mean<-apply(ACE,1,mean)
sd<-apply(ACE,1,sd)
measure<-rep('ACE',nsamp)
ACE_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Simpson)
mean<-apply(Simpson,1,mean)
sd<-apply(Simpson,1,sd)
measure<-rep('Simpson',nsamp)
Simpson_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Fisher)
mean<-apply(Fisher,1,mean)
sd<-apply(Fisher,1,sd)
measure<-rep('Fisher',nsamp)
Fisher_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Chao1)
mean<-apply(Chao1,1,mean)
sd<-apply(Chao1,1,sd)
measure<-rep('Chao1',nsamp)
Chao1_stats<-data.frame(SampleID,mean,sd,measure)

alphadiv<-cbind(as.character(rich_stats$SampleID),rich_stats$mean,even_stats$mean,Shannon_stats$mean,ACE_stats$mean,Simpson_stats$mean,Fisher_stats$mean,Chao1_stats$mean)

colnames(alphadiv)<-c("SampleID","richness","evenness","Shannon","ACE","Simpson","Fisher","Chao1")
write.csv(alphadiv,"alphadiv.csv")
