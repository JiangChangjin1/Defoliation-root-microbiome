rm(list = ls())
library(phyloseq)
library(ggsignif)
library(tidyr)
library(Rmisc)
library(ggplot2)
IntData = readRDS("IntegratedBacteriaData.rds")
IntDataRarefy<-rarefy_even_depth(IntData,rngseed=1,sample.size= min(colSums(otu_table(IntData))), replace=T)
IntDataRarefyPlant<- subset_samples(IntDataRarefy, !Compartment %in% "BS")



SampleData<-sample_data(IntDataRarefyPlant)
AllBetaTable<-c()
for(Exp in c("Exp1","Exp2")){
    for(Com in c("Root","Rhizosphere")){
        for(Var in c("NIP","MH63")){
            if(Exp == "Exp1"){
               time_all<-c("day1","day3","day5","day7")
               }
            if(Exp == "Exp2"){
               time_all<-c("day2","day3","day4","day5")  
               }
            for(Tim in time_all){
                IntDataRarefyPlantE<- subset_samples(IntDataRarefyPlant, Experiment %in% Exp)
                IntDataRarefyPlantEC<- subset_samples(IntDataRarefyPlantE, Compartment %in% Com)
                IntDataRarefyPlantECV<- subset_samples(IntDataRarefyPlantEC, Variety %in% Var)
                IntDataRarefyPlantECVT<- subset_samples(IntDataRarefyPlantECV, Timepoint %in% Tim)
					
                Otutable<-t(otu_table(IntDataRarefyPlantECVT))
				BrayCurtis = vegan::vegdist(Otutable, method = "bray") 
				BrayCurtis= as.matrix(BrayCurtis)
				IntDataRarefyPlantECVT_Int<- subset_samples(IntDataRarefyPlantECVT, Treatment %in% "Intact")
				IntName<-rownames(sample_data(IntDataRarefyPlantECVT_Int))
				IntDataRarefyPlantECVT_Def<- subset_samples(IntDataRarefyPlantECVT, Treatment %in% "Defoliation")
					
				DefName<-rownames(sample_data(IntDataRarefyPlantECVT_Def))
				IntNameUniq<-unique(as.vector(BrayCurtis[IntName,IntName]))
				IntNameUniq<-IntNameUniq[-which(IntNameUniq == 0)]
				DefNameUniq<-unique(as.vector(BrayCurtis[IntName,DefName]))
					
				BetaTable<-data.frame(distance=c(IntNameUniq,DefNameUniq),group=paste(Exp,Com,Var,Tim,sep="_"),BetaGroup=rep(c("Within treatments","Between treatments"),c(length(IntNameUniq) ,length(DefNameUniq) )))
				AllBetaTable<-rbind(AllBetaTable,BetaTable)
				}
            }
		}
	}

AllBetaTable<- separate(AllBetaTable, group, into = c("Experiment", "Compartment", "Variety","Timepoint"), sep = "_") 
AllBetaTable$MergeTimeBeta<-paste(AllBetaTable$Timepoint,AllBetaTable$BetaGroup,sep=" ")
p<-list()
i=1
for(Exp in c("Exp1","Exp2")){ 
	for(Var in c("NIP","MH63")){ 
	    for(Com in c("Root","Rhizosphere")){ 
	        AllBetaTableE<- subset(AllBetaTable, Experiment == Exp)
	        AllBetaTableEV<- subset(AllBetaTableE, Variety == Var)
            AllBetaTableEVC<- subset(AllBetaTableEV, Compartment == Com)
            if(Exp == "Exp1"){
               AllBetaTableEVC$MergeTimeBeta<-factor(AllBetaTableEVC$MergeTimeBeta,levels = c("day1 Within treatments","day1 Between treatments","day3 Within treatments","day3 Between treatments","day5 Within treatments","day5 Between treatments","day7 Within treatments","day7 Between treatments"))
               }
            if(Exp == "Exp2"){
               AllBetaTableEVC$MergeTimeBeta<-factor(AllBetaTableEVC$MergeTimeBeta,levels =  c("day2 Within treatments","day2 Between treatments","day3 Within treatments","day3 Between treatments","day4 Within treatments","day4 Between treatments","day5 Within treatments","day5 Between treatments"))
               }
            if(Exp == "Exp1"){
               errorbar_up<-function(x){ 
               mean(x)+sd(x)
               }
               errorbar_down<-function(x){
               mean(x)-sd(x)
               }
               p[[i]]<-ggplot(AllBetaTableEVC,aes(x = MergeTimeBeta, y = distance, color = BetaGroup,fill=BetaGroup)) +
                       geom_bar(stat = 'summary',fun=mean)+
					   stat_summary(geom="errorbar",               
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.6,color = "black")+
                       geom_signif(
                       textsize  = 4,
                       comparisons = list(
                          c("day1 Within treatments","day1 Between treatments"),
                          c("day3 Within treatments","day3 Between treatments"),
                          c("day5 Within treatments","day5 Between treatments"),
                          c("day7 Within treatments","day7 Between treatments")), 
                       map_signif_level = T, 
                       test = wilcox.test, 
                       vjust=0.2, 
                       tip_length = 0.05,col="black", y_position =rep(c(0.7,0.7,0.75),6) 
                       )+
                       theme_bw() + 
                       xlab("") +
                       ylab("Bray-Curtis dissimilarity")+
                       theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
                       coord_cartesian(clip = "off")+ 
                       scale_color_manual(values =c("#E37933","#406CB4"))+
                       scale_fill_manual(values =c("#E37933","#406CB4"))+
                       theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                       theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")+
                       geom_jitter(alpha = 0.6,color = "black")
               }
            if(Exp == "Exp2"){
               errorbar_up<-function(x){ 
               mean(x)+sd(x)
               }
               errorbar_down<-function(x){
               mean(x)-sd(x)
               }
               p[[i]]<-ggplot(AllBetaTableEVC,aes(x = MergeTimeBeta, y = distance, color = BetaGroup,fill=BetaGroup)) +
                      geom_bar(stat = 'summary',fun=mean)+
                      stat_summary(geom="errorbar",               
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.6,color = "black")+
                        geom_signif(
                        textsize  = 4,
                        comparisons = list(c("day2 Within treatments","day2 Between treatments"),
                                           c("day3 Within treatments","day3 Between treatments"),
                                           c("day4 Within treatments","day4 Between treatments"),
                                           c("day5 Within treatments","day5 Between treatments")), 
                        map_signif_level = T, 
                        test = wilcox.test, 
                        vjust=0.2, 
                        tip_length = 0.05,col="black", y_position =rep(c(0.7,0.7,0.75),6) 
                      )+
                      theme_bw() + 
                      xlab("") +
                      ylab("Bray-Curtis dissimilarity")+
                      theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
                      coord_cartesian(clip = "off")+ 
                      scale_color_manual(values =c("#E37933","#406CB4"))+
                      scale_fill_manual(values =c("#E37933","#406CB4"))+
                      theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())+
                      geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")+
                      geom_jitter(alpha = 0.6,color = "black")
                }
            i=i+1
			}
	    }
	}        

pdf("Beta_Bar_plot8-replaceT-wt2-f.pdf",width=15,height=8)
multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],cols =4)
dev.off()
            
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			beta_data_all_tre<-AllBetaTable
            beta_data_all_tre$beta_group<-beta_data_all_tre$BetaGroup
            beta_data_all_tre0<- separate(beta_data_all_tre, group, into = c("Experiment", "Compartment", "Variety","Timepoint"), sep = "_") 
            
            
            beta_data_all_tre0$group_comb<-paste(beta_data_all_tre0$Timepoint,beta_data_all_tre0$beta_group,sep=" ")
            
            
            
            p<-list()
            i=1
            for(Exp in c("Exp1","Exp2")){ 
              for(Var in c("NIP","MH63")){ 
                for(Com in c("Root","Rhizosphere")){ 
                  #Exp = "Exp1"
                  #Var = "NIP"
                  #Com = "Root"
                  beta_data_all_tre1<- subset(beta_data_all_tre0, Experiment == Exp)
                  beta_data_all_tre2<- subset(beta_data_all_tre1, Variety == Var)
                  beta_data_all_tre3<- subset(beta_data_all_tre2, Compartment == Com)
                  #beta_data_all_tre3$group_comb<-factor(beta_data_all_tre3$group_comb,levels = c("day1 Within treatments","day1 Between treatments","day3 Within treatments","day3 Between treatments","day5 Within treatments","day5 Between treatments","day7 Within treatments","day7 Between treatments"))
                  
                  #beta_data_all_tre3$group_m<-paste(beta_data_all_tre3$Timepoint,beta_data_all_tre3$group_comb,sep="_")
                  beta_data_all_tre3$group_m<-beta_data_all_tre3$group_comb
                  if(Exp == "Exp1"){
                    beta_data_all_tre3$group_m<-factor(beta_data_all_tre3$group_m,levels = c("day1 Within treatments","day1 Between treatments","day3 Within treatments","day3 Between treatments","day5 Within treatments","day5 Between treatments","day7 Within treatments","day7 Between treatments")
                    )
                  }
                  
                  if(Exp == "Exp2"){
                    beta_data_all_tre3$group_m<-factor(beta_data_all_tre3$group_m,levels =  c("day2 Within treatments","day2 Between treatments","day3 Within treatments","day3 Between treatments","day4 Within treatments","day4 Between treatments","day5 Within treatments","day5 Between treatments")
                    )
                  }
                  
                  
                  
                  
                  if(Exp == "Exp1"){
                    errorbar_up<-function(x){ #误差线上下限公式 mean-sd到mena+sd
                      mean(x)+sd(x)
                    }
                    errorbar_down<-function(x){
                      mean(x)-sd(x)
                    }
                    p[[i]]<-ggplot(beta_data_all_tre3,aes(x = group_m, y = distance, color = beta_group,fill=beta_group)) +
                      #geom_rect(xmin = 0.4, xmax = 2.5,
                      #          ymin = -Inf, ymax = Inf,
                      #          fill ='#d2dbdf',
                      #          inherit.aes = F)+
                      #geom_point(alpha = 0.6)+
                      
                      geom_bar(stat = 'summary',fun=mean)+
                      stat_summary(geom="errorbar",               #绘制误差线
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.6,color = "black")+
                      
                      #geom_boxplot(alpha = 0.5) +
                      #scale_color_manual(values = color)+
                      # 先算一下显著性差异，再手动添加
                      # geom_signif(comparisons = list(c("Within treatments intact","Between treatments intact"),c("Within treatments defoliation","Between treatments defoliation")),
                      #            test = "t.test",
                      #            map_signif_level = F)+
                      #annotate("text", x = 0.5, y = 0.6, label ="***",size = 4)+
                      geom_signif(
                        textsize  = 4,
                        comparisons = list(#c("day2day3_Within treatments intact","day2day3_Between treatments intact"),c("day2day3_Within treatments defoliation","day2day3_Between treatments defoliation"),
                          # c("day2day4_Within treatments intact","day2day4_Between treatments intact"),c("day2day4_Within treatments defoliation","day2day4_Between treatments defoliation"),
                          #c("day2day5_Within treatments intact","day2day5_Between treatments intact"),c("day2day5_Within treatments defoliation","day2day5_Between treatments defoliation"),
                          #c("day3day4_Within treatments intact","day3day4_Between treatments intact"),c("day3day4_Within treatments defoliation","day3day4_Between treatments defoliation"),
                          #c("day3day5_Within treatments intact","day3day5_Between treatments intact"),c("day3day5_Within treatments defoliation","day3day5_Between treatments defoliation"),
                          #c("day4day5_Within treatments intact","day4day5_Between treatments intact"),c("day4day5_Within treatments defoliation","day4day5_Between treatments defoliation")
                          c("day1 Within treatments","day1 Between treatments"),
                          c("day3 Within treatments","day3 Between treatments"),
                          c("day5 Within treatments","day5 Between treatments"),
                          c("day7 Within treatments","day7 Between treatments")), #检测两者之间的差异显著性
                        map_signif_level = T, #添加星号标记
                        test = wilcox.test, #检测方法
                        vjust=0.2, #标注和横线的距离
                        tip_length = 0.05,col="black", y_position =rep(c(0.7,0.7,0.75),6) #两端短竖线的长度
                      )+
                      theme_bw() + 
                      xlab("") +
                      ylab("Bray-Curtis dissimilarity")+
                      theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+#+#+
                      #ylim(0.1,0.8)+
                      coord_cartesian(clip = "off")+ ###打破线条只能在框内添加的限制
                      scale_color_manual(values =c("#E37933","#406CB4"))+
                      scale_fill_manual(values =c("#E37933","#406CB4"))+
                      #annotate(geom="text",x=1.5,y=0.9,label="day1")+#,vjust=-16)+
                      #geom_segment(x=4,xend=4,y=0.85,yend=0.9)+
                      theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                      
                      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())+
                      geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")+
                      geom_jitter(alpha = 0.6,color = "black")
                    
                    
                    
                  }
                  
                  
                  
                  if(Exp == "Exp2"){
                    errorbar_up<-function(x){ #误差线上下限公式 mean-sd到mena+sd
                      mean(x)+sd(x)
                    }
                    errorbar_down<-function(x){
                      mean(x)-sd(x)
                    }
                    p[[i]]<-ggplot(beta_data_all_tre3,aes(x = group_m, y = distance, color = beta_group,fill=beta_group)) +
                      #geom_rect(xmin = 0.4, xmax = 2.5,
                      #          ymin = -Inf, ymax = Inf,
                      #          fill ='#d2dbdf',
                      #          inherit.aes = F)+
                      #geom_point(alpha = 0.6)+
                     
                      geom_bar(stat = 'summary',fun=mean)+
                      stat_summary(geom="errorbar",               #绘制误差线
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.6,color = "black")+
                    
                      #geom_boxplot(alpha = 0.5) +
                      #scale_color_manual(values = color)+
                      # 先算一下显著性差异，再手动添加
                      # geom_signif(comparisons = list(c("Within treatments intact","Between treatments intact"),c("Within treatments defoliation","Between treatments defoliation")),
                      #            test = "t.test",
                      #            map_signif_level = F)+
                      #annotate("text", x = 0.5, y = 0.6, label ="***",size = 4)+
                      geom_signif(
                        textsize  = 4,
                        comparisons = list(c("day2 Within treatments","day2 Between treatments"),
                                           c("day3 Within treatments","day3 Between treatments"),
                                           c("day4 Within treatments","day4 Between treatments"),
                                           c("day5 Within treatments","day5 Between treatments")), #检测两者之间的差异显著性), #检测两者之间的差异显著性
                        map_signif_level = T, #添加星号标记
                        test = wilcox.test, #检测方法
                        vjust=0.2, #标注和横线的距离
                        tip_length = 0.05,col="black", y_position =rep(c(0.7,0.7,0.75),6) #两端短竖线的长度
                      )+
                      theme_bw() + 
                      xlab("") +
                      ylab("Bray-Curtis dissimilarity")+
                      theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+#+#+
                      #ylim(0.1,0.8)+
                      coord_cartesian(clip = "off")+ ###打破线条只能在框内添加的限制
                      scale_color_manual(values =c("#E37933","#406CB4"))+
                      scale_fill_manual(values =c("#E37933","#406CB4"))+
                      #annotate(geom="text",x=1.5,y=0.9,label="day1")+#,vjust=-16)+
                      #geom_segment(x=4,xend=4,y=0.85,yend=0.9)+
                      theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                      
                      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())+
                      geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")+
                      geom_jitter(alpha = 0.6,color = "black")
                    
                    
                    
                  }
                  
                  
                  
                  i=i+1
                  
                }
              }
            }
      
            
            library(Rmisc)
            pdf("test_01_15x8_beta_intact_and_def_leave_within_intact_color3.pdf",width=15,height=8)
            multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],cols =4)
            
            dev.off()
            
			