par(family="Times") ##SOOOOOS TIMES NEW ROMAN
require(arules)
require(arulesViz)
colfunc <- colorRampPalette(c("red", "lightgreen"))

##CCLP_Drug
cclp.drug<-read.csv("~/Projects/MBA/data/CCLP_Drug_Categorical.tsv",sep='\t')
cclp.drug[,1]<-gsub(" ","",cclp.drug[,1])
cclp.drug[,1]<-gsub("-","",cclp.drug[,1])
cclp.drug[,1]<-gsub("\\.","",cclp.drug[,1])
cclp.drug[,1]<-gsub("/","",cclp.drug[,1])
cclp.drug[,1]<-toupper(cclp.drug[,1])


load("/home/kvougas/Projects/MBA/new_data/RData/RULES_SIGNIFICANT_Sup3_Conf3in1001_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("/home/kvougas/Projects/MBA/new_data/RData/RULES_TWO-WAY_SIGNIFICANT_Sup8in1001_Conf8in1001.RData")

#load("/home/kvougas/Projects/MBA/data/RULES_2WAY_SIGNIFICANT_SUP7_CONF01_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("/home/kvougas/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")



##One-way
##Top Support

rules.sig<-sort(rules.sig,by="support")
r.sens<-subset(rules.sig,subset=rhs%pin%'Sensitive')
r.res<-subset(rules.sig,subset=rhs%pin%'Resistant')
plot(rules.sig[1:10000],col=colfunc(100))
plot(r.sens[1:1000], method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))
plot(r.res[1:1000], method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))

#plot(rules.sig[1:10], method="graph")

##Top Confidence
rules.sig<-sort(rules.sig,by="confidence")
plot(rules.sig[1:10000],col=colfunc(100))

##Top Lift
rules.sig<-sort(rules.sig,by="lift")
plot(rules.sig[1:10000],col=colfunc(100))
##


##Two-way
##Top Support
rules.2way.sig<-sort(rules.2way.sig,by="support")
plot(rules.2way.sig[1:100000],col=colfunc(100))

##Top Confidence
rules.2way.sig<-sort(rules.2way.sig,by="confidence")
plot(rules.2way.sig[1:100000],col=colfunc(100))

##Top Lift
rules.2way.sig<-sort(rules.2way.sig,by="lift")
plot(rules.2way.sig[1:100000],col=colfunc(100))
##


##All
rules.sig.res<- subset(rules.sig, subset = (rhs %pin% "Resistant"))
rules.sig.sens<- subset(rules.sig, subset = (rhs %pin% "Sensi"))
rules.sig.sens<-sort(rules.sig.sens,by="support")
#rules.sig.sens<-sort(rules.sig.sens,by="lift")
#rules.sig.sens<-sort(rules.sig.sens,by="confidence")
plot(rules.sig.sens[1:2000], method="grouped", interactive = F,control=list(k=200))

rules.sig.res<-sort(rules.sig.res,by="support")
#rules.sig.res<-sort(rules.sig.res,by="lift")
#rules.sig.res<-sort(rules.sig.res,by="confidence")
plot(rules.sig.res[1:2000], method="grouped", interactive = F,control=list(k=200))
##


##Tissue
rules.sig.res<- sort(subset(rules.sig, subset = (lhs %pin% "Tissue=" & rhs %pin% "Resistant")),by="confidence")
rules.sig.sens<- sort(subset(rules.sig, subset = (lhs %pin% "Tissue=" & rhs %pin% "Sensi")),by="confidence")
rules.sig.sens<-sort(rules.sig.sens,by="support")
#rules.sig.sens<-sort(rules.sig.sens,by="lift")
#rules.sig.sens<-sort(rules.sig.sens,by="confidence")
plot(rules.sig.sens, method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))

rules.sig.res<-sort(rules.sig.res,by="support")
#rules.sig.res<-sort(rules.sig.res,by="lift")
#rules.sig.res<-sort(rules.sig.res,by="confidence")
plot(rules.sig.res, method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))
##




#BRAF-MEK
r.sub<-subset(rules.sig,subset=
                rhs%pin%'HG-6-64-1_IC50=Sensitive'|
                rhs%pin%'TL-2-105_IC50=Sensitive'|
                rhs%pin%'VX-11e_IC50=Sensitive'|
                rhs%pin%'FR-180204_IC50=Sensitive'|
                rhs%pin%'RDEA119_IC50=Sensitive'|
                rhs%pin%'CI-1040_IC50=Sensitive'|
                rhs%pin%'PLX4720_IC50=Sensitive'|
                rhs%pin%'SL 0101-1_IC50=Sensitive'|
                rhs%pin%'PD-0325901_IC50=Sensitive'|
                rhs%pin%'SB590885_IC50=Sensitive'|
                rhs%pin%'selumetinib_IC50=Sensitive'|
                rhs%pin%'Trametinib_IC50=Sensitive'|
                rhs%pin%'Dabrafenib_IC50=Sensitive'|
                rhs%pin%'selumetinib_IC50=Sensitive'
                )
r.sub<-sort(r.sub,by="support",decreasing = T)
temp<-plot(r.sub[1:1000], method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))

r.sub<-sort(r.sub,by="lift",decreasing = T)
temp<-plot(r.sub[1:1000], method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))

#PI3K
r.sub<-subset(rules.sig,subset=
                rhs%pin%'AZD6482_IC50=Sensitive'|
                rhs%pin%'OSU-03012_IC50=Sensitive'|
                rhs%pin%'AKT inhibitor VIII_IC50=Sensitive'|
                rhs%pin%'BX-912_IC50=Sensitive'|
                rhs%pin%'ZSTK474_IC50=Sensitive'|
                rhs%pin%'AS605240_IC50=Sensitive'|
                rhs%pin%'KIN001-102_IC50=Sensitive'|
                rhs%pin%'CAL-101_IC50=Sensitive'|
                rhs%pin%'GSK2126458_IC50=Sensitive'|
                rhs%pin%'KIN001-244_IC50=Sensitive'|
                rhs%pin%'PI-103_IC50=Sensitive'|
                rhs%pin%'PIK-93_IC50=Sensitive'|
                rhs%pin%'GSK690693_IC50=Sensitive'|
                rhs%pin%'MK-2206_IC50=Sensitive'|
                rhs%pin%'BEZ235_IC50=Sensitive'|
                rhs%pin%'GDC0941_IC50=Sensitive'|
                rhs%pin%'AZD6482_IC50=Sensitive'|
                rhs%pin%'GDC0941_IC50=Sensitive'
)
r.sub<-sort(r.sub,by="support",decreasing = T)
temp<-plot(r.sub[1:1000], method="grouped", interactive = T,control=list(k=50,main="",col=colfunc(100)))


r.sub<-sort(r.sub,by="lift",decreasing = T)
temp<-plot(r.sub[1:1000], method="grouped", interactive = F,control=list(k=50,main="",col=colfunc(100)))



#########################################
#IC50-Heatmaps
########################################
require(arules)
require(gplots)
load("/home/kvougas/Projects/MBA/new_data/RData/RULES_SIGNIFICANT_Sup3_Conf3in1001_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("/home/kvougas/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")
load("/home/kvougas/Projects/MBA/new_data/RData/drugs_ic50_matrix.RData")
drugs.m.s<-scale(drugs.m,scale = T,center = T)
final.drugs<-read.csv("~/Projects/MBA/new_data/RData/final_drugs.txt",sep='\t')


###Dox si-targets###

#MAGI3#
idx<-which(m.final[,'MAGI3_GE']=="over")
temp.m.final.drugs<-drugs.m.s[idx,]
m.idx<-match(labels(temp.m.final.drugs)[[2]],final.drugs$DRUG.NAME)
temp.m.final.drugs<-temp.m.final.drugs[,m.idx[!is.na(m.idx)]]
#m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
#temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
#temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/new_data/paper_prep/plots/MAGI3over_Heatmap.jpg",width=6000,height=4550,pointsize = 16)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=1,cexRow=2,trace="none",key=F,srtCol = 45,srtRow = 75,dendrogram="none")
dev.off()



#BRAF-MEK BRAF=Mut
drugs<-c(
  'HG-6-64-1',
  'TL-2-105',
  'VX-11e',
  'FR-180204',
  'RDEA119',
  'CI-1040',
  'PLX4720',
  'SL 0101-1',
  'PD-0325901',
  'SB590885',
  'selumetinib',
  'Trametinib',
  'Dabrafenib'
  )

drugs<-unique(drugs)
idx<-which(m.final[,'BRAF']=="Mut")
temp.m.final.drugs<-drugs.m.s[idx,drugs]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
#dimnames(temp.m.final.drugs)<-list(paste(labels(temp.m.final.drugs)[[1]],m.final[idx,1],sep="_"),gsub("_IC50","",labels(temp.m.final.drugs)[[2]]))
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC50","",labels(temp.m.final.drugs)[[2]]))
#jpeg("~/Projects/MBA/data/Viz/BRAF_Mut_Heatmap.jpg",width=2400,height=1820,pointsize = 12)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=0.75,cexRow=0.5,trace="none",key=F,srtCol = 90,srtRow = 0,dendrogram="none")
#dev.off()

#Melanoma
idx<-which(m.final[,'Tissue']=="melanoma")
temp.m.final.drugs<-drugs.m.s[idx,drugs]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
#dimnames(temp.m.final.drugs)<-list(paste(labels(temp.m.final.drugs)[[1]],m.final[idx,1],sep="_"),gsub("_IC50","",labels(temp.m.final.drugs)[[2]]))
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC50","",labels(temp.m.final.drugs)[[2]]))
#jpeg("~/Projects/MBA/data/Viz/BRAF_Mut_Heatmap.jpg",width=2400,height=1820,pointsize = 12)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=0.75,cexRow=0.5,trace="none",key=F,srtCol = 90,srtRow = 0,dendrogram="none")


####PI3K Pathway
drugs<-c(
  'AZD6482',
  'OSU-03012',
  'AKT inhibitor VIII',
  'BX-912',
  'ZSTK474',
  'AS605240',
  'KIN001-102',
  'CAL-101',
  'GSK2126458',
  'KIN001-244',
  'PI-103',
  'PIK-93',
  'GSK690693',
  'MK-2206',
  'BEZ235',
  'GDC0941',
  'AZD6482',
  'GDC0941'
)


drugs<-unique(drugs)
  
idx<-scan("~/Projects/MBA/data/PI3K_SENS_RULES4_IDX.txt")
rules.sub<-rules.sig[idx]
rules.sub<- subset(rules.sub, subset = (lhs %pin% "PIK3CA=Mut" & rhs %pin% "Sensitive"))
temp.drugs<-gsub("\\{|\\}|=Sensitive","",as.character(unlist(inspect(rhs(rules.sub)))))
idx<-which(m.final[,'PIK3CA']=="Mut")
temp.m.final.drugs<-m.final.drugs[idx,]
m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/data/Viz/PIK3CA_Mut_Heatmap.jpg",width=2400,height=1820,pointsize = 12)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=3,cexRow=1,trace="none",key=F,srtCol = 10,srtRow = 33,dendrogram="none")
dev.off()

#Tissue=lung_small_cell_carcinoma
load("~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
rules.sub<- subset(rules.sig, subset = (lhs %pin% "Tissue=lung_small_cell_carcinoma" & rhs %pin% "Resistant"))
#write(rules.sub,file="/home/kvougas/Projects/MBA/data/FMO5_Res_1way.tsv",sep='\t')
temp.drugs<-gsub("\\{|\\}|=Resistant","",as.character(unlist(inspect(rhs(rules.sub)))))
idx<-which(m.final[,'Tissue']=="lung_small_cell_carcinoma")
temp.m.final.drugs<-m.final.drugs[idx,]
m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/data/Viz/Tissue_Small_Cell_Lung_Carcinoma.jpg",width=2400,height=1820,pointsize = 15)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=0.9,cexRow=0.9,trace="none",key=F,srtCol = 20,srtRow = 33,dendrogram="none")
#heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=3,cexRow=1,trace="none",key=F,srtCol = 10,srtRow = 33,dendrogram="none")
dev.off()
#




















##Boxplots for prediction performance
require(ggplot2)
dl.rf.comb<-read.csv("~/Projects/MBA/data/DL_RF_Comb_for_ggplot2.txt",sep='\t')
#Sens
p <- ggplot(dl.rf.comb, aes(x=Classifier, y=Sensitivity)) 
 p + 
  geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  theme(panel.background = element_rect(fill='white'))

#ACC
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=ACC)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #NPV
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=NPV)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #Specificity
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=Specificity)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #PPV
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=PPV)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #FPR
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=FPR)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
