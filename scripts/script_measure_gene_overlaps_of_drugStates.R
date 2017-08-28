rules4<-read.csv("~/Projects/MBA/new_data/paper_prep/data/Drugs_predictorGenes.tsv",sep='\t',colClasses = "character")


drugStatus<-paste(rules4$Drug,rules4$Drug_Status,sep="=")
m<-matrix(NA,nrow=length(drugStatus),ncol=length(drugStatus),dimnames = list(drugStatus,drugStatus))
for (i in 1:length(drugStatus)){
  m[drugStatus[i],drugStatus[i]]<-as.numeric(rules4$Predictor_Genes_No[i])
}

comb<-combn(drugStatus,2)
for (i in 1:ncol(comb)){
  print(i)
  drug1<-unlist(strsplit(comb[1,i],"="))
  drug2<-unlist(strsplit(comb[2,i],"="))
  drug1.idx<-which(rules4$Drug==drug1[1]&rules4$Drug_Status==drug1[2])
  drug2.idx<-which(rules4$Drug==drug2[1]&rules4$Drug_Status==drug2[2])
  pred.genes.drug1<-unlist(strsplit(rules4$Predictor_Genes[drug1.idx],","))
  pred.genes.drug1<-unique(gsub("_GE","",pred.genes.drug1))
  pred.genes.drug1<-unique(gsub("_CNV","",pred.genes.drug1))
  pred.genes.drug2<-unlist(strsplit(rules4$Predictor_Genes[drug2.idx],","))
  pred.genes.drug2<-unique(gsub("_GE","",pred.genes.drug2))
  pred.genes.drug2<-unique(gsub("_CNV","",pred.genes.drug2))
  m.idx<-match(pred.genes.drug1,pred.genes.drug2)
  overlap<-length(m.idx[!is.na(m.idx)])
  m[comb[2,i],comb[1,i]]<-overlap
  rm(m.idx)
}
write.table(m,file="~/Projects/MBA/new_data/paper_prep/data/Predictor_Genes_Overlap_rules4.csv",sep='\t')
save(m,file="~/Projects/MBA/new_data/paper_prep/data/Predictor_Genes_Overlap_rules4.RData")

#####################################
####Calculate overlap probability####
#####################################
load("~/Projects/MBA/new_data/paper_prep/data/Predictor_Genes_Overlap_rules4.RData")
rules4<-read.csv("~/Projects/MBA/new_data/paper_prep/data/Drugs_predictorGenes.tsv",sep='\t',colClasses = "character")
drugStatus<-paste(rules4$Drug,rules4$Drug_Status,sep="=")
comb<-combn(drugStatus,2)
m.p<-matrix(NA,nrow=length(drugStatus),ncol=length(drugStatus),dimnames = list(drugStatus,drugStatus))
overlap<-NA

#Take all genes in the rules
for (i in 1:nrow(rules4)){
  temp.genes<-unlist(strsplit(rules4$Predictor_Genes[i],","))
  temp.genes<-unique(gsub("_GE","",temp.genes))
  temp.genes<-unique(gsub("_CNV","",temp.genes))
  if (i==1) {all.genes<-temp.genes} else {all.genes<-c(all.genes,temp.genes)}
}
all.genes.table<-table(all.genes)
all.genes.unique<-names(all.genes.table)

for (i in 1:ncol(comb)){
print(i)
drug1.gene.no<-m[comb[1,i],comb[1,i]]
drug2.gene.no<-m[comb[2,i],comb[2,i]]
for (f in 1:50){
rand.genes.drug1<-sample(all.genes.unique,drug1.gene.no,prob=all.genes.table)
rand.genes.drug2<-sample(all.genes.unique,drug2.gene.no,prob=all.genes.table)
m.idx<-match(rand.genes.drug1,rand.genes.drug2)
overlap[f]<-length(m.idx[!is.na(m.idx)])  
}
if (m[comb[2,i],comb[1,i]]>=mean(overlap)) {
  m.p[comb[2,i],comb[1,i]]<-(1-pnorm(m[comb[2,i],comb[1,i]],mean=mean(overlap),sd=sd(overlap)))*2  
  m.p[comb[1,i],comb[2,i]]<-"over"
} else {
  m.p[comb[2,i],comb[1,i]]<-(pnorm(m[comb[2,i],comb[1,i]],mean=mean(overlap),sd=sd(overlap)))*2  
  m.p[comb[1,i],comb[2,i]]<-"under"
  }

}

write.table(m.p,file="~/Projects/MBA/new_data/paper_prep/data/Predictor_Genes_Overlap_rules4_p_values.csv",sep='\t')
save(m.p,file="~/Projects/MBA/new_data/paper_prep/data/Predictor_Genes_Overlap_rules4_p_values.RData")

