##LOAD DATA##
cell.lines<-read.csv("~/Projects/MBA/new_data/Cell_Line_Details.tsv",sep='\t')
mutations<-read.csv("~/Projects/MBA/new_data/CosmicCLP_MutantExport.tsv",sep='\t')
#gene.expression<-read.csv("~/Projects/MBA/new_data/CosmicCLP_CompleteGeneExpression.tsv",sep='\t')
gene.expression<-read.csv("~/Projects/MBA/new_data/CosmicCLP_RawGeneExpression.tsv",sep='\t')
CNV<-read.csv("~/Projects/MBA/new_data/CosmicCLP_CompleteCNA.tsv",sep='\t')
drugs<-read.csv("~/Projects/MBA/new_data/drugs_ic50.tsv",sep='\t')
##

#######################
##CREATE SUB-MATRICES##
#######################

#Mutation Matrix#
rm.idx<-which(mutations$Mutation.Description=="Substitution - coding silent")
mutations<-mutations[-rm.idx,]
rm(rm.idx)
mutations$Gene.name<-toupper(as.character(gsub("_.*","",mutations$Gene.name)))
genes.mut<-unique(mutations$Gene.name)
mutation.m<-matrix(NA,nrow=nrow(cell.lines),ncol=length(genes.mut),dimnames=list(cell.lines$COSMIC.identifier,genes.mut))
for (i in 1:nrow(mutations)){
#for (i in 1:20){
  error.check<-NULL
  error.check<-try(mutation.m[as.character(mutations$ID_sample[i]),as.character(mutations$Gene.name[i])],silent=T)
  if (inherits(error.check, "try-error")) next
  print(i)
  if (is.na(mutation.m[as.character(mutations$ID_sample[i]),as.character(mutations$Gene.name[i])])){
    mutation.m[as.character(mutations$ID_sample[i]),as.character(mutations$Gene.name[i])]<-as.character(mutations$Mutation.AA[i])
  } else {
    mutation.m[as.character(mutations$ID_sample[i]),as.character(mutations$Gene.name[i])]<-paste(mutation.m[as.character(mutations$ID_sample[i]),as.character(mutations$Gene.name[i])],as.character(mutations$Mutation.AA[i]),sep=",")
  }
}
#save(mutation.m,file="~/Projects/MBA/new_data/RData/mutation_matrix.RData")


#Gene expression Matrix#
gene.expression$GENE_NAME<-toupper(gene.expression$GENE_NAME)
genes<-unique(as.character(gene.expression$GENE_NAME))
gene.expr.m<-matrix(NA,nrow=nrow(cell.lines),ncol=length(genes),dimnames=list(cell.lines$COSMIC.identifier,genes))
for (i in 1:nrow(gene.expression)){
  print (i)
  error.check<-NULL
  error.check<-try(gene.expr.m[as.character(gene.expression$SAMPLE_ID[i]),as.character(gene.expression$GENE_NAME[i])],silent=T)
  if (inherits(error.check, "try-error")) next
  gene.expr.m[as.character(gene.expression$SAMPLE_ID[i]),as.character(gene.expression$GENE_NAME[i])]<-gene.expression$GENE_EXPRESSION[i]
}
save(gene.expr.m,file="~/Projects/MBA/new_data/RData/raw_gene_expression_matrix.RData")


#CNV Matrix#
CNV$gene_name<-toupper(CNV$gene_name)
genes<-unique(as.character(CNV$gene_name))
cnv.m<-matrix(NA,nrow=nrow(cell.lines),ncol=length(genes),dimnames=list(cell.lines$COSMIC.identifier,genes))
for (i in 1:nrow(CNV)){
  #for (i in 1:20){
  error.check<-NULL
  error.check<-try(cnv.m[as.character(CNV$ID_SAMPLE[i]),as.character(CNV$gene_name[i])],silent=T)
  if (inherits(error.check, "try-error")) next
  print(i)
    cnv.m[as.character(CNV$ID_SAMPLE[i]),as.character(CNV$gene_name[i])]<-as.character(CNV$MUT_TYPE[i])
  }
save(cnv.m,file="~/Projects/MBA/new_data/RData/cnv_matrix.RData")


#Drug Matrix#
drug.names<-unique(as.character(drugs$DRUG_NAME))
drugs.m<-matrix(NA,nrow=nrow(cell.lines),ncol=length(drug.names),dimnames=list(cell.lines$COSMIC.identifier,drug.names))
for (i in 1:nrow(drugs)){
   print(i)
   try(drugs.m[as.character(drugs$COSMIC_ID[i]),as.character(drugs$DRUG_NAME[i])]<-drugs$LN_IC50[i],silent=T)
}
#save(drugs.m,file="~/Projects/MBA/new_data/RData/drugs_ic50_matrix.RData")


########################
##CREATE MASTER MATRIX##
########################
load("~/Projects/MBA/new_data/RData/mutation_matrix.RData")
load("~/Projects/MBA/new_data/RData/gene_expression_matrix.RData")
load("~/Projects/MBA/new_data/RData/cnv_matrix.RData")
load("~/Projects/MBA/new_data/RData/drugs_ic50_matrix.RData")
cell.lines<-read.csv("~/Projects/MBA/new_data/Cell_Line_Details.tsv",sep='\t')

#Change Gene 2 Gene_CNV for CNV
dimnames(cnv.m)<-list(labels(cnv.m)[[1]],paste(labels(cnv.m)[[2]],"CNV",sep="_"))

#Transform gene expression from numerical to categorical. Over>=2 & Under<=-2. Change Gene 2 Gene_GE
#for _GE_NUM version do not run the next 3 lines
gene.expr.m[gene.expr.m>(-2)&gene.expr.m<2]<-NA
gene.expr.m[gene.expr.m>=2]<-'over'
gene.expr.m[gene.expr.m!="over"&!is.na(gene.expr.m)]<-'under'
dimnames(gene.expr.m)<-list(labels(gene.expr.m)[[1]],paste(labels(gene.expr.m)[[2]],"GE",sep="_"))

#Simplify Mutation. All mutations are 'Mut'
mutation.m[!is.na(mutation.m)]<-'Mut'

#Column-wise Z-score drugs matrix and transform to categorical. Resistant>=1 & Sensitive<=-1. Change Drug 2 Drug_IC50
drugs.m<-apply(drugs.m,2,scale,scale=T,center=T)
drugs.m[drugs.m>(-1)&drugs.m<1]<-NA
drugs.m[drugs.m>=1]<-'Resistant'
drugs.m[drugs.m!="Resistant"&!is.na(drugs.m)]<-'Sensitive'
dimnames(drugs.m)<-list(labels(drugs.m)[[1]],paste(labels(drugs.m)[[2]],"IC50",sep="_"))

##Compile master matrix
m.final<-cbind(as.matrix(cell.lines[,9]),mutation.m,gene.expr.m,cnv.m,drugs.m)
dimnames(m.final)<-list(as.character(cell.lines$Sample.Name),c("Tissue",labels(m.final)[[2]][-1]))
#save(m.final,file="~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")
#save(m.final,file="~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL_GE_NUM.RData")


#Remove empty rows
c<-NA
for (i in 1:nrow(m.final)){
  print(i)
  c[i]<-length(which(is.na(m.final[i,1:max(grep("_CNV",labels(m.final)[[2]],ignore.case = F))])))
}
rm.idx<-which(c==(max(grep("_CNV",labels(m.final)[[2]],ignore.case = F))-1))
if (length(rm.idx)>0) m.final<-m.final[-rm.idx,]

#Remove empty columns
c<-NA
for (i in 1:ncol(m.final)){
  print(i)
  c[i]<-length(which(is.na(m.final[,i])))
}
rm.idx<-which(c==nrow(m.final))
m.final<-m.final[,-rm.idx]
#save(m.final,file="~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")
save(m.final,file="~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL_GE_NUM.RData")

###Do permutations (Parallel Computing) on CATEGORICAL (NOT GE_NUM)
require(parallel)
load("~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")
cl<-makeCluster(7)
clusterExport(cl, "m.final")
m.final.perm<-parSapply (cl=cl, 1:ncol(m.final), function (col) m.final[,col]<-sample(m.final[,col]))
dimnames(m.final.perm)<-list(labels(m.final)[[1]],labels(m.final)[[2]])
stopCluster(cl)
save(m.final.perm,file="~/Projects/MBA/new_data/RData/MASTER_MATRIX_PERMUTED.RData")




###Make Train & Test, normal & permuted ###
load("~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")

tissue.table<-sort(table(m.final[,1]),decreasing=T)
m.final.row.names<-labels(m.final)[[1]]
for (i in 1:length(tissue.table)){
  temp.idx<-which(m.final[,1]==names(tissue.table[i])) 
  temp.m.final<-m.final[temp.idx,]
  temp.m.final.row.names<-m.final.row.names[which(m.final[,1]==names(tissue.table[i])) ]
  train.no<-round((tissue.table[i]/3)*2,digits = 0)
  train.random.idx<-sample(c(1:tissue.table[i]),train.no)
  print(paste(tissue.table[i],train.no,paste(train.random.idx,collapse="-"),sep=","))
  if (i==1){
    train<-temp.m.final[train.random.idx,]
    test<-temp.m.final[-train.random.idx,]
    train.row.names<-temp.m.final.row.names[train.random.idx]
    test.row.names<-temp.m.final.row.names[-train.random.idx]
    next
  }
  if (tissue.table[i]==1){
    train<-rbind(train,temp.m.final)
    train.row.names<-c(train.row.names,temp.m.final.row.names)
    next
  } else {
    train.no<-round((tissue.table[i]/3)*2,digits = 0)
    train.random.idx<-sample(c(1:tissue.table[i]),train.no)
    train<-rbind(train,temp.m.final[train.random.idx,])
    test<-rbind(test,temp.m.final[-train.random.idx,])
    train.row.names<-c(train.row.names,temp.m.final.row.names[train.random.idx])
    test.row.names<-c(test.row.names,temp.m.final.row.names[-train.random.idx])
  }
}

dimnames(train)<-list(train.row.names,labels(m.final)[[2]])
dimnames(test)<-list(test.row.names,labels(m.final)[[2]])

save(train,file="~/Projects/MBA/new_data/RData/TRAIN.RData")
save(test,file="~/Projects/MBA/new_data/RData/TEST.RData")

#Create _GE_NUM version of TRAIN & TEST
load("~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL_GE_NUM.RData")
load("~/Projects/MBA/new_data/RData/TRAIN.RData")
load("~/Projects/MBA/new_data/RData/TEST.RData")
train.row.names<-labels(train)[[1]]
test.row.names<-labels(test)[[1]]

train[,c(19386:35830)]<-m.final[train.row.names,c(19386:35830)]
test[,c(19386:35830)]<-m.final[test.row.names,c(19386:35830)]
save(train,file="~/Projects/MBA/new_data/RData/TRAIN_GE_NUM.RData")
save(test,file="~/Projects/MBA/new_data/RData/TEST_GE_NUM.RData")



###Do permutations (Parallel Computing) on TRAIN
require(parallel)
load("~/Projects/MBA/new_data/RData/TRAIN.RData")
cl<-makeCluster(7)
clusterExport(cl, "train")
train.perm<-parSapply (cl=cl, 1:ncol(train), function (col) train[,col]<-sample(train[,col]))
dimnames(train.perm)<-list(labels(train)[[1]],labels(train)[[2]])
stopCluster(cl)
save(train.perm,file="~/Projects/MBA/new_data/RData/TRAIN_PERMUTED.RData")



###Do permutations (Parallel Computing) on TRAIN & TEST GE_NUM
require(parallel)
load("~/Projects/MBA/new_data/RData/TRAIN_GE_NUM.RData")
load("~/Projects/MBA/new_data/RData/TEST_GE_NUM.RData")
cl<-makeCluster(7)
clusterExport(cl, "train")
clusterExport(cl, "test")
train.perm<-parSapply (cl=cl, 1:ncol(train), function (col) train[,col]<-sample(train[,col]))
dimnames(train.perm)<-list(labels(train)[[1]],labels(train)[[2]])
test.perm<-parSapply (cl=cl, 1:ncol(test), function (col) test[,col]<-sample(test[,col]))
dimnames(test.perm)<-list(labels(test)[[1]],labels(test)[[2]])
stopCluster(cl)
save(train.perm,file="~/Projects/MBA/new_data/RData/TRAIN_GE_NUM_PERMUTED.RData")
save(test.perm,file="~/Projects/MBA/new_data/RData/TEST_GE_NUM_PERMUTED.RData")
