###############################################################
################START FROM HERE
###############################################################

par(family="Times") ##SOOOOOS TIMES NEW ROMAN
load("~/Projects/MBA/new_data/RData/TRAIN_GE_NUM.RData")
load("~/Projects/MBA/new_data/RData/TEST_GE_NUM.RData")

train.ge.idx<-grep("_GE",labels(train)[[2]])
test.ge.idx<-grep("_GE",labels(test)[[2]])

train.drug.idx<-grep("IC50",labels(train)[[2]])
test.drug.idx<-grep("IC50",labels(test)[[2]])

temp.train<-train[,train.ge.idx]
temp.test<-test[,test.ge.idx]
temp.train[is.na(temp.train)]<-"0"
temp.test[is.na(temp.test)]<-"0"
train[,train.ge.idx]<-temp.train
test[,test.ge.idx]<-temp.test
temp.test<-test[,2:(min(test.drug.idx)-1)] 
temp.train<-train[,2:(min(train.drug.idx)-1)] 
temp.test[is.na(temp.test)]<-"Normal"
temp.train[is.na(temp.train)]<-"Normal"
test[,2:(min(test.drug.idx)-1)]<-temp.test 
train[,2:(min(train.drug.idx)-1)]<-temp.train 
train<-as.data.frame(train)
test<-as.data.frame(test)
train.ge.idx<-grep("_GE",names(train))
test.ge.idx<-grep("_GE",names(test))
for (i in train.ge.idx){
  train[,i]<-as.numeric(levels(train[,i]))[train[,i]]
}
for (i in test.ge.idx){
  test[,i]<-as.numeric(levels(test[,i]))[test[,i]]
}

#Remove rows which do not contain gene expression info
train.rm.idx<-which(train[,28000]==0)
test.rm.idx<-which(test[,28000]==0)
train<-train[-train.rm.idx,]
test<-test[-test.rm.idx,]

##################
##################

##########################
###Start actual prediction
##########################
load("~/Projects/MBA/new_data/RData/TRAIN_RULES_SIGNIFICANT_Sup3_Conf3in1001_DYNAMIC_THRESHOLD_FDR0.05.RData")
require(arules)
require(ROCR)
require(stringr)
require(h2o)
require(RSQLite)

con<-dbConnect(SQLite(),"~/Projects/MBA/new_data/ARM.db")

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, nthreads=6, max_mem_size="12G")


drugs<-unique(inspect(rhs(rules.sig))$items)
drugs<-gsub("\\{|\\}","",drugs)

result<-dbGetQuery(con,"SELECT Drug,Drug_Status FROM random_forests_prediction_metrics ORDER BY id DESC LIMIT 1")

if (nrow(result)>0) {
temp.drug<-paste(result$Drug,result$Drug_Status,sep="=")
idx<-which(drugs==temp.drug)
drugs<-drugs[-c(1:idx)]
rm(temp.drug)
rm(idx)
}


for (drug in drugs){
  
  temp.rules<-subset(rules.sig,subset=rhs%pin%drug)
  temp.rules<-subset(temp.rules,subset=support>quantile(quality(temp.rules)$support)[2])#25% Quantile support
  predictor.genes.0<-inspect(lhs(temp.rules))
  predictor.genes.0<-predictor.genes.0$items
  predictor.genes.0<-unique(gsub("\\{|\\}|=.*","",predictor.genes.0))
  rm.idx<-which(predictor.genes.0=="Tissue")
  drug2<-unlist(strsplit(drug,"="))[1]
  drug.status<-unlist(strsplit(drug,"="))[2]
  if (length(rm.idx)>0) predictor.genes.0<-predictor.genes.0[-rm.idx]
  
  if (length(grep("=",predictor.genes.0))>0) predictor.genes.0<-predictor.genes.0[-grep("=",predictor.genes.0)]
  if (length(predictor.genes.0)==0) {
    query<-paste("INSERT INTO random_forests_prediction_metrics (Drug,Drug_Status,Predictor_Genes_No) VALUES ('",
                  drug2,"','",drug.status,"','0')",sep="")
    dbSendQuery(con,query)
    rm(query)
    next
    }
  no.of.rounds<-round(length(predictor.genes.0)/200,digits=0)*3
  if (no.of.rounds<1) no.of.rounds<-1  
  print(paste(drug,no.of.rounds,length(predictor.genes.0),sep=","))
  temp.auc<-NA
  temp.auc.counter<-1
  for (bagging in 1:no.of.rounds){
    if (length(predictor.genes.0)>=200) {
      predictor.genes<-sample(predictor.genes.0,200)
    } else {
      predictor.genes<-predictor.genes.0
    }
    temp.train<-train[,c(drug2,predictor.genes,"Tissue")]
    temp.test<-test[,c(drug2,predictor.genes,"Tissue")]
    
    if(drug.status=="Sensitive"){
      temp.train[,1]<-as.character(temp.train[,1])
      temp.test[,1]<-as.character(temp.test[,1])
      temp.train[is.na(temp.train[,1]),1]<-"Non-Sensitive"
      temp.train[temp.train[,1]=="Resistant",1]<-"Non-Sensitive"
      temp.test[is.na(temp.test[,1]),1]<-"Non-Sensitive"
      temp.test[temp.test[,1]=="Resistant",1]<-"Non-Sensitive"
      temp.train[,1]<-as.factor(temp.train[,1])
      temp.test[,1]<-as.factor(temp.test[,1])
    } else {
      temp.train[,1]<-as.character(temp.train[,1])
      temp.test[,1]<-as.character(temp.test[,1])
      temp.train[is.na(temp.train[,1]),1]<-"Non-Resistant"
      temp.train[temp.train[,1]=="Sensitive",1]<-"Non-Resistant"
      temp.test[is.na(temp.test[,1]),1]<-"Non-Resistant"
      temp.test[temp.test[,1]=="Sensitive",1]<-"Non-Resistant"
      temp.train[,1]<-as.factor(temp.train[,1])
      temp.test[,1]<-as.factor(temp.test[,1])
    }
    
    ge.idx<-grep("_GE",labels(temp.train)[[2]])
    right.columns<-c(1:ncol(temp.train))
    right.columns<-right.columns[-ge.idx]
    right.columns<-right.columns[-1]
    right.columns<-right.columns[-length(right.columns)]
    
    for(i in right.columns){
      temp.train[,i]<-as.character(temp.train[,i])
      temp.test[,i]<-as.character(temp.test[,i])
      train.levels<-sort(unique(temp.train[,i]))
      test.levels<-sort(unique(temp.test[,i]))
      total.levels<-sort(unique(c(train.levels,test.levels)))
      m.idx.train<-match(total.levels,train.levels)
      missing.train<-total.levels[is.na(m.idx.train)]
      if (length(missing.train)>0){
        for (impute.value in missing.train){
          normal.idx<-which(temp.train[,i]=="Normal")
          temp.train[sample(normal.idx,1),i]<-impute.value
          #temp.train[sample(c(1:nrow(temp.train)),1),i]<-impute.value
        }
      }
      temp.train[,i]<-as.factor(temp.train[,i])
      temp.test[,i]<-as.factor(temp.test[,i])  
    }
    
    
    temp.train.h2o<-as.h2o(temp.train,destination_frame = "temp.train.h2o")
    temp.test.h2o<-as.h2o(temp.test,destination_frame = "temp.test.h2o")
    
    model2 <-
      h2o.randomForest(x = 2:ncol(temp.train),  # column numbers for predictors
                       y = 1,   # column number for label
                       training_frame = temp.train.h2o,
                       ntrees=round(length(predictor.genes)/2,digits=0),
                       balance_classes=T,
                       nfolds = 5 
      ) 
    temp.auc[temp.auc.counter]<-model2@model$cross_validation_metrics@metrics$AUC
    temp.auc.counter=temp.auc.counter+1
    
    p.h2o.2 <- h2o.predict(model2, temp.test.h2o)
    p2<-as.data.frame(p.h2o.2)
    p2<-p2[,-1]
    p<-p2
    if (bagging==1) {p.comb<-p} else {p.comb<-cbind(p.comb,p)}
  }
  p.final<-p.comb[,seq(2,ncol(p.comb),2)]
  p.final<-apply(p.final,1,weighted.mean,temp.auc)
  
  
  
  if (sum(p.final,na.rm=T)==0) next
  
  pred <- prediction(p.final, temp.test[,1])
  perf <- performance(pred,"tpr","fpr")
  
  print("CHECK")
  drug3<-drug
  drug3<-gsub("\\/","_",drug3)
  filename<-paste("~/Projects/MBA/new_data/Prediction_h2o_rf_quantile_support_filter_GE_NUM/plots/",drug3,".png",sep="")
  png(filename,width=640,height=640,pointsize=12)
  plot(perf,colorize=T,main=gsub("\\.","-",gsub("_.*="," => ",drug)))
  lines(x=c(0,1),y=c(0,1),lty=2)
  text(0.03,1,paste("AUC",try(round(unlist(performance(pred,"auc")@y.values),digits = 2),silent=T),sep="="))
  dev.off()
  
  temp.auc2<-NULL
  temp.sens<-NULL
  temp.spec<-NULL
  temp.ppv<-NULL
  temp.npv<-NULL
  temp.acc<-NULL
  temp.fpr<-NULL
  
  try(temp.auc2<-unlist(performance(pred,"auc")@y.values),silent=T)
  try(temp.sens<-unlist(performance(pred,"sens")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(temp.spec<-unlist(performance(pred,"spec")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(temp.ppv<-unlist(performance(pred,"ppv")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(temp.npv<-unlist(performance(pred,"npv")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(temp.acc<-unlist(performance(pred,"acc")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(temp.fpr<-unlist(performance(pred,"fpr")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)

  query<-paste("INSERT INTO random_forests_prediction_metrics (Drug,Drug_Status,Predictor_Genes,Predictor_Genes_No,AUC,Sensitivity,Specificity,PPV,NPV,ACC,FPR) VALUES ('",
               drug2,"','",
               drug.status,"','",
               paste(predictor.genes.0,collapse=","),"','",
               length(predictor.genes.0),"','",
               temp.auc2,"','",
               temp.sens,"','",
               temp.spec,"','",
               temp.ppv,"','",
               temp.npv,"','",
               temp.acc,"','",
               temp.fpr,"')",
               sep=""
  )
  
  dbSendQuery(con,query)
  rm(query)
  
    predictor.genes<-NA
  p.comb<-NA
  p.final<-NA
  
  ls_temp <- h2o.ls()
  for (n_ls in 1:nrow(ls_temp)) {
    h2o.rm(as.character(ls_temp[n_ls, 1]))
    }
}
dbDisconnect(con)
h2o.shutdown(prompt = TRUE)
