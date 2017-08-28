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

con<-dbConnect(SQLite(),"~/Projects/MBA/new_data/ARM_dlrftuned.db")
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, nthreads=23, max_mem_size="55G")


#drug<-"A-443654_IC50=Resistant"
#drugs<-unique(inspect(rhs(rules.sig))$items)
#drugs<-gsub("\\{|\\}","",drugs)
drugs<-scan("~/Projects/MBA/new_data/RData/GOOD_DRUGS_BATCH.txt",what="character",sep='\n')


for (drug in drugs){
  drug2<-unlist(strsplit(drug,"="))[1]
  drug.status<-unlist(strsplit(drug,"="))[2]
  
  temp.rules<-subset(rules.sig,subset=rhs%pin%drug)
  temp.rules<-subset(temp.rules,subset=support>quantile(quality(temp.rules)$support)[2])#25% Quantile support
  predictor.genes.0<-inspect(lhs(temp.rules))
  predictor.genes.0<-predictor.genes.0$items
  predictor.genes.0<-unique(gsub("\\{|\\}|=.*","",predictor.genes.0))
  rm.idx<-which(predictor.genes.0=="Tissue")
  if (length(rm.idx)>0) predictor.genes.0<-predictor.genes.0[-rm.idx]
  
  
  temp.train<-train[,c(drug2,predictor.genes.0,"Tissue")]
  temp.test<-test[,c(drug2,predictor.genes.0,"Tissue")]
  
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
  if (length(ge.idx)>0) right.columns<-right.columns[-ge.idx]
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
    m.idx.test<-match(total.levels,test.levels)
    missing.test<-total.levels[is.na(m.idx.test)]
    if (length(missing.test)>0){
      for (impute.value in missing.test){
        normal.idx<-which(temp.test[,i]=="Normal")
        temp.test[sample(normal.idx,1),i]<-impute.value
      }
    }
    temp.train[,i]<-as.factor(temp.train[,i])
    temp.test[,i]<-as.factor(temp.test[,i])  
  }
  
  
  temp.train.h2o<-as.h2o(temp.train,destination_frame = "temp.train.h2o")
  temp.test.h2o<-as.h2o(temp.test,destination_frame = "temp.test.h2o")
  
  hyprerparameters<-list(
    activation = c("MaxoutWithDropout","RectifierWithDropout","TanhWithDropout"),
    input_dropout_ratio=c(0.5,0.1),
    hidden_dropout_ratios=list(c(0.1,0.1,0.1),c(0.5,0.5,0.5),c(0.1,0.5,0.1),c(0.5,0.1,0.5)),
    #hidden=list(c(50,50),c(100,100),c(200,200),c(300,300),c(30,100,300),c(300,100,30),c(100,30,100),c(300,30,300),c(300,100,300),c(100,10,100)),
    hidden=list(c(100,100,100),c(200,200,200),c(300,300,300),c(100,200,300),c(300,200,100),c(300,200,300),c(200,300,200)),
    balance_classes=TRUE,
    epochs=c(100,200)
  )
  
  rf.hyprerparameters<-list(
    ntrees=c(round(ncol(temp.train)/2,digits=0),
             round(ncol(temp.train)/4,digits=0),
             (3*round(ncol(temp.train)/4,digits=0)),
             round(ncol(temp.train),digits=0),
             (2*round(ncol(temp.train),digits=0)),
             1000,3000,5000
    ),
    balance_classes=TRUE,
    mtries=c(
      round(sqrt(ncol(temp.train)),digits=0),
      round(sqrt(ncol(temp.train))/2,digits=0),
      round(sqrt(ncol(temp.train))*2,digits=0)
    )
  )
  
  search_criteria = list(strategy = "RandomDiscrete", 
                         max_models = 150,#CHANGE TO 1000!!!!!!! 
                         seed=1234567) 
  
  
  
  
  query<-paste("SELECT Drug,Drug_Status FROM deep_learning_prediction_metrics where Drug='",
               drug2,"' AND Drug_Status='",drug.status,"'",sep="")
  result<-dbGetQuery(con,query)
  rm(query)
  
  if (nrow(result)==0) {
  
    
    ###Deep Learning Grid Search
  
    grid1 <- h2o.grid("deeplearning",
                     x = 2:ncol(temp.train),  
                     y = 1,   
                     training_frame = temp.train.h2o,
                     validation_frame = temp.test.h2o,
                     hyper_params = hyprerparameters,
                     stopping_metric = "AUC",
                     search_criteria = search_criteria
                     )
    
    grid <- h2o.getGrid(grid1@grid_id,sort_by="AUC",decreasing=TRUE)
    #grid
    
    #grid@summary_table[1,]
    best_model <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
    #best_model
    p.h2o.2 <- h2o.predict(best_model, temp.test.h2o)
    p2<-as.data.frame(p.h2o.2)
    p.final<-p2[,3]
    
    pred <- prediction(p.final, temp.test[,1])
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
    
    query<-paste("INSERT INTO deep_learning_prediction_metrics (Drug,Drug_Status,AUC,Sensitivity,Specificity,PPV,NPV,ACC,FPR) VALUES ('",
                 drug2,"','",
                 drug.status,"','",
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
    
  p.final<-NA
  }
    #######################
    
    
    
    
    
    ####Random Fotests Grid Search######
    
  query<-paste("SELECT Drug,Drug_Status FROM random_forests_prediction_metrics where Drug='",
               drug2,"' AND Drug_Status='",drug.status,"'",sep="")
  result<-dbGetQuery(con,query)
  rm(query)
  
  
  if (nrow(result)==0) {
    
    grid2 <- h2o.grid("randomForest",
                      x = 2:ncol(temp.train),  
                      y = 1,   
                      training_frame = temp.train.h2o,
                      validation_frame = temp.test.h2o,
                      hyper_params = rf.hyprerparameters,
                      stopping_metric = "AUC",
                      search_criteria = search_criteria
    )
    
    grid <- h2o.getGrid(grid2@grid_id,sort_by="AUC",decreasing=TRUE)
    #grid
    
    #grid@summary_table[1,]
    best_model <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
    #best_model
    p.h2o.2 <- h2o.predict(best_model, temp.test.h2o)
    p2<-as.data.frame(p.h2o.2)
    p.final<-p2[,3]
    
    
    pred <- prediction(p.final, temp.test[,1])
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
    
    query<-paste("INSERT INTO random_forests_prediction_metrics (Drug,Drug_Status,AUC,Sensitivity,Specificity,PPV,NPV,ACC,FPR) VALUES ('",
                 drug2,"','",
                 drug.status,"','",
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
    
    p.final<-NA
  }
  
    ls_temp <- h2o.ls()
    for (n_ls in 1:nrow(ls_temp)) {
      h2o.rm(as.character(ls_temp[n_ls, 1]))
    }
}



dbDisconnect(con)
h2o.shutdown(prompt = F)



 