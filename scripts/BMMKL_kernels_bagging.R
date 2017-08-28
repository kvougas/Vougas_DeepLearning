# load library for gaussian kernel calculation (For real-valued genomic views)
library("KRLS")
library("proxy")

# package for rules extraction
require(arules)
require(ROCR)

#BMMKL classification source: https://github.com/mehmetgonen/bemkl
source("bemkl_supervised_classification_variational_train.R")
source("bemkl_supervised_classification_variational_test.R")

# load rules, training and test sets
load("TRAIN_RULES_SIGNIFICANT_Sup3_Conf3in1001_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("TEST_GE_NUM.RData")
load("TRAIN_GE_NUM.RData")

train.set <- train
test.set <- test


drugs<-unique(inspect(rhs(rules.sig))$items)
drugs<-gsub("\\{|\\}","",drugs)

# drugs_res <- drugs[grepl("=Resistant",drugs)]
# drug_res<-as.character(gsub("\\=.*","",drugs_res))
# drug_res_files<- gsub("/", "_", drug_res)
# 
drugs_sens <- drugs[grepl("=Sensitive",drugs)]
drug_sens<-as.character(gsub("\\=.*","",drugs_sens))
drug_sens_files<- gsub("/", "_", drug_sens)

drugs<-gsub("\\=.*","",drugs)
drugs<- unique(drugs)
predictor.genes.sens <- as.character()
genes.count <- as.numeric()
auc_final_bagging <- data.frame()

for (i in 229:length(drug_sens_files)){
  temp.rules<-subset(rules.sig,subset=rhs%pin%drugs_sens[i])
  temp.rules<-subset(temp.rules,subset=support>quantile(quality(temp.rules)$support)[2])#25% Quantile support
  predictor.genes.0<-inspect(lhs(temp.rules))
  predictor.genes.0<-predictor.genes.0$items
  predictor.genes.0<-unique(gsub("\\{|\\}|=.*","",predictor.genes.0))
  rm.idx<-which(predictor.genes.0=="Tissue")
  if (length(rm.idx)>0) predictor.genes.0<-predictor.genes.0[-rm.idx]
  genes.count <- c(genes.count,length(predictor.genes.0))
  predictor.genes.sens <- as.vector(predictor.genes.0)
  
### bagging
  no.of.rounds<-round(length(predictor.genes.sens)/200,digits=0)*3
  if (no.of.rounds<1) no.of.rounds<-4  
  print(paste(drug_sens[i],no.of.rounds,length(predictor.genes.sens),sep=","))
  
  if (length(predictor.genes.sens)>=200) {
    predictor.genes <- matrix(0, no.of.rounds,200)
  } else {
    predictor.genes.s<- sample(predictor.genes.sens,round(length(predictor.genes.sens)*0.75,digits=0))
    predictor.genes <- matrix(0, no.of.rounds,length(predictor.genes.s))
  }
  
  for (b in 1:(no.of.rounds)){
    if (length(predictor.genes.sens)>=200) {
      predictor.genes.smpl<-sample(predictor.genes.sens,200)
    } else {
      predictor.genes.smpl<- sample(predictor.genes.sens,round(length(predictor.genes.sens)*0.75,digits=0))
    }
  predictor.genes[b,]<- predictor.genes.smpl
  }
  
  
  temp.auc<-NA
  temp.auc.counter<-1
  result.cv <- data.frame()
  mean.auc.cv <- vector()
  probs_sensitive_bagging <- data.frame()
  result <- data.frame()
  for (bagging in 1:(no.of.rounds)){
    print(predictor.genes[bagging,])
    pg <-predictor.genes[bagging,]
    
    if (length(pg)>0) {
    rules.train <-train.set[,c(predictor.genes[bagging,],drug_sens[i])]
    train <- rules.train 
    train.t <- rules.train
    
    rules.test<-test.set[,c(predictor.genes[bagging,],drug_sens[i])]
    test <- rules.test

###############  3-fold crossvalidation of training set 
    folds <- cut(seq(1,nrow(train.t)),breaks=3,labels=FALSE)
    
    for (r in 1:3){
      train.cv <- train.t[!(folds==r),]
      test.cv <- train.t[folds==r,]
      
      # preparation of kernels for classification task - gaussian kernel for values and jaccard similarity coef. for binary matrices
      # TRAINING SET
      #should be an Ntra x Ntra x P matrix containing similarity values between training samples
      gene_expr <- as.matrix(train.cv[,grepl("_GE",colnames(train.cv))])
      CNV <- train.cv[,grepl("_CNV",colnames(train.cv))]
      IC50 <- train.cv[,grepl("_IC50",colnames(train.cv))]
      Mut_stat <- as.matrix(train.cv[,!grepl("_IC50",colnames(train.cv))])
      if (!all(is.na(gene_expr))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_GE",colnames(Mut_stat))])}
      if (!all(is.na(CNV))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_CNV",colnames(Mut_stat))])}

      # training set - gene expression
      if (!all(is.na(gene_expr))){
        gene_expr[is.na(gene_expr)] <- 0
        noise <- diag(x=1, nrow(train.cv),nrow(train.cv))
        a <- ncol(gene_expr) 
        gene_expr <- as.matrix(gene_expr)
        gene_expr<- apply(gene_expr,2,as.numeric)
        gene_expr_kernel <- gausskernel(X=gene_expr,sigma=a)
        gene_expr_kernel <- gene_expr_kernel+noise
      }
      
      # training set - CNV GAIN
      if (!all(is.na(CNV))){
        CNV_gain <- as.matrix(CNV)
        CNV_gain[as.character(CNV_gain)=="loss"] <- 0
        CNV_gain[as.character(CNV_gain)==",loss"] <- 0
        CNV_gain[as.character(CNV_gain)=="gain"] <- 1
        CNV_gain[as.character(CNV_gain)==",gain"] <- 1
        CNV_gain[is.na(CNV_gain)] <- 0
        CNV_gain <- as.matrix(CNV_gain)
        CNV_gain <- apply(CNV_gain,2,as.numeric)
        
        jacc_gain<- simil(CNV_gain, method = "Jaccard")
        jacc_gain <- as.matrix(jacc_gain)
        
        diag(jacc_gain) <- 2
        rows.0 <- which(rowSums(CNV_gain)==0)
        if (length(rows.0)>0){
          min.jacc <-min(apply(jacc_gain, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_gain[row,c(which(jacc_gain[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_gain) <- row.names(CNV)
        colnames(jacc_gain) <- row.names(CNV)
        
        # training set - CNV GAIN
        CNV_loss <- as.matrix(CNV)
        CNV_loss[as.character(CNV_loss)=="loss"] <- 1
        CNV_loss[as.character(CNV_loss)==",loss"] <- 1
        CNV_loss[as.character(CNV_loss)=="gain"] <- 0
        CNV_loss[as.character(CNV_loss)==",gain"] <- 0
        CNV_loss[is.na(CNV_loss)] <- 0
        CNV_loss <- as.matrix(CNV_loss)
        CNV_loss <- apply(CNV_loss,2,as.numeric)
        jacc_loss<- simil(CNV_loss, method = "Jaccard")
        jacc_loss <- as.matrix(jacc_loss)
        
        diag(jacc_loss) <- 2
        rows.0 <- which(rowSums(CNV_loss)==0)
        if (length(rows.0)>0){
          min.jacc <-min(apply(jacc_loss, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_loss[row,c(which(jacc_loss[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_loss) <- row.names(CNV)
        colnames(jacc_loss) <- row.names(CNV)
      }
      
      # training set - Mutation status
      if (!all(is.na(Mut_stat))){
        Mut_status <- Mut_stat
        Mut_status[as.character(Mut_status)=="Mut"] <- 1
        Mut_status[is.na(Mut_status)] <- 0
        Mut_status <- as.matrix(Mut_status)
        Mut_status <- apply(Mut_status,2,as.numeric)
        jacc_mut<- simil(Mut_status, method = "Jaccard")
        jacc_mut <- as.matrix(jacc_mut)
        
        diag(jacc_mut) <- 2
        rows.0 <- which(rowSums(Mut_status)==0)
        if (length(rows.0)>0) {
          min.jacc <-min(apply(jacc_mut, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_mut[row,c(which(jacc_mut[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_mut) <- row.names(Mut_stat)
        colnames(jacc_mut) <- row.names(Mut_stat)
      }

      # compile training set composed of constructed kernels
      m <- nrow(train.cv)
      if (exists("gene_expr_kernel") & exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,4))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_gain
        training_set[,,3] <- jacc_loss
        training_set[,,4] <- jacc_mut
      } else if (exists("gene_expr_kernel") & !exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,2))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_mut
      } else if (exists("gene_expr_kernel") & exists("jacc_gain") & !exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,3))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_gain
        training_set[,,3] <- jacc_loss
      }else if (!exists("gene_expr_kernel") & exists("jacc_gain") & !exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,2))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
      }else if (!exists("gene_expr_kernel") & exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,3))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
        training_set[,,3] <- jacc_mut
      }else if (!exists("gene_expr_kernel") & !exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,1))
        training_set[,,1] <- jacc_mut
      }else {
        training_set <- array(,dim=c(m,m,1))
        training_set[,,1] <- gene_expr_kernel
      }
      
#TEST SET
#should be an Ntra x Ntest x P matrix containing similarity values between training and test samples  
      
#training
b <- ncol(test.cv) 
      
      #training
      GE <- as.matrix(train.cv[,grepl("_GE",colnames(train.cv))])
      CNV <- train.cv[,grepl("_CNV",colnames(train.cv))]
      IC50 <- train.cv[,grepl("_IC50",colnames(train.cv))]
      Mut_stat <- as.matrix(train.cv[,!grepl("_IC50",colnames(train.cv))])
      if (!all(is.na(GE))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_GE",colnames(Mut_stat))])}
      if (!all(is.na(CNV))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_CNV",colnames(Mut_stat))]) }
      
      #test
      GE_test <- as.matrix(test.cv[,grepl("_GE",colnames(test.cv))])
      CNV_test <- test.cv[,grepl("_CNV",colnames(test.cv))]
      IC50_test <- as.matrix(test.cv[,grepl("_IC50",colnames(test.cv))])
      colnames(IC50_test) <- grep("_IC50",colnames(test.cv),value = T)
      Mut_stat_test <- as.matrix(test.cv[,!grepl("_IC50",colnames(test.cv))])
      if (!all(is.na(GE_test))){Mut_stat_test <- as.matrix(Mut_stat_test[,!grepl("_GE",colnames(Mut_stat_test))]) }
      if (!all(is.na(CNV_test))){Mut_stat_test <- as.matrix(Mut_stat_test[,!grepl("_CNV",colnames(Mut_stat_test))])}
      
      # test set - gene expression
      if (!all(is.na(GE))){
        gene_expr <- GE
        gene_expr[is.na(gene_expr)] <- 0
        gene_expr <- as.matrix(gene_expr)
        gene_expr<- apply(gene_expr,2,as.numeric)
        row.names(gene_expr) <- row.names(GE)
        gene.row.train <- row.names(GE)
        #test
        gene_expr_test <- GE_test
        gene_expr_test[is.na(gene_expr_test)] <- 0
        gene_expr_test <- as.matrix(gene_expr_test)     
        gene_expr_test<- apply(gene_expr_test,2,as.numeric)
        row.names(gene_expr_test) <- row.names(GE_test)
        gene.col.test <- row.names(GE_test)
        #gaussian kernel
        all_genes<- rbind(gene_expr,gene_expr_test)
        gauss_kernel <- gausskernel(all_genes,sigma=ncol(all_genes))
        gauss_kernel_tr <-gauss_kernel[c(gene.row.train),]
        gauss_kernel_final <- gauss_kernel_tr[,c(gene.col.test)]
      }
      
      # test set -CNV gain
      if (!all(is.na(CNV))){
        CNV_gain <- as.matrix(CNV)
        CNV_gain[as.character(CNV_gain)=="loss"] <- 0
        CNV_gain[as.character(CNV_gain)==",loss"] <- 0
        CNV_gain[as.character(CNV_gain)=="gain"] <- 1
        CNV_gain[as.character(CNV_gain)==",gain"] <- 1
        CNV_gain[is.na(CNV_gain)] <- 0
        CNV_gain <- as.matrix(CNV_gain)
        CNV_gain <- apply(CNV_gain,2,as.numeric)
        
        CNV_gain_test <- as.matrix(CNV_test)
        CNV_gain_test[as.character(CNV_gain_test)=="loss"] <- 0
        CNV_gain_test[as.character(CNV_gain_test)==",loss"] <- 0
        CNV_gain_test[as.character(CNV_gain_test)=="gain"] <- 1
        CNV_gain_test[as.character(CNV_gain_test)==",gain"] <- 1
        CNV_gain_test[is.na(CNV_gain_test)] <- 0
        CNV_gain_test <- as.matrix(CNV_gain_test)
        CNV_gain_test <- apply(CNV_gain_test,2,as.numeric)
        
        jacc_gain<- simil(CNV_gain, CNV_gain_test, method = "Jaccard")
        jacc_gain <- as.matrix(jacc_gain)
        
        rows.0.train <- which(rowSums(CNV_gain)==0)
        rows.0.test <- which(rowSums(CNV_gain_test)==0)
        
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_gain, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_gain[row,rowt] <- min.jacc
            }
          }
        }
        
        # test set - CNV loss
        CNV_loss <- as.matrix(CNV)
        CNV_loss[as.character(CNV_loss)=="loss"] <- 1
        CNV_loss[as.character(CNV_loss)==",loss"] <- 1
        CNV_loss[as.character(CNV_loss)=="gain"] <- 0
        CNV_loss[as.character(CNV_loss)==",gain"] <- 0
        CNV_loss[is.na(CNV_loss)] <- 0
        CNV_loss <- as.matrix(CNV_loss)
        CNV_loss <- apply(CNV_loss,2,as.numeric)
        
        # CNV loss test
        CNV_loss_test <- as.matrix(CNV_test)
        CNV_loss_test[as.character(CNV_loss_test)=="loss"] <- 1
        CNV_loss_test[as.character(CNV_loss_test)==",loss"] <- 1
        CNV_loss_test[as.character(CNV_loss_test)=="gain"] <- 0
        CNV_loss_test[as.character(CNV_loss_test)==",gain"] <- 0
        CNV_loss_test[is.na(CNV_loss_test)] <- 0
        CNV_loss_test <- as.matrix(CNV_loss_test)
        CNV_loss_test <- apply(CNV_loss_test,2,as.numeric)
        
        jacc_loss<- simil(CNV_loss,CNV_loss_test, method = "Jaccard")
        jacc_loss <- as.matrix(jacc_loss)
        
        rows.0.train <- which(rowSums(CNV_loss)==0)
        rows.0.test <- which(rowSums(CNV_loss_test)==0)
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_loss, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_loss[row,rowt] <- min.jacc
            }
          }
        }
        
      }
      
      # test set - Mutation status
      
      if (!all(is.na(Mut_stat))){
        Mut_status <- Mut_stat
        Mut_status[as.character(Mut_status)=="Mut"] <- 1
        Mut_status[is.na(Mut_status)] <- 0
        Mut_status <- as.matrix(Mut_status)
        Mut_status <- apply(Mut_status,2,as.numeric)
        
        #Mut test
        Mut_status_test <- Mut_stat_test
        Mut_status_test[as.character(Mut_status_test)=="Mut"] <- 1
        Mut_status_test[is.na(Mut_status_test)] <- 0
        Mut_status_test <- as.matrix(Mut_status_test)
        Mut_status_test <- apply(Mut_status_test,2,as.numeric)
        
        jacc_mut<- simil(Mut_status,Mut_status_test, method = "Jaccard")
        jacc_mut <- as.matrix(jacc_mut)
        
        rows.0.train <- which(rowSums(Mut_status)==0)
        rows.0.test <- which(rowSums(Mut_status_test)==0)
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_mut, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_mut[row,rowt] <- min.jacc
            }
          }
        }
        
      }

      # compile test set composed of constructed kernels
      n <- nrow(train.cv)
      m <- nrow(test.cv)
      dnames <- rownames(test.cv)
      
      
      if (exists("gauss_kernel_final") & exists("jacc_gain") & exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,4))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_gain
        test_set[,,3] <- jacc_loss
        test_set[,,4] <- jacc_mut
      } else if (exists("gauss_kernel_final") & !exists("jacc_gain") & exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,2))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_mut
      } else if (exists("gauss_kernel_final") & exists("jacc_gain") & !exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,3))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_gain
        test_set[,,3] <- jacc_loss
      } else if (!exists("gauss_kernel_final") & exists("jacc_gain") & !exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,2))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
      }else if (!exists("gauss_kernel_final") & exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,3))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
        training_set[,,3] <- jacc_mut
      }else if (!exists("gauss_kernel_final") & !exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,1))
        training_set[,,1] <- jacc_mut
      } else {
        test_set <- array(,dim=c(n,m,1))
        test_set[,,1] <- gauss_kernel_final
      }
      

      # TRAINING OUTPUT
      #should be an Ntra x 1 matrix containing class labels         
      # IC-50 sensitive
      IC50_sensitive <- as.matrix(IC50)
      IC50_sensitive[as.character(IC50_sensitive)=="Sensitive"] <- 1
      IC50_sensitive[as.character(IC50_sensitive)=="Resistant"] <- (-1)
      IC50_sensitive[is.na(IC50_sensitive)] <- (-1)
      IC50_sensitive <- as.data.frame(IC50_sensitive, stringsAsFactors = FALSE)
      IC50_sensitive <- apply(IC50_sensitive,2,as.numeric)
      row.names(IC50_sensitive) <- row.names(IC50)
      colnames(IC50_sensitive) <- colnames(IC50)
      
      training_output_sensitive <- IC50_sensitive

      # IC-50 resistant
      # IC50_resistant <- as.matrix(IC50)
      # IC50_resistant[as.character(IC50_resistant)=="Resistant"] <- 1
      # IC50_resistant[as.character(IC50_resistant)=="Sensitive"] <- (-1)
      # IC50_resistant[is.na(IC50_resistant)] <- (-1)
      # IC50_resistant <- as.data.frame(IC50_resistant, stringsAsFactors = FALSE)
      # IC50_resistant <- apply(IC50_resistant,2,as.numeric)
      # row.names(IC50_resistant) <- row.names(IC50)
      # colnames(IC50_resistant) <- colnames(IC50)
      # 
      # 
      # training_output_resistant <- IC50_resistant
      # save(training_output_resistant, file=paste(f,"_training_output_resistant.RData",sep = ""))
      
      ########## test labels for ROC curves
      
      # IC-50 sensitive
      IC50_sensitive <- as.matrix(IC50_test)
      IC50_sensitive[as.character(IC50_sensitive)=="Sensitive"] <- 1
      IC50_sensitive[as.character(IC50_sensitive)=="Resistant"] <- (-1)
      IC50_sensitive[is.na(IC50_sensitive)] <- (-1)
      IC50_sensitive_ROC_labels <- as.data.frame(IC50_sensitive, stringsAsFactors = FALSE)
      IC50_sensitive_ROC_labels <- apply(IC50_sensitive_ROC_labels,2,as.numeric)
      row.names(IC50_sensitive_ROC_labels) <- row.names(IC50_test)
      colnames(IC50_sensitive_ROC_labels) <- colnames(IC50_test)
      
      # IC-50 resistant
      # IC50_resistant <- as.matrix(IC50_test)
      # IC50_resistant[as.character(IC50_resistant)=="Resistant"] <- 1
      # IC50_resistant[as.character(IC50_resistant)=="Sensitive"] <- (-1)
      # IC50_resistant[is.na(IC50_resistant)] <- (-1)
      # IC50_resistant_ROC_labels <- as.data.frame(IC50_resistant, stringsAsFactors = FALSE)
      # IC50_resistant_ROC_labels <- apply(IC50_resistant_ROC_labels,2,as.numeric)
      # row.names(IC50_resistant_ROC_labels) <- row.names(IC50_test)
      # colnames(IC50_resistant_ROC_labels) <- colnames(IC50_test)
      # save(IC50_resistant_ROC_labels, file=paste(f,"_IC50_resistant_ROC_labels.RData",sep = ""))
      
      remove(gene_expr_kernel)
      remove(gauss_kernel_final)
      remove(jacc_gain)
      remove(jacc_loss)
      remove(jacc_mut)
      
      # CLASSIFICATION
      
      T<- ncol(training_output_sensitive)
      N<- ncol(test_set[,,1])
      d <- dim(training_set)[3]
      weights <- matrix(0,T,4)
      probs_sensitive <- matrix(0,N,T)
      
      
      #initalize the parameters of the algorithm
      parameters <- list()
      
      #set the hyperparameters of gamma prior used for sample weights
      parameters$alpha_lambda <- 1
      parameters$beta_lambda <- 1
      
      #set the hyperparameters of gamma prior used for bias
      parameters$alpha_gamma <- 1
      parameters$beta_gamma <- 1
      
      #set the hyperparameters of gamma prior used for kernel weights
      parameters$alpha_omega <- 1
      parameters$beta_omega <- 1
      
      ### IMPORTANT ###
      #For gamma priors, you can experiment with three different (alpha, beta) values
      #(1, 1) => default priors
      #(1e-10, 1e+10) => good for obtaining sparsity
      #(1e-10, 1e-10) => good for small sample size problems
      
      #set the number of iterations
      parameters$iteration <- 200
      
      #set the margin parameter
      parameters$margin <- 1
      
      #determine whether you want to store the lower bound values
      parameters$progress <- 0
      
      #set the seed for random number generator used to initalize random variables
      parameters$seed <- 1606
      
      #set the standard deviation of intermediate representations
      parameters$sigma_g <- 0.1
      
      #initialize the kernels and class labels for training
      Ktrain <- training_set #should be an Ntra x Ntra x P matrix containing similarity values between training samples
      ytrain <- as.numeric(training_output_sensitive) #should be an Ntra x 1 matrix containing class labels (contains only -1s and +1s)
      
      #perform training
      state <- bemkl_supervised_classification_variational_train(Ktrain, ytrain, parameters)
      
      #display the kernel weights
      # print(state$be$mu[-1])
      weights <- state$be$mu[-1]
      #save(weights, file=paste(bagging,"_",drug_sens_files[i],"_kernel_weights_sensitive.RData",sep=""))
      
      #initialize the kernels for testing
      Ktest <- test_set #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples
      
      #perform prediction
      prediction <- bemkl_supervised_classification_variational_test(Ktest, state)
      
      #display the predicted probabilities
      #print(prediction$p)
      probs_sensitive <- as.matrix(prediction$p)
      # save(probs_sensitive, file=paste(drug_sens[i],"_probs_sensitive.RData",sep= ""))
      
      
      IC50_sensitive_ROC_labels <-as.matrix(IC50_sensitive_ROC_labels)
      pred <- prediction(probs_sensitive[,1],IC50_sensitive_ROC_labels[,1])
      result.cv[r,"Drug"] <- paste(bagging,"_",r,"_",drug_sens_files[i],sep="")
      result.cv[r,"AUC"] <- unlist(performance(pred,"auc")@y.values)
      
    }
    #save(result.cv,file = paste(drug_sens_files[i],"_AUCs_crossvalidation.RData",sep = ""))
    print(result.cv)
    mean.auc.cv[bagging] <- mean(result.cv$AUC)
    print(mean.auc.cv)
    
    
###############  end of  3-fold crossvalidation of training set 
    
    
############### construction of test and training set for classification
    # TRAINING SET
    #should be an Ntra x Ntra x P matrix containing similarity values between training samples
      gene_expr <- as.matrix(train[,grepl("_GE",colnames(train))])
      CNV <- train[,grepl("_CNV",colnames(train))]
      IC50 <- train[,grepl("_IC50",colnames(train))]
      Mut_stat <- as.matrix(train[,!grepl("_IC50",colnames(train))])
      if (!all(is.na(gene_expr))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_GE",colnames(Mut_stat))])}
      if (!all(is.na(CNV))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_CNV",colnames(Mut_stat))])}
      #Mut_stat <- Mut_stat[,!grepl("Tissue",colnames(Mut_stat))]
      
      # training set - gene expression
      if (!all(is.na(gene_expr))){
      gene_expr[is.na(gene_expr)] <- 0
      noise <- diag(x=1, nrow(train),nrow(train))
      a <- ncol(gene_expr) 
      gene_expr <- as.matrix(gene_expr)
      gene_expr<- apply(gene_expr,2,as.numeric)
      gene_expr_kernel <- gausskernel(X=gene_expr,sigma=a)
      gene_expr_kernel <- gene_expr_kernel+noise
      }
      
    # training set - CNV GAIN
      if (!all(is.na(CNV))){
        CNV_gain <- as.matrix(CNV)
        CNV_gain[as.character(CNV_gain)=="loss"] <- 0
        CNV_gain[as.character(CNV_gain)==",loss"] <- 0
        CNV_gain[as.character(CNV_gain)=="gain"] <- 1
        CNV_gain[as.character(CNV_gain)==",gain"] <- 1
        CNV_gain[is.na(CNV_gain)] <- 0
        CNV_gain <- as.matrix(CNV_gain)
        CNV_gain <- apply(CNV_gain,2,as.numeric)
        
        jacc_gain<- simil(CNV_gain, method = "Jaccard")
        jacc_gain <- as.matrix(jacc_gain)
        
        diag(jacc_gain) <- 2
        rows.0 <- which(rowSums(CNV_gain)==0)
        if (length(rows.0)>0){
          min.jacc <-min(apply(jacc_gain, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_gain[row,c(which(jacc_gain[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_gain) <- row.names(CNV)
        colnames(jacc_gain) <- row.names(CNV)
        
 # training set - CNV GAIN
        CNV_loss <- as.matrix(CNV)
        CNV_loss[as.character(CNV_loss)=="loss"] <- 1
        CNV_loss[as.character(CNV_loss)==",loss"] <- 1
        CNV_loss[as.character(CNV_loss)=="gain"] <- 0
        CNV_loss[as.character(CNV_loss)==",gain"] <- 0
        CNV_loss[is.na(CNV_loss)] <- 0
        CNV_loss <- as.matrix(CNV_loss)
        CNV_loss <- apply(CNV_loss,2,as.numeric)
        jacc_loss<- simil(CNV_loss, method = "Jaccard")
        jacc_loss <- as.matrix(jacc_loss)
        
        diag(jacc_loss) <- 2
        rows.0 <- which(rowSums(CNV_loss)==0)
        if (length(rows.0)>0){
          min.jacc <-min(apply(jacc_loss, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_loss[row,c(which(jacc_loss[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_loss) <- row.names(CNV)
        colnames(jacc_loss) <- row.names(CNV)
      }
      
      # training set - Mutation status
      if (!all(is.na(Mut_stat))){
        Mut_status <- Mut_stat
        Mut_status[as.character(Mut_status)=="Mut"] <- 1
        Mut_status[is.na(Mut_status)] <- 0
        Mut_status <- as.matrix(Mut_status)
        Mut_status <- apply(Mut_status,2,as.numeric)
        jacc_mut<- simil(Mut_status, method = "Jaccard")
        jacc_mut <- as.matrix(jacc_mut)
        
        diag(jacc_mut) <- 2
        rows.0 <- which(rowSums(Mut_status)==0)
        if (length(rows.0)>0) {
          min.jacc <-min(apply(jacc_mut, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0){
            jacc_mut[row,c(which(jacc_mut[row,]==1))] <- min.jacc
          }
        }
        row.names(jacc_mut) <- row.names(Mut_stat)
        colnames(jacc_mut) <- row.names(Mut_stat)
      }

      
      # compile training set composed of constructed kernels     
      m <- nrow(train)
      if (exists("gene_expr_kernel") & exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,4))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_gain
        training_set[,,3] <- jacc_loss
        training_set[,,4] <- jacc_mut
      } else if (exists("gene_expr_kernel") & !exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,2))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_mut
      } else if (exists("gene_expr_kernel") & exists("jacc_gain") & !exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,3))
        training_set[,,1] <- gene_expr_kernel
        training_set[,,2] <- jacc_gain
        training_set[,,3] <- jacc_loss
      }else if (!exists("gene_expr_kernel") & exists("jacc_gain") & !exists("jacc_mut")){
      training_set <- array(,dim=c(m,m,2))
      training_set[,,1] <- jacc_gain
      training_set[,,2] <- jacc_loss
      }else if (!exists("gene_expr_kernel") & exists("jacc_gain") & exists("jacc_mut")){
      training_set <- array(,dim=c(m,m,3))
      training_set[,,1] <- jacc_gain
      training_set[,,2] <- jacc_loss
      training_set[,,3] <- jacc_mut
    }else if (!exists("gene_expr_kernel") & !exists("jacc_gain") & exists("jacc_mut")){
      training_set <- array(,dim=c(m,m,1))
      training_set[,,1] <- jacc_mut
    }else {
        training_set <- array(,dim=c(m,m,1))
        training_set[,,1] <- gene_expr_kernel
      }

      
#TEST SET
#should be an Ntra x Ntest x P matrix containing similarity values between training and test samples      

      #training
      GE <- as.matrix(train[,grepl("_GE",colnames(train))])
      CNV <- train[,grepl("_CNV",colnames(train))]
      IC50 <- train[,grepl("_IC50",colnames(train))]
      Mut_stat <- as.matrix(train[,!grepl("_IC50",colnames(train))])
      if (!all(is.na(GE))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_GE",colnames(Mut_stat))])}
      if (!all(is.na(CNV))){Mut_stat <- as.matrix(Mut_stat[,!grepl("_CNV",colnames(Mut_stat))]) }
      #Mut_stat <- Mut_stat[,!grepl("Tissue",colnames(Mut_stat))]
      
      #test
      GE_test <- as.matrix(test[,grepl("_GE",colnames(test))])
      CNV_test <- test[,grepl("_CNV",colnames(test))]
      IC50_test <- as.matrix(test[,grepl("_IC50",colnames(test))])
      colnames(IC50_test) <- grep("_IC50",colnames(test),value = T)
      Mut_stat_test <- as.matrix(test[,!grepl("_IC50",colnames(test))])
      if (!all(is.na(GE_test))){Mut_stat_test <- as.matrix(Mut_stat_test[,!grepl("_GE",colnames(Mut_stat_test))]) }
      if (!all(is.na(CNV_test))){Mut_stat_test <- as.matrix(Mut_stat_test[,!grepl("_CNV",colnames(Mut_stat_test))])}
      #Mut_stat_test <- Mut_stat_test[,!grepl("Tissue",colnames(Mut_stat_test))]
      
      
      # test set - gene expression
      if (!all(is.na(GE))){
      gene_expr <- GE
      gene_expr[is.na(gene_expr)] <- 0
      gene_expr <- as.matrix(gene_expr)
      gene_expr<- apply(gene_expr,2,as.numeric)
      row.names(gene_expr) <- row.names(GE)
      gene.row.train <- row.names(GE)
      #test
      gene_expr_test <- GE_test
      gene_expr_test[is.na(gene_expr_test)] <- 0
      gene_expr_test <- as.matrix(gene_expr_test)     
      gene_expr_test<- apply(gene_expr_test,2,as.numeric)
      row.names(gene_expr_test) <- row.names(GE_test)
      gene.col.test <- row.names(GE_test)
      #gaussian kernel
      all_genes<- rbind(gene_expr,gene_expr_test)
      gauss_kernel <- gausskernel(all_genes,sigma=ncol(all_genes))
      gauss_kernel_tr <-gauss_kernel[c(gene.row.train),]
      gauss_kernel_final <- gauss_kernel_tr[,c(gene.col.test)]
      }
      
      
      # test set -CNV gain
      if (!all(is.na(CNV))){
        CNV_gain <- as.matrix(CNV)
        CNV_gain[as.character(CNV_gain)=="loss"] <- 0
        CNV_gain[as.character(CNV_gain)==",loss"] <- 0
        CNV_gain[as.character(CNV_gain)=="gain"] <- 1
        CNV_gain[as.character(CNV_gain)==",gain"] <- 1
        CNV_gain[is.na(CNV_gain)] <- 0
        CNV_gain <- as.matrix(CNV_gain)
        CNV_gain <- apply(CNV_gain,2,as.numeric)
        
        CNV_gain_test <- as.matrix(CNV_test)
        CNV_gain_test[as.character(CNV_gain_test)=="loss"] <- 0
        CNV_gain_test[as.character(CNV_gain_test)==",loss"] <- 0
        CNV_gain_test[as.character(CNV_gain_test)=="gain"] <- 1
        CNV_gain_test[as.character(CNV_gain_test)==",gain"] <- 1
        CNV_gain_test[is.na(CNV_gain_test)] <- 0
        CNV_gain_test <- as.matrix(CNV_gain_test)
        CNV_gain_test <- apply(CNV_gain_test,2,as.numeric)
        
        jacc_gain<- simil(CNV_gain, CNV_gain_test, method = "Jaccard")
        jacc_gain <- as.matrix(jacc_gain)
        
        rows.0.train <- which(rowSums(CNV_gain)==0)
        rows.0.test <- which(rowSums(CNV_gain_test)==0)
        
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_gain, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_gain[row,rowt] <- min.jacc
            }
          }
        }
        
        
        # test set - CNV loss
        CNV_loss <- as.matrix(CNV)
        CNV_loss[as.character(CNV_loss)=="loss"] <- 1
        CNV_loss[as.character(CNV_loss)==",loss"] <- 1
        CNV_loss[as.character(CNV_loss)=="gain"] <- 0
        CNV_loss[as.character(CNV_loss)==",gain"] <- 0
        CNV_loss[is.na(CNV_loss)] <- 0
        CNV_loss <- as.matrix(CNV_loss)
        CNV_loss <- apply(CNV_loss,2,as.numeric)
        # CNV loss test
        CNV_loss_test <- as.matrix(CNV_test)
        CNV_loss_test[as.character(CNV_loss_test)=="loss"] <- 1
        CNV_loss_test[as.character(CNV_loss_test)==",loss"] <- 1
        CNV_loss_test[as.character(CNV_loss_test)=="gain"] <- 0
        CNV_loss_test[as.character(CNV_loss_test)==",gain"] <- 0
        CNV_loss_test[is.na(CNV_loss_test)] <- 0
        CNV_loss_test <- as.matrix(CNV_loss_test)
        CNV_loss_test <- apply(CNV_loss_test,2,as.numeric)
        
        
        jacc_loss<- simil(CNV_loss,CNV_loss_test, method = "Jaccard")
        jacc_loss <- as.matrix(jacc_loss)
        
        
        rows.0.train <- which(rowSums(CNV_loss)==0)
        rows.0.test <- which(rowSums(CNV_loss_test)==0)
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_loss, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_loss[row,rowt] <- min.jacc
            }
          }
        }
        
      } 
      
      # test set - Mutation status
      
      if (!all(is.na(Mut_stat))){
        Mut_status <- Mut_stat
        Mut_status[as.character(Mut_status)=="Mut"] <- 1
        Mut_status[is.na(Mut_status)] <- 0
        Mut_status <- as.matrix(Mut_status)
        Mut_status <- apply(Mut_status,2,as.numeric)
        #Mut test
        Mut_status_test <- Mut_stat_test
        Mut_status_test[as.character(Mut_status_test)=="Mut"] <- 1
        Mut_status_test[is.na(Mut_status_test)] <- 0
        Mut_status_test <- as.matrix(Mut_status_test)
        Mut_status_test <- apply(Mut_status_test,2,as.numeric)
        
        jacc_mut<- simil(Mut_status,Mut_status_test, method = "Jaccard")
        jacc_mut <- as.matrix(jacc_mut)
        
        rows.0.train <- which(rowSums(Mut_status)==0)
        rows.0.test <- which(rowSums(Mut_status_test)==0)
        if (length(rows.0.train)>0 & length(rows.0.test)>0){
          min.jacc <-min(apply(jacc_mut, 1, FUN = function(x) {min(x[x > 0])}))
          for (row in rows.0.train){
            for (rowt in rows.0.test){
              jacc_mut[row,rowt] <- min.jacc
            }
          }
        }
        
      } 
      

# compile test set composed of constructed kernels    
      ### combine all into training matrix
      n <- nrow(train)
      m <- nrow(test)
      dnames <- rownames(test)
      
      
      if (exists("gauss_kernel_final") & exists("jacc_gain") & exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,4))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_gain
        test_set[,,3] <- jacc_loss
        test_set[,,4] <- jacc_mut
      } else if (exists("gauss_kernel_final") & !exists("jacc_gain") & exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,2))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_mut
        
      } else if (exists("gauss_kernel_final") & exists("jacc_gain") & !exists("jacc_mut")){
        test_set <- array(,dim=c(n,m,3))
        test_set[,,1] <- gauss_kernel_final
        test_set[,,2] <- jacc_gain
        test_set[,,3] <- jacc_loss
      } else if (!exists("gauss_kernel_final") & exists("jacc_gain") & !exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,2))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
      }else if (!exists("gauss_kernel_final") & exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,3))
        training_set[,,1] <- jacc_gain
        training_set[,,2] <- jacc_loss
        training_set[,,3] <- jacc_mut
      }else if (!exists("gauss_kernel_final") & !exists("jacc_gain") & exists("jacc_mut")){
        training_set <- array(,dim=c(m,m,1))
        training_set[,,1] <- jacc_mut
      } else {
        test_set <- array(,dim=c(n,m,1))
        test_set[,,1] <- gauss_kernel_final
      }

      # TRAINING OUTPUT
      #should be an Ntra x 1 matrix containing class labels   
      
      # # IC-50 sensitive
      IC50_sensitive <- as.matrix(IC50)
      IC50_sensitive[as.character(IC50_sensitive)=="Sensitive"] <- 1
      IC50_sensitive[as.character(IC50_sensitive)=="Resistant"] <- (-1)
      IC50_sensitive[is.na(IC50_sensitive)] <- (-1)
      IC50_sensitive <- as.data.frame(IC50_sensitive, stringsAsFactors = FALSE)
      IC50_sensitive <- apply(IC50_sensitive,2,as.numeric)
      row.names(IC50_sensitive) <- row.names(IC50)
      colnames(IC50_sensitive) <- colnames(IC50)
      
      training_output_sensitive <- IC50_sensitive

    
      # IC-50 resistant
      # IC50_resistant <- as.matrix(IC50)
      # IC50_resistant[as.character(IC50_resistant)=="Resistant"] <- 1
      # IC50_resistant[as.character(IC50_resistant)=="Sensitive"] <- (-1)
      # IC50_resistant[is.na(IC50_resistant)] <- (-1)
      # IC50_resistant <- as.data.frame(IC50_resistant, stringsAsFactors = FALSE)
      # IC50_resistant <- apply(IC50_resistant,2,as.numeric)
      # row.names(IC50_resistant) <- row.names(IC50)
      # colnames(IC50_resistant) <- colnames(IC50)
      # 
      # 
      # training_output_resistant <- IC50_resistant
      
      ########## test labels for ROC curves
      
      # IC-50 sensitive
      IC50_sensitive <- as.matrix(IC50_test)
      IC50_sensitive[as.character(IC50_sensitive)=="Sensitive"] <- 1
      IC50_sensitive[as.character(IC50_sensitive)=="Resistant"] <- (-1)
      IC50_sensitive[is.na(IC50_sensitive)] <- (-1)
      IC50_sensitive_ROC_labels <- as.data.frame(IC50_sensitive, stringsAsFactors = FALSE)
      IC50_sensitive_ROC_labels <- apply(IC50_sensitive_ROC_labels,2,as.numeric)
      row.names(IC50_sensitive_ROC_labels) <- row.names(IC50_test)
      colnames(IC50_sensitive_ROC_labels) <- colnames(IC50_test)
      # save(IC50_sensitive_ROC_labels, file=paste(drug_res[f],"_IC50_sensitive_ROC_labels.RData",sep=""))
      
      # IC-50 resistant
      # IC50_resistant <- as.matrix(IC50_test)
      # IC50_resistant[as.character(IC50_resistant)=="Resistant"] <- 1
      # IC50_resistant[as.character(IC50_resistant)=="Sensitive"] <- (-1)
      # IC50_resistant[is.na(IC50_resistant)] <- (-1)
      # IC50_resistant_ROC_labels <- as.data.frame(IC50_resistant, stringsAsFactors = FALSE)
      # IC50_resistant_ROC_labels <- apply(IC50_resistant_ROC_labels,2,as.numeric)
      # row.names(IC50_resistant_ROC_labels) <- row.names(IC50_test)
      # colnames(IC50_resistant_ROC_labels) <- colnames(IC50_test)


    remove(gene_expr_kernel)
    remove(gauss_kernel_final)
    remove(jacc_gain)
    remove(jacc_loss)
    remove(jacc_mut)
    
    # CLASSIFICATION
    
    T<- ncol(training_output_sensitive)
    N<- ncol(test_set[,,1])
    d <- dim(training_set)[3]
    weights <- matrix(0,T,4)
    probs_sensitive <- matrix(0,N,T)
    
      
      #initalize the parameters of the algorithm
      parameters <- list()
      
      #set the hyperparameters of gamma prior used for sample weights
      parameters$alpha_lambda <- 1
      parameters$beta_lambda <- 1
      
      #set the hyperparameters of gamma prior used for bias
      parameters$alpha_gamma <- 1
      parameters$beta_gamma <- 1
      
      #set the hyperparameters of gamma prior used for kernel weights
      parameters$alpha_omega <- 1
      parameters$beta_omega <- 1
      
      ### IMPORTANT ###
      #For gamma priors, you can experiment with three different (alpha, beta) values
      #(1, 1) => default priors
      #(1e-10, 1e+10) => good for obtaining sparsity
      #(1e-10, 1e-10) => good for small sample size problems
      
      #set the number of iterations
      parameters$iteration <- 200
      
      #set the margin parameter
      parameters$margin <- 1
      
      #determine whether you want to store the lower bound values
      parameters$progress <- 0
      
      #set the seed for random number generator used to initalize random variables
      parameters$seed <- 1606
      
      #set the standard deviation of intermediate representations
      parameters$sigma_g <- 0.1
      
      #initialize the kernels and class labels for training
      Ktrain <- training_set #should be an Ntra x Ntra x P matrix containing similarity values between training samples
      ytrain <- as.numeric(training_output_sensitive) #should be an Ntra x 1 matrix containing class labels (contains only -1s and +1s)
      
      #perform training
      state <- bemkl_supervised_classification_variational_train(Ktrain, ytrain, parameters)
      
      #display the kernel weights
     # print(state$be$mu[-1])
      weights <- state$be$mu[-1]
      #save(weights, file=paste(bagging,"_",drug_sens_files[i],"_kernel_weights_sensitive.RData",sep=""))
      
      #initialize the kernels for testing
      Ktest <- test_set #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples
      
      #perform prediction
      prediction <- bemkl_supervised_classification_variational_test(Ktest, state)
      
      #display the predicted probabilities
      #print(prediction$p)
      probs_sensitive <- as.matrix(prediction$p)

      if (bagging==1) {probs_sensitive_bagging <- probs_sensitive[,1]} else { probs_sensitive_bagging <- cbind(probs_sensitive_bagging,probs_sensitive[,1])}

    }
  }
  
  save(probs_sensitive_bagging, file=paste(bagging,"_",drug_sens_files[i],"_probs_sensitive.RData",sep= ""))
  save(mean.auc.cv, file=paste(bagging,"_",drug_sens_files[i],"_mean_auc_cv.RData",sep= ""))
  print(probs_sensitive_bagging)
  print(mean.auc.cv)
  
  #### weighted mean
  if (ncol(probs_sensitive_bagging)>0){
  p.final<-apply(probs_sensitive_bagging,1,weighted.mean,mean.auc.cv)
  IC50_sensitive_ROC_labels <-as.matrix(IC50_sensitive_ROC_labels)
  pred <- prediction(p.final,IC50_sensitive_ROC_labels[,1])
  auc_final_bagging[i,"Drug"] <- paste(drug_sens[i],sep="")
  auc_final_bagging[i,"AUC"] <- unlist(performance(pred,"auc")@y.values)
  save(auc_final_bagging,file = "auc_final_bagging_temp.RData")
  }
}

save(auc_final_bagging,file = "auc_final_bagging.RData")
write.table(auc_final_bagging, file="auc_final_bagging.txt", sep="\t")

