#############################################
###############Whole DataSet#################
#############################################

require(arules)

load("~/Projects/MBA/new_data/RData/MASTER_MATRIX_CATEGORICAL.RData")
load("~/Projects/MBA/new_data/RData/MASTER_MATRIX_PERMUTED.RData")

b.t<-as(as.data.frame(m.final),"transactions")
b.t.perm<-as(as.data.frame(m.final.perm),"transactions")
drug.idx<-grep("IC50",labels(m.final)[[2]])
drugs<-labels(m.final)[[2]][drug.idx]
drugs.total<-c(paste(drugs,"Sensitive",sep="="),paste(drugs,"Resistant",sep="="))

rules<-apriori(b.t,parameter = list(support = (3/length(b.t)), confidence = 3/length(b.t), minlen=2, maxlen=2),appearance = list(rhs=drugs.total,default="lhs"))



temp.quality<-quality(rules)
temp.quality[,3]<-NA
temp.quality<-temp.quality[!duplicated(temp.quality), ]
gc()

for (i in 1:nrow(temp.quality)){
  print(i)
  rules.perm<-apriori(b.t.perm,parameter = list(support = temp.quality$support[i], confidence = temp.quality$confidence[i], minlen=2, maxlen=2))
  if (length(rules.perm)==0) next
  temp<-hist(log(quality(rules.perm)[,3]),breaks=80,plot=F)
  lift.estim.0.05<-exp(temp$breaks[which(cumsum(temp$counts)>sum(temp$counts)*0.95)[1]])
  temp.quality$lift[i]<-lift.estim.0.05
  rm(rules.perm)
  gc()
}
temp.quality$lift[temp.quality$lift<=1]<-NA
temp.quality$lift[is.na(temp.quality$lift)]<-min(temp.quality$lift,na.rm=T)


gc()


status<-rep(NA,length(rules))
quality.rules<-quality(rules)
require(arules)
for (i in 1:nrow(temp.quality)){
  print(i)
  temp.idx<-which(quality.rules$support==temp.quality$support[i] & quality.rules$confidence==temp.quality$confidence[i])
  status[temp.idx]<-sapply (quality.rules$lift[temp.idx], function(lift) if(lift>temp.quality$lift[i]){temp.out<-"pass"}else{temp.out<-"fail"})
}
gc()

rules.sig<-rules[which(status=="pass")]
save(rules.sig,file="~/Projects/MBA/new_data/RData/RULES_SIGNIFICANT_Sup3_Conf3in1001_DYNAMIC_THRESHOLD_FDR0.05.RData")
rm(rules)
gc()

