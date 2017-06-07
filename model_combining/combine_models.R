# we import the tuned models and training data again...

load("chrom.models.Rdata")
load("seq.models.Rdata")
load("ML_data_summary_train_scaled.Rdata")
load("clustp_scores_train_scaled.Rdata")
load("velements.Rdata")

require(tunetest)
require(data.table)
require(foreach)
require(doMC)
require(reshape2)
require(stringr)
require(e1071)
require(randomForest)
require(glmnet)
require(caret)
require(PRROC)
require(plyr)
require(ggplot2)
require(RColorBrewer)

# predicting the models on their 'leave one out' test sets:
# i.e. call to tunetest::predict()
predictions.chromatin<-lapply( tunedmodels.chrom, function(model){ lapply(model, function(feature){lapply(feature, function(fold){predict(fold, type="prob")} )} )} )
predictions.seq<-lapply( tunedmodels.seq, function(model){ lapply(model, function(feature){lapply(feature, function(fold){predict(fold, type="prob")} )} )} )

# combining the predictions on the test sets...
preds.chromatin<-lapply( predictions.chromatin, function(model){ lapply(model, function(feature){unlist( lapply(feature, function(fold){fold$test.pred[,1]}) ) } )} )
preds.seq<-lapply( predictions.seq, function(model){ lapply(model, function(feature){ unlist( lapply(feature, function(fold){fold$test.pred[,1]}) ) } )} )

# merging into one big data.frame
preds<-data.frame( "chromatin.LASSO"=preds.chromatin$lasso$all, "chromatin.SVMlin"=preds.chromatin$svm.lin$all, "chromatin.SVMrad"=preds.chromatin$svm.rad$all, "chromatin.RF"=preds.chromatin$rf$all, "sequence.LASSO"=preds.seq$lasso$all, "sequence.SVMlin"=preds.seq$svm.lin$all, "sequence.RF"=preds.seq$rf$all)
preds<-data.frame(apply(preds,2,unlist))

# assign fold information
# 'n' simply contains the number of observations per fold
n<-as.vector( lapply( predictions.chromatin, function(model){ lapply(model, function(feature){do.call("rbind", lapply(feature, function(fold){nrow(fold$test.pred)}) ) } )} )$rf$all )
preds$n<-rep(1:10,n)
preds$y<-unlist( lapply(predictions.chromatin$lasso$all, function(x){x$test.pred[,2]}))

# forming all possible 9 v 1 splits :
combs<-combn(x = 10, 9)

# The function below acts on the columns of 'combs'
# it combines the corresponding folds using Ridge regularized regression and the sum of ranks
# performances for all models are returned
# it operates on the objects generated above and the indexes are hard-coded (so watch out what you are doing...)
set.seed(80085)

fit_ridge<-function(i){
  
  data.train<-preds[ preds$n %in% combs[,i], ]
  # remove the 'n' column, subset to the folds of interest:
  data.test<-preds[ !( preds$n %in% combs[,i]) , c(1:7,9) ]
  
  y <- data.train$y
  x <- as.matrix(data.train[,1:7])
  
  # setting observations weights for the Ridge-GLM (balance weights for positive and negative classes)
  obsweights<-rep(1, length(y))
  tab<-table(y)
  wgt<-as.numeric(tab["negative"]/tab["positive"])
  obsweights[ y=="positive" ]<-wgt
  
  # fitting ridge regularized models
  # remember to set "nfolds" to a sensible value depending on the size of the set...
  model.chromatin<-cv.glmnet(x[,1:4],y,type.measure="mse",nfolds=10,alpha=0,family="binomial",weights = obsweights)
  model.sequence<-cv.glmnet(x[,5:7],y,type.measure="mse",nfolds=10,alpha=0,family="binomial", weights = obsweights)
  model.mixed<-cv.glmnet(x[,1:7],y,type.measure="mse",nfolds=10,alpha=0,family="binomial", weights = obsweights)
  
  # this function retrieves the original auprcs and aurocs of single models (see below)
  single.perfs<-get_perf.test(data.test)
  
  pred.chromatin<-as.vector( predict(model.chromatin, newx=as.matrix(data.test[,-c(5:9)]), type="response",s="lambda.1se") )
  pred.sequence<-as.vector( predict(model.sequence, newx=as.matrix(data.test[,-c(1:4,8:ncol(data.test))]), type="response",s="lambda.1se") )
  pred.mixed<-as.vector( predict(model.mixed, newx=as.matrix(data.test[,-ncol(data.test)]), type="response",s="lambda.1se") )
  
  # computing the sum of ranks...
  pred.sor.sequence<-apply( apply( data.test[,-c(1:4,8:ncol(data.test))], 2, rank ), 1, sum )
  pred.sor.chromatin<-apply( apply( data.test[,-c(5:9)], 2, rank ), 1, sum )
  pred.sor.combined<-pred.sor.sequence*0.3+pred.sor.chromatin
  
  preds_all<-data.frame("sequence.Ridge"=pred.sequence,"chromatin.Ridge"=pred.chromatin,"combined.Ridge"=pred.mixed,"sequence.SOR"=pred.sor.sequence,"chromatin.SOR"=pred.sor.chromatin,"combined.SOR"=pred.sor.combined,"y"=preds$y[ !(preds$n %in% combs[,i]) ])
  
  ind.pos<-data.test$y=="positive"
  ind.neg<-data.test$y=="negative"
  
  # pr.* contain the pr-curves for the ridge models
  pr.chromatin<-pr.curve(pred.chromatin[ind.pos],pred.chromatin[ind.neg],curve=T)
  pr.sequence<-pr.curve(pred.sequence[ind.pos],pred.sequence[ind.neg],curve=T)
  pr.mixed<-pr.curve(pred.mixed[ind.pos],pred.mixed[ind.neg],curve=T)
  # pr.sor.* contain pr-curves for the sum of ranks
  pr.sor.sequence<-pr.curve(pred.sor.sequence[ind.pos],pred.sor.sequence[ind.neg],curve=T)
  pr.sor.chromatin<-pr.curve(pred.sor.chromatin[ind.pos],pred.sor.chromatin[ind.neg],curve=T)
  pr.sor.combined<-pr.curve(pred.sor.combined[ind.pos],pred.sor.combined[ind.neg],curve=T)
  
  # roc.* contain the roc-curves for the ridge models  
  roc.chromatin<-roc.curve(pred.chromatin[ind.pos],pred.chromatin[ind.neg],curve=T)
  roc.sequence<-roc.curve(pred.sequence[ind.pos],pred.sequence[ind.neg],curve=T)
  roc.mixed<-roc.curve(pred.mixed[ind.pos],pred.mixed[ind.neg],curve=T)
  # roc.sor.* contain roc-curves for the sum of ranks
  roc.sor.sequence<-roc.curve(pred.sor.sequence[ind.pos],pred.sor.sequence[ind.neg],curve=T)
  roc.sor.chromatin<-roc.curve(pred.sor.chromatin[ind.pos],pred.sor.chromatin[ind.neg],curve=T)
  roc.sor.combined<-roc.curve(pred.sor.combined[ind.pos],pred.sor.combined[ind.neg],curve=T)
  
  # getting coefficients ("weights")
  coefs<-coef(model.mixed, "lambda.1se")
  
  combined.res<-matrix(c(pr.chromatin$auc.integral, roc.chromatin$auc, pr.sequence$auc.integral, roc.sequence$auc, pr.sor.sequence$auc.integral, roc.sor.sequence$auc, pr.sor.chromatin$auc.integral, roc.sor.chromatin$auc, pr.sor.combined$auc.integral, roc.sor.combined$auc, pr.mixed$auc.integral, roc.mixed$auc), nrow=2 )
  colnames(combined.res)<-c("chromatin.comb.Ridge","sequence.comb.Ridge","sequence.comb.SOR","chromatin.comb.SOR","seq.chr.comb.SOR","seq.chr.comb.Ridge")
  
  result<-list( "performances"=cbind( single.perfs, combined.res ), "coef.mixed"=coefs, "model"=model.mixed, "roc"=roc.mixed, "pr"=pr.mixed, "single.combined"=list("chr"=list("pr"=pr.chromatin, "roc"=roc.chromatin), "seq"=list("pr"=pr.sequence,"roc"=roc.sequence)), "predictions"=preds_all)
  
  return(result)
  
}

get_perf.test<-function(data.test){
  
  ind.pos<-data.test$y=="positive"
  ind.neg<-!ind.pos
  
  prs<-sapply(data.test[,-ncol(data.test)], function(x){ pr.curve(x[ind.pos], x[ind.neg])$auc.integral })
  rocs<-sapply(data.test[,-ncol(data.test)], function(x){ roc.curve(x[ind.pos], x[ind.neg])$auc })
  
  return(rbind(prs, rocs))
  
}

res<-vector("list",ncol(combs))

for( r in 1:ncol(combs) ){
  
  res[[r]]<-fit_ridge(r)
  print(r)
  
}

# PR and ROC boxplots
perfs_summary<-do.call(rbind, lapply(res, function(x)x[[1]]) )
perfs_summary<-data.frame(perfs_summary)
perfs_summary$performance.measure<-c("PR","ROC")
perfs_summary<-melt(perfs_summary)
perfs_summary$feature.set<-str_extract("chromatin|sequence|seq.chr.comb",string = perfs_summary$variable)
is.comb<-str_extract(pattern="comb", as.character(perfs_summary$variable))
is.comb[ is.na(is.comb) ]<-as.character(perfs_summary$feature.set)[ is.na(is.comb) ]
perfs_summary$is.comb<-factor(is.comb)
perfs_summary$feature.set<-factor(perfs_summary$feature.set, levels=c("sequence","chromatin","seq.chr.comb"))
perfs_summary$variable<-factor(perfs_summary$variable, levels=levels(perfs_summary$variable)[c(1:11,13,12)])
ggplot(perfs_summary[ perfs_summary$performance.measure=="PR", ], aes(x=variable, y=value, fill=feature.set:is.comb))+theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1), panel.border=element_rect(fill=NA,linetype = 1,colour = "black"), legend.position="bottom")+geom_boxplot(outlier.colour = "white")+geom_jitter(size=0.5, width=0.5, alpha=0.7)+facet_grid(.~performance.measure+feature.set, space = "free_x", scales="free_x")+scale_fill_manual(values=brewer.pal(6,"Paired")[c(2:1,3:4,6)])+ylab("AUPRC")+xlab("Method")
ggplot(perfs_summary[ perfs_summary$performance.measure=="ROC", ], aes(x=variable, y=value, fill=feature.set:is.comb))+theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1), panel.border=element_rect(fill=NA,linetype = 1,colour = "black"), legend.position="bottom")+geom_boxplot(outlier.colour = "white")+geom_jitter(size=0.5, width=0.5, alpha=0.7)+facet_grid(.~performance.measure+feature.set, space = "free_x", scales="free_x")+scale_fill_manual(values=brewer.pal(6,"Paired")[c(2:1,3:4,6)])+ylab("AUROC")+xlab("Method")

# finally, fitting a single ridge model for genome-wide predictions:
preds<-cbind( "chromatin.LASSO"=preds.chromatin$lasso$all, "chromatin.SVMlin"=preds.chromatin$svm.lin$all, "chromatin.SVMrad"=preds.chromatin$svm.rad$all, "chromatin.RF"=preds.chromatin$rf$all, "sequence.LASSO"=preds.seq$lasso$all, "sequence.SVMlin"=preds.seq$svm.lin$all, "sequence.RF"=preds.seq$rf$all)
y<-unlist( lapply(predictions.chromatin$lasso$all, function(x){x$test.pred[,2]}))

weights<-rep(1,length(y))
tab<-table(y)
weights[ y=="positive" ]<-as.numeric(tab["negative"]/tab["positive"])

ridge<-cv.glmnet( as.matrix(preds), y = y, weights = weights, type.measure = "mse", family="binomial", alpha=0 )
coef(ridge)

# done!
