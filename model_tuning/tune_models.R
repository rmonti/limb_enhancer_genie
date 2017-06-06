

# if you haven't done so already, install  the tunetest package from source (located inside this folder)
# install.packages("./tunetest", repos=NULL, type="source")

require(tunetest)
data(rf)
data(lasso)
data(svm.lin)
data(svm.rad)

# load required packages
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
require(doMC)

# adjust the number of cores used to tune the models:
registerDoMC(cores = 32)

# loading the data
# chromatin:
load("ML_data_summary_train_scaled.Rdata")
# -> ML_data_summary.train
# these are the Dnase-I, Dnase-I Enrichment, Histone modifications, p300, cohesin, and methylation features 

# sequence:
load("clustp_scores_train_scaled.Rdata")
# -> clustp_scores.train
# these are the features calculated by identifying conserved TF-binding motif-clusters within the elements

# VISTA elements and associated metadata:
load("velements.Rdata")
# -> velements
# velements is a GRanges object that contains the metadata associated with the rows of ML_data_summary.train and clustp_scores.train

# sanity checks
# stopifnot(all(row.names(ML_data_summary.train) == row.names(clustp_scores.train) ))
# stopifnot(all(names(velements) == row.names(ML_data_summary.train) ))

# labels for the data ( positive / negative -> active / inactive in limb)
y<-factor(velements$limb, labels=c("negative","positive"))

# creating folds for model tuning/testing:
set.seed(1991)
folds<-createFolds(y = y, k = 10, list = TRUE, returnTrain = TRUE )

# setting cross-validataion scheme and summary function (see caret package for details)
# defines the scheme to use for tuning parameters (e.g. lambda for LASSO, or cost/gamma for radial SVM) 
fitControl<-trainControl(method="cv", number=10, classProbs = TRUE, summaryFunction = twoClassSummary.2, returnData = FALSE)

# defining the models that should be fit on the chromatin / sequence features
models.chromatin<-c("lasso","svm.lin","svm.rad","rf")
models.sequence<-c("lasso","svm.lin","rf")

# defining the feature sets (columns) we want to use for the chromatin / sequence models
# we select all features:
# we could also limit ourselves to subsets of features or single features, just add them to the list!
features.chromatin<-list('all'=1:ncol(ML_data_summary.train))
features.seq<-list("all"=1:ncol(clustp_scores.train))

# we construct the fit.info objects used to train the models from the objects defined above,
# this will generate nested lists defining which data, observations, models and features to use
fit.infos.chrom<-lapply(models.chromatin, function(model){ lapply(features.chromatin, function(features){ lapply(folds, function(fold){ list(model=model, data="ML_data_summary.train", y=y, features=features, observations=fold, fitControl=fitControl)}) })})
names(fit.infos.chrom)<-models.chromatin

fit.infos.seq<-lapply(models.sequence, function(model){ lapply(features.seq, function(features){ lapply(folds, function(fold){ list(model=model, data="clustp_scores.train", y=y, features=features, observations=fold, fitControl=fitControl)}) })})
names(fit.infos.seq)<-models.sequence

# the tuning grids for the SVM sequence-models are slightly adjusted:
# a tuning grid this fine will need a lot of computational power so I would suggest using rougher one for testing this "at home" :P
grid.svm.seq<-expand.grid(cost = 10^seq(0, -7, length=60), gamma = 1, weight = 1)
for ( i in 1:length(fit.infos.seq$svm.lin$all) ) fit.infos.seq$svm.lin$all[[i]]$grid <- grid.svm.seq

tunedmodels.chrom<-vector("list", length(fit.infos.chrom))
tunedmodels.chrom<-lapply(tunedmodels.chrom, function(x)vector("list",length(features.chromatin)))

# the two nested loops below fit the models.
# most of the code is hidden within the call to "fit_model"
# see the fit_model() function
# it will take care of balancing the class weights etc, and wraps caret::train() using custom settings
fit_model()

# this should run on any computer in a reasonable amount of time:
for ( i in 1:length(tunedmodels.chrom) ){
  
  for ( j in 1:length(tunedmodels.chrom[[i]]) ){
    
      t<-proc.time()
      tunedmodels.chrom[[i]][[j]]<-foreach( m = 1:length(fit.infos.chrom[[i]][[j]]) ) %dopar% fit_model(fit.infos.chrom[[i]][[j]][[m]])
      names(tunedmodels.chrom[[i]])[j]<-names(fit.infos.chrom[[i]])[j]
      print(proc.time()-t)
      print( paste("chromatin models fit for: model =", names(fit.infos.chrom)[i],", features =", names(fit.infos.chrom[[i]])[j]) )
  }
  
  names(tunedmodels.chrom)[i]<-names(fit.infos.chrom)[i]

}

save(tunedmodels.chrom, file = "chrom.models.Rdata")

tunedmodels.seq<-vector("list", length(fit.infos.seq))
tunedmodels.seq<-lapply(tunedmodels.seq, function(x)vector("list",length(features.seq)))

# this can take very long:
for ( i in 1:length(tunedmodels.seq) ){
  
  for ( j in 1:length(tunedmodels.seq[[i]]) ){
    
    t<-proc.time()
    tunedmodels.seq[[i]][[j]]<-foreach( m = 1:length(fit.infos.seq[[i]][[j]]) ) %dopar% fit_model(fit.infos.seq[[i]][[j]][[m]])
    names(tunedmodels.seq[[i]])[j]<-names(fit.infos.seq[[i]])[j]
    print(proc.time()-t)
    print( paste("sequence models fit for: model =", names(fit.infos.seq)[i],", features =", names(fit.infos.seq[[i]])[j]) )
    
  }
  
  names(tunedmodels.seq)[i]<-names(fit.infos.seq)[i]
  
}

save(tunedmodels.seq, file = "seq.models.Rdata")

# q(save = "no")
