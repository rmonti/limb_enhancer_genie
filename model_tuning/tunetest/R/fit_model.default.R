fit_model.default <-
function(fit.info=NULL, ...){
  
  if( length(fit.info)==0 ){ cat("need fit.info object (list) with following named elements: \n$model:\tone of either 'svm.lin', 'svm.rad', 'lasso', or 'rf' \n$data:\tname of the data X to use (not the actual data object) \n$y:\tthe classes (factor variable, labels have to be 'positive' and 'negative') \n$features:\tthe feature subset to use \n$observations:\tthe observations to use for training (all observations NOT in this set will be used for testing when 'predict(x.fit)' is called)\n$fitControl:\tan instance of trainControl defining the CV-scheme \n$grid:\toptional tuning grid")
                             stop("please define fit.info") }
  
  if ( is.character(fit.info$features) ){ stopifnot( all(fit.info$features %in% colnames(eval(parse(text=fit.info$data)))))
  }else(stopifnot(all(fit.info$features %in% 1:ncol(eval(parse(text=fit.info$data))))))
                                          
  stopifnot( fit.info$model %in% c("rf","svm.lin","svm.rad","lasso") )
  
  model<-vector( "list", length = 2 )
  names(model)<-c("fit","args" )
  model$args<-fit.info
  
  if (length(fit.info$features)==1){ model$args$single.feature<-TRUE
  }else model$args$single.feature<-FALSE
  
  if ( fit.info$model == "rf" ){
    
    if( length(fit.info$grid)==0){ grid<-expand.grid( mtry=sqrt(length(fit.info$features)), ntree=1000, strat=1 ) 
                                   }else  grid<-fit.info$grid
    
    model$fit<-train( x = eval(parse( text=fit.info$data ))[ fit.info$observations, fit.info$features, drop=FALSE ], y = fit.info$y[ fit.info$observations ], method=eval(parse(text=fit.info$model)), ... , importance=TRUE, tuneGrid = grid, trControl = fit.info$fitControl, metric = "ROC")
    
  }else if ( fit.info$model == "svm.lin" ){
    
    # balancing class weights:
    n.obs<-table( fit.info$y[ fit.info$observations ] )
    wgt<-as.numeric( n.obs["negative"] / n.obs[ "positive" ] )
    
    # ordering such that positive overvations come first, this prevents the decision values from getting messed up
    x = eval(parse( text=fit.info$data ))[ fit.info$observations, fit.info$features, drop=FALSE ]
    y = fit.info$y[ fit.info$observations ]
    
    x<-x[ order(y, decreasing = TRUE),, drop=FALSE ]
    y<-y[ order(y, decreasing = TRUE) ]
    
    if( length(fit.info$grid)==0 ){ grid <- expand.grid( cost=10^seq(1,-5,by=-0.2),gamma=1, weight=wgt )
                                    }else grid <- fit.info$grid
    grid$weight<-wgt
    
    model$fit<-train( x = x, y = y, method=eval(parse(text=fit.info$model)), ... , tuneGrid = grid, trControl = fit.info$fitControl, metric = "ROC")
    
  }else if ( fit.info$model == "svm.rad" ){
    
    # balancing class weights:
    n.obs<-table( fit.info$y[ fit.info$observations ] )
    wgt<-as.numeric( n.obs["negative"] / n.obs[ "positive" ] )
    
    # ordering such that positive overvations come first, this prevents the decision values from getting messed up
    x = eval(parse( text=fit.info$data ))[ fit.info$observations, fit.info$features, drop=FALSE ]
    y = fit.info$y[ fit.info$observations ]
    
    x<-x[ order(y, decreasing = TRUE),, drop=FALSE ]
    y<-y[ order(y, decreasing = TRUE) ]
    
    if( length(fit.info$grid)==0){ #grid<-expand.grid( cost=10^seq(1,-4,by=-0.5), gamma=seq(1,0.0001,length=10), weight=wgt )
                                    grid<-expand.grid( cost=10^seq(2,-3,by=0.5), gamma=10^seq(-1,-7,length=10), weight=wgt )
                                  }else grid<-fit.info$grid
    grid$weight<-wgt
    
    model$fit<-train( x = x, y = y, method=eval(parse(text=fit.info$model)), ... , tuneGrid = grid, trControl = fit.info$fitControl, metric = "ROC")
    
    
  }else if ( fit.info$model == "lasso" ){
    
    # balancing class weights:
    n.obs<-table( fit.info$y[ fit.info$observations ] )
    wgt<-as.numeric( n.obs["negative"] / n.obs[ "positive" ] )
    
    # ordering such that positive overvations come first, this prevents the decision values from getting messed up
    weights<-rep(1, length(fit.info$observations) )
    weights[ fit.info$y[ fit.info$observations ] == "positive" ] <- wgt
    
    if( length(fit.info$grid)==0 ){ grid<-expand.grid(alpha=1, lambda=10^seq(0,-6,by=-0.1))
                                    }else grid<-fit.info$grid
    
    if( length(fit.info$features)==1 ){ x <- eval(parse( text=fit.info$data ))[ fit.info$observations, fit.info$features, drop=FALSE ]
                                        x <- as.matrix( cbind(x, dummy=rep(1,nrow(x))) )
                                        
    } else x<-eval(parse( text=fit.info$data ))[ fit.info$observations, fit.info$features, drop=FALSE ]
    
    model$fit<-train( x = x, y = fit.info$y[ fit.info$observations ], method=eval(parse(text=fit.info$model)), weights = weights, ... , tuneGrid = grid, trControl = fit.info$fitControl, metric = "ROC")
    
  }
  
  class(model)<-paste0( fit.info$model, ".fit" )
  return(model)
  
}
