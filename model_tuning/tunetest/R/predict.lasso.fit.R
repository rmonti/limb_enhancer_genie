predict.lasso.fit <-
function(x, plot=FALSE, s="lambda.min",... ){
  
  require(PRROC)
  
  if ( s == "lambda.min" ){
    
    lambda<-x$fit$bestTune$lambda
    
  }else if( s == "lambda.1se" ){
    
    res<-x$fit$results
    lambda.min<-x$fit$bestTune$lambda
    
    res<-res[ order(res$lambda, decreasing=TRUE), ]
    lambda<-res$lambda[ which( res$ROC + res$ROCSD >= res$ROC[ res$lambda == lambda.min ][1] )[1] ]
    
  }else  lambda<-s
  
  x.test<-as.matrix( eval(parse(text=x$args$data))[ -1*(x$args$observations), x$args$features, drop=FALSE ])
  
  if ( x$args$single.feature ) x.test <- cbind(x.test, dummy=rep(1,nrow(x.test)) )
  
  test.predictions <- predict( x$fit$finalModel, newx=x.test , s=lambda, type="response" )
  test.y<-x$args$y[ -1*(x$args$observations) ]
  
  pr<-pr.curve( scores.class0 = test.predictions[ test.y == "positive" ], scores.class1 = test.predictions[ test.y == "negative" ], curve = TRUE, rand.compute = TRUE)
  roc<-roc.curve( scores.class0 = test.predictions[ test.y == "positive" ], scores.class1 = test.predictions[ test.y == "negative" ], curve = TRUE, rand.compute = TRUE)
  
  preds<-list("test.pred"= data.frame("predictions"=test.predictions, "y"=test.y), "pr.curve"=pr, "roc.curve"=roc )
  class(preds)<-"test.prediction"
  
  return(preds)
  
}
