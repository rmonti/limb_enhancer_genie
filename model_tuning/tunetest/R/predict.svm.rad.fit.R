predict.svm.rad.fit <-
function(x, plot=FALSE, type="dec.values", ... ){
  
  stopifnot(type %in% c("dec.values","prob") )
  
  require(PRROC)
  if( type=="dec.values") test.predictions<-unlist( attributes( predict( x$fit$finalModel, newdata=eval(parse(text=x$args$data))[ -1*(x$args$observations), x$args$features ], decision.values=TRUE ))$decision.values )
  if (type=="prob") test.predictions<-unlist( attributes( predict( x$fit$finalModel, newdata=eval(parse(text=x$args$data))[ -1*(x$args$observations), x$args$features ], probability=TRUE ))$probabilities[,"positive"] )
  test.y<-x$args$y[ -1*(x$args$observations) ]
  pr<-pr.curve( scores.class0 = test.predictions[ test.y == "positive"], scores.class1 = test.predictions[ test.y == "negative" ], curve = TRUE, rand.compute = TRUE)
  roc<-roc.curve( scores.class0 = test.predictions[ test.y == "positive"], scores.class1 = test.predictions[ test.y == "negative" ], curve = TRUE, rand.compute = TRUE)
  
  preds<-list("test.pred"= data.frame("predictions"=test.predictions, "y"=test.y), "pr.curve"=pr, "roc.curve"=roc )
  class(preds)<-"test.prediction"
  
  return(preds)
  
}
