coef.lasso.fit <-
function(object, s="lambda.min", ...){
  
  if ( s == "lambda.min" ){
    
    return(coef(object$fit$finalModel, s=object$fit$bestTune$lambda))
    
  }else if( s == "lambda.1se" ){
    
    res<-object$fit$results
    lambda.min<-object$fit$bestTune$lambda
    
    res<-res[ order(res$lambda, decreasing=TRUE), ]
    lambda.1se<-res$lambda[ which( res$ROC + res$ROCSD >= res$ROC[ res$lambda == lambda.min ][1] )[1] ]
    
    return(coef(object$fit$finalModel, s=lambda.1se))
    
  }else  return(coef(object$fit$finalModel, s=s))
  
}
