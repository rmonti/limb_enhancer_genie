plot.test.prediction <-
function(x, type="roc", ... ){
  
  if ( type=="roc" ) plot(x$roc.curve, ... )
  if ( type=="pr"  ) plot(x$pr.curve, ... )
  
}
