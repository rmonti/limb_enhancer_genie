plot_predictionlist<-function(predlist, type="ROC", return.data=FALSE, r=0.05, fun.min="min", fun.max="max", col="black", ribbon=FALSE ){
  
  stopifnot( type %in% c("ROC","AUPRC") )
  
  if ( type=="ROC" ){
    
    curves<-lapply(predlist, function(x){data.frame(x$roc$curve)})    
    curves<-lapply(curves, setNames, c("FPR","TPR","Cutoff") )
    
    n <- sapply(curves, nrow)
    
    var.x<-"FPR"
    var.y<-"TPR"
    
    curves<-data.frame( do.call(rbind, curves) )
    
  }
  
  if ( type=="AUPRC" ){
    
    curves<-lapply(predlist, function(x){data.frame(x$pr$curve)})
    curves<-lapply(curves, setNames, c("Recall","Precision","Cutoff") )
    
    n <- sapply(curves, nrow)
    
    var.x<-"Recall"
    var.y<-"Precision"
    
    curves<-data.frame( do.call(rbind, curves) )
    
  }
  
  # curves[,var.x]<-floor(curves[,var.x]/r)*r
  # curves[,var.x]<-round(curves[,var.x]/r)*r
  
  curves[,var.x]<-cut(curves[,var.x], seq(0,1,by=r), include.lowest=TRUE, labels = rollmean(seq(0,1,by=r),k=2), ordered_result = TRUE ) 
  curves[,var.x]<-as.numeric(as.character(curves[,var.x]))
  
  if ( type == "AUPRC" ){
    
    curves.summary<-data.table(curves)[, list("mean"=mean(Precision),"sd"=sd(Precision),"min"=min(Precision),"max"=max(Precision),"mean.cutoff"=mean(Cutoff),"median.cutoff"=median(Cutoff),"sd.cutoff"=sd(Cutoff),"max.cutoff"=max(Cutoff),"min.cutoff"=min(Cutoff)), by=Recall]
    setnames(curves.summary, "Recall","x")
    
  }
  
  
  if ( type == "ROC" ){
    
    curves.summary<-data.table(curves)[, list("mean"=mean(TPR),"sd"=sd(TPR),"min"=min(TPR),"max"=max(TPR),"mean.cutoff"=mean(Cutoff),"median.cutoff"=median(Cutoff),"sd.cutoff"=sd(Cutoff),"max.cutoff"=max(Cutoff),"min.cutoff"=min(Cutoff)), by=FPR]
    setnames(curves.summary, "FPR", "x")
    
  }
  
  if (ribbon){ print( ggplot()+
          geom_ribbon(data=data.frame(curves.summary),aes(ymin=curves.summary$min, ymax=curves.summary$max, x=x), alpha=0.075, fill=col, col=col)+
          geom_line(data=data.frame(curves.summary), aes(x=x, y=mean), alpha=0.75, linetype=3)+
          stat_summary(data=curves, aes_string(x=var.x, y=var.y), fun.y=mean, fun.ymax = fun.max, fun.ymin=fun.min)+
          xlab(var.x)+
          ylab(names(curves)[2])+
          theme_minimal())

  } else{print( ggplot()+
          geom_line(data=data.frame(curves.summary), aes(x=x, y=mean), alpha=0.75, linetype=3)+
          stat_summary(data=curves, aes_string(x=var.x, y=var.y), fun.y=mean, fun.ymax = fun.max, fun.ymin=fun.min)+
          geom_point(data=data.frame(curves.summary), aes(col=mean.cutoff, x=x, y=mean))+
          xlab(var.x)+
          ylab(names(curves)[2])+
          theme_minimal()+
          scale_colour_gradientn(colours=rev(brewer.pal(9,"Spectral"))))
  
          p<-list( geom_line(data=data.frame(curves.summary), aes(x=x, y=mean), alpha=0.75, linetype=3),
            stat_summary(data=curves, aes_string(x=var.x, y=var.y), fun.y=mean, fun.ymax = fun.max, fun.ymin=fun.min),
            geom_point(data=data.frame(curves.summary), aes(col=mean.cutoff, x=x, y=mean)),
            xlab(var.x),
            ylab(names(curves)[2]),
            theme_minimal(),
            scale_colour_gradientn(colours=rev(brewer.pal(9,"Spectral"))) )
  
  }
  
  
  curves.summary<-data.frame(curves.summary)
  colnames(curves.summary)<-c(var.x, paste(names(curves)[2],colnames(curves.summary)[2:5],sep="."),"mean.cutoff","median.cutoff","sd.cutoff","max.cutoff","min.cutoff")
  
  if ( return.data==TRUE ){ return(curves.summary)
  }else invisible(p)
  
}

.datatable.aware=TRUE



