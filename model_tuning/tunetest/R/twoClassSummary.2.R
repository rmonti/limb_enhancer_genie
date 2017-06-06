twoClassSummary.2 <-
function (data, lev = NULL, model = NULL){
    if (length(levels(data$obs)) > 2) 
      stop(paste("Your outcome has", length(levels(data$obs)), 
                 "levels. The twoClassSummary() function isn't appropriate."))
    requireNamespace("pROC")
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
      stop("levels of observed and predicted data do not match")
    rocObject <- try(pROC::roc(data$obs, data[, lev[1]]), silent = TRUE)
    rocAUC <- if (class(rocObject)[1] == "try-error") 
      NA
    else rocObject$auc
    
    pos<-data[,"positive"][data[,"obs"]=="positive"]
    neg<-data[,"positive"][data[,"obs"]=="negative"]
    
    pr<-as.vector(pr.curve(scores.class0 = pos, scores.class1 = neg)$auc.integral)
    
    out <- c(rocAUC, pr, sensitivity(data[, "pred"], data[, "obs"], 
                                     lev[1]), specificity(data[, "pred"], data[, "obs"], lev[2]))
    names(out) <- c("ROC","AUPRC", "Sens", "Spec")
    out
}
