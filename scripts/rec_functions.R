

k_validate <- function(seed){
  set.seed(seed)
  
  
  out <- list()
  
  out[[1]] <- getROC(train.seq, train.outcomes, test.seq, test.outcomes)
  
  
  avg_AUCPR <- lapply(out, function(x){ x[[1]]})
  message("\nAccuracy average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}


# Function create and test RF model #
getROC <- function(train.seq, train.outcomes, test.seq, test.outcomes){
  if(all(train.seq$Sample == train.outcomes$Sample) == FALSE){
    stop("Training Sample_IDs do not match")
  }
  if(all(test.seq$Sample == test.outcomes$Sample) == FALSE){
    stop("Testing Sample_IDs do not match")
  }
  model.training <- randomForest(x = train.seq[,-1, drop=FALSE],
                                 y = as.factor(train.outcomes$Reason.for.discontinuation))
  
  test.preds <- predict(model.training, test.seq[,-1, drop=FALSE])
  
  # print(model.training)
  # Check variable importance
  #varImpPlot(model.training)
  # Prediction confusion matrix
  print(table(observed = test.outcomes$Reason.for.discontinuation,
        predicted = test.preds))
  
  prediction_for_roc_curve <- predict(model.training, 
                                      test.seq[,-1, drop=FALSE],
                                      type="prob")
  
  pred <- prediction(prediction_for_roc_curve[,2], 
    test.outcomes$Reason.for.discontinuation, 
    label.ordering = 
    c("Therapy complete", "Progression"))
  #for making aucpr curves
  # scores <- data.frame(prediction_for_roc_curve[,2], test.outcomes$Reason.for.discontinuation)
  # 
  # scores <- scores %>%
  #   mutate(score = case_when((test.outcomes.Reason.for.discontinuation == "Therapy complete") ~ 1,
  #                            TRUE ~ 0))
  # print(prediction_for_roc_curve[,2])
  # aucpr <- pr.curve(scores.class0=scores[scores$score=="0",]$`prediction_for_roc_curve...2.`,
  #                   scores.class1=scores[scores$score=="1",]$`prediction_for_roc_curve...2.`,
  #                   curve=T)
  # 
  # y <- as.data.frame(aucpr$curve)
  # print(ggplot(y, aes(V1, V2))+geom_path()+ylim(0,1))
  
  perf <- performance(pred, "tpr", "fpr")
  
  AUCPR_ROCR <- performance(pred, measure = "aucpr")
  print(AUCPR_ROCR)
  AUCPR <- AUCPR_ROCR@y.values[[1]]
  print(AUCPR)
  
  df <- data.frame(FalsePositive=c(perf@x.values[[1]]),
                   TruePositive=c(perf@y.values[[1]]))
  out <- list(AUCPR, df)
  
  return(out)
}


# Grabbing AUROC values from input models #
grabVals <- function(input){
  out <- list()
  for(i in 1:length(seed_list)){
    out <- append(out, input[[i]][[1]][[1]])
  }
  out
}


# Main function call to generate RF models using 25 seeds #
kTest <- function(seed_list){
  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i])
  }
  out <- grabVals(out)
  
  # CW edit to function
  out <- do.call(rbind.data.frame, out)
  colnames(out) <- c("AUCPR")
  out
}

