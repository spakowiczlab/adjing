
sorting <- function(data, cohort, outcome1, outcome2){
  filter(data, 
         Cohort == cohort, 
         Reason.for.discontinuation == outcome1 | 
           Reason.for.discontinuation == outcome2)  
}



k_validate <- function(seed, neg.outcome, pos.outcome, modelversion, train.seq, test.seq, train.outcomes, test.outcomes){
  set.seed(seed)
  

  
  out <- list()
  
  out[[1]] <- getROC(seed, train.seq, train.outcomes, test.seq, test.outcomes, neg.outcome, pos.outcome, modelversion)
  
  
  avg_AUCPR <- lapply(out, function(x){ x[[1]]})
  message("\nAccuracy average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}


# Function create and test RF model #
getROC <- function(seed, train.seq, train.outcomes, test.seq, test.outcomes, neg.outcome, pos.outcome, modelversion){
  if (modelversion == "Boruta") {
    ###################### Boruta for selection to the model ######################
    set.seed(seed)
    
    # Drop Patient Id for modeling
    x_train <- train.seq[ , -1, drop = FALSE ]
    y_train <- as.factor(train.outcomes[["Reason.for.discontinuation"]])
    
    # Combine for Boruta formula interface
    boruta_data <- cbind(y = y_train, x_train)
    
    # Run Boruta
    bor <- Boruta(y ~ ., data = boruta_data, doTrace = 2, maxRuns = 2500)
    
    # Resolve tentative features
    bor_fixed <- TentativeRoughFix(bor)
    
    # Get confirmed features only
    selected_vars <- getSelectedAttributes(bor_fixed, withTentative = FALSE)
    
    min_vars <- 1   # set your desired minimum number of predictors
    
    if(length(selected_vars) < min_vars){
      warning("Boruta confirmed fewer than min_vars — including tentative features")
      selected_vars <- getSelectedAttributes(bor_fixed, withTentative = TRUE)
    }
    
    # rank based RF selection
    if(length(selected_vars) < min_vars){
      warning("Boruta still selected fewer than min_vars — filling with top RF importance features")
      
      # Fit quick RF on all predictors
      rf_temp <- randomForest(x = x_train, y = y_train)#, ntree = max(500, tree))
      
      # Rank features by importance
      imp <- importance(rf_temp)
      top_vars <- rownames(imp)[order(-imp[,1])]
      
      # Add top features until min_vars is met
      selected_vars <- unique(c(selected_vars, top_vars[1:min_vars]))
    }
    
    message("\nBoruta selected ", length(selected_vars), " variables.")
    
    print(selected_vars)
    
    # Restrict train/test to selected features
    train.seq <- train.seq[ , c("Sample", selected_vars), drop = FALSE ]
    test.seq  <- test.seq [ , c("Sample", selected_vars), drop = FALSE ]
    ###############################################################################
  } 
  
  if(all(train.seq$Sample == train.outcomes$Sample) == FALSE){
    stop("Training Sample_IDs do not match")
  }
  if(all(test.seq$Sample == test.outcomes$Sample) == FALSE){
    stop("Testing Sample_IDs do not match")
  }
  set.seed(seed)
  model.training <- randomForest(x = train.seq[,-1, drop=FALSE], y = 
                          as.factor(train.outcomes$Reason.for.discontinuation), importance=TRUE)
  
  test.preds <- predict(model.training, 
                        test.seq[,-1, drop=FALSE])

  
  # Prediction confusion matrix
  pred_cm <- table(observed = test.outcomes$Reason.for.discontinuation,
                   predicted = test.preds)
  
  print(pred_cm)
  
  
  prediction_for_roc_curve <- predict(model.training, 
                                      test.seq[,-1, drop=FALSE],
                                      type="prob")
  
  pred <- prediction(prediction_for_roc_curve[,2], 
                     test.outcomes$Reason.for.discontinuation, 
                     label.ordering = 
                       c(neg.outcome, pos.outcome))

  perf <- performance(pred, "tpr", "fpr")
  
  AUCPR_ROCR <- performance(pred, measure = "aucpr")
  print(AUCPR_ROCR)
  AUCPR <- AUCPR_ROCR@y.values[[1]]
  print(AUCPR)
  
  df <- data.frame(FalsePositive=c(perf@x.values[[1]]),
                   TruePositive=c(perf@y.values[[1]]))
  
  if (modelversion == "Boruta") {
    out <- list(AUCPR, df, pred_cm, selected_vars)
  }
  else
  {
    out <- list(AUCPR, df, pred_cm)
  }
  return(out)
}


# Grabbing AUROC values from input models #
grabVals <- function(output, input, seed_list){
  output <- list()
  
  for(i in 1:length(seed_list)){
    output <- append(output, input[[i]][[1]][[1]])
  }
  
  output <- do.call(rbind.data.frame, output)
  colnames(output) <- c("AUCPR")
  output
}


# Main function call to generate RF models 
kTest <- function(seed_list, neg.outcome, pos.outcome, modelversion, outcomeobject, dataobject, traincohort, testcohort, genelist){
  
  if (deparse(substitute(genelist)) == "clin_list") {
    
    dataobject <- dataobject %>%
    mutate(
      across(
        .cols = Stage:Sex,
        .fns = as.factor
      )
    )
  }
    train.outcomes <- sorting(outcomeobject,
                            traincohort,
                            neg.outcome,
                            pos.outcome)
  
  train.outcomes <- arrange(train.outcomes, desc(Sample))
  
  
  test.outcomes <- sorting(outcomeobject,
                           testcohort,
                           neg.outcome,
                           pos.outcome)
  
  test.outcomes <- arrange(test.outcomes, desc(Sample))
  
  train.outcomes <-
    train.outcomes %>%
    select(Sample,
           `Reason.for.discontinuation`)
  
  
  test.outcomes <-
    test.outcomes %>%
    select(Sample,
           `Reason.for.discontinuation`)
  
  if (modelversion == "Random") {
    train.outcomes$Reason.for.discontinuation<-
      (sample(train.outcomes$Reason.for.discontinuation))
    
    test.outcomes$Reason.for.discontinuation<-
      (sample(test.outcomes$Reason.for.discontinuation))
    
  }
  

  train.seq <- filter(dataobject,
                      (Sample %in% train.outcomes$Sample))
  train.seq <- arrange(train.seq, desc(Sample))
  
  
  test.seq <- filter(dataobject,
                     (Sample %in% test.outcomes$Sample))
  test.seq <- arrange(test.seq, desc(Sample))
  
  
  
  train.seq <- train.seq %>% select(Sample, 
                                    genelist
  )
  
  test.seq <- test.seq %>% select(Sample, 
                                  genelist
  )
  
    

  

  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i], neg.outcome, pos.outcome, modelversion, train.seq, test.seq, train.outcomes, test.outcomes)
  }
  out
}




p.calc <- function(data, random){
  avg <- (data)
  p <- ecdf(random)
  p(avg)
  
}




add.metrics <- function(randomdata, output, pos.outcome){
  test_conmat2 <- data.frame(run = 1:25000,
                             Accuracy = NA,
                             Specificity = NA,
                             Precision = NA,
                             Recall = NA,
                             F1 = NA)
  
  for (i in 1:25000) {
    test_conmat <- confusionMatrix(randomdata[[i]][[1]][[3]],
                                   positive = pos.outcome)
    
    test_conmat2$Accuracy[i] <- test_conmat[[3]][1]
    test_conmat2$Specificity[i] <- test_conmat[[4]][2]
    test_conmat2$Precision[i] <- test_conmat[[4]][5]
    test_conmat2$Recall[i] <- test_conmat[[4]][6]
    test_conmat2$F1[i] <- test_conmat[[4]][7]
    
    
  }
  test_conmat2[is.na(test_conmat2)] <- 0
  x = 1
  output <- data.frame(run = 1:1000,
                       Avg.Accuracy = NA,
                       Avg.Specificity = NA,
                       Avg.Precision = NA,
                       Avg.Recall = NA,
                       Avg.F1 = NA)
  
  for (i in 1:1000) {
    y = x+24
    
    output$Avg.Accuracy[i] <- mean(test_conmat2$Accuracy[x:y])
    output$Avg.Specificity[i] <- mean(test_conmat2$Specificity[x:y])
    output$Avg.Precision[i] <- mean(test_conmat2$Precision[x:y])
    output$Avg.Recall[i] <- mean(test_conmat2$Recall[x:y])
    output$Avg.F1[i] <- mean(test_conmat2$F1[x:y])
    
    x = x+25
  }
  output
}


get_vars <- function(myList) {
  unique(unlist(lapply(myList, function(x) rownames(x[[1]][[4]]))))
}

grabImp <- function(input, seed_list){
  
  out <- list()
  
  for (i in seq_along(seed_list)) {
    
    # Extract getROC() output for this seed
    res <- input[[i]][[1]]
    
    impdf <- res[[4]] %>%   # Genes selected by boruta
      as.data.frame() 
    
    impdf$seed  <- seed_list[i]
    impdf$run_number  <- i
    
    out[[i]] <- impdf
  }
  
  # bind into one big data frame
  dplyr::bind_rows(out)
}



#running random model
shuffled_rf <- function(seed_list, neg.outcome, pos.outcome, outcomeobject, dataobject, traincohort, testcohort, genelist)
{
output <- list()

for (j in 1:length(run_number)) {
  set.seed(seed_list_r[j])
  

  
  shuffled_results <- kTest(seed_list,
                            neg.outcome,
                            pos.outcome,
                            "Random",
                            outcomeobject,
                            dataobject,
                            traincohort,
                            testcohort,
                            genelist) 
  
  output <- append(output, shuffled_results)

  output
}
output
}



met.p.calc <- function(randomdata, pos.outcome, dataobject){
  #function that generates averages of metrics for randomized run
  random_met <- add.metrics(randomdata,
                            random_met,
                            pos.outcome)
  
  
  #storing experimental run additional metrics   
  temp_conmat2 <- data.frame(run = 1:25,
                             Accuracy = NA,
                             Specificity = NA,
                             Precision = NA,
                             Recall = NA,
                             F1 = NA)
  
  #generating additional metrics for experimental run    
  for (i in 1:25) {
    temp_conmat <- confusionMatrix(dataobject[[i]][[1]][[3]],
                                   positive = pos.outcome)
    
    temp_conmat2$Accuracy[i] <- temp_conmat[[3]][1]
    
    temp_conmat2$Specificity[i] <- temp_conmat[[4]][2]
    
    temp_conmat2$Precision[i] <- temp_conmat[[4]][5]
    
    temp_conmat2$Recall[i] <- temp_conmat[[4]][6]
    
    temp_conmat2$F1[i] <- temp_conmat[[4]][7]
    
  }
  #storing percentiles    
  output <- data.frame(Accuracy = NA,
                       Specificity = NA,
                       Precision = NA,
                       Recall = NA,
                       F1 = NA)
  
  
  avg.acc <- mean(temp_conmat2$Accuracy)
  output$Accuracy[1] <- p.calc(avg.acc,
                               random_met$Avg.Accuracy)
  
  avg.s <- mean(temp_conmat2$Specificity)
  output$Specificity[1] <- p.calc(avg.s,
                                  random_met$Avg.Specificity)
  
  avg.prec <- mean(temp_conmat2$Precision)
  output$Precision[1] <- p.calc(avg.prec,
                                random_met$Avg.Precision)
  
  avg.recall <- mean(temp_conmat2$Recall)
  output$Recall[1] <- p.calc(avg.recall,
                             random_met$Avg.Recall)
  
  avg.F1 <- mean(temp_conmat2$F1)
  avg.F1[is.na(avg.F1)] <- 0
  
  output$F1[1] <- p.calc(avg.F1,
                         random_met$Avg.F1)
  
  
  print(output$F1)
  
  
  output
}




metric_graph <- function(randomdata, pos.outcome, dataobject){
  
  metricsdf <- add.metrics(randomdata,
                           metricsdf,
                           pos.outcome)
  
  modelmetricsdf <- data.frame(run = 1:25,
                               Accuracy = NA,
                               Specificity = NA,
                               Precision = NA,
                               Recall = NA,
                               F1 = NA)
  
  
  for (i in 1:25) {
    temp_conmat <- confusionMatrix(dataobject[[i]][[1]][[3]],
                                   positive = pos.outcome)
    
    modelmetricsdf$Accuracy[i] <- temp_conmat[[3]][1]
    
    modelmetricsdf$Specificity[i] <- temp_conmat[[4]][2]
    
    modelmetricsdf$Precision[i] <- temp_conmat[[4]][5]
    
    modelmetricsdf$Recall[i] <- temp_conmat[[4]][6]
    
    modelmetricsdf$F1[i] <- temp_conmat[[4]][7]
  }
  
  modelmetricsdf.avg.F1<- mean(modelmetricsdf$F1)
  modelmetricsdf.avg.F1[is.na(modelmetricsdf.avg.F1)] <- 0
  

  metrics_plot <- ggplot(data = metricsdf, aes(x = Avg.F1)) + 
    geom_density(fill = "skyblue", alpha = 0.4, linewidth = 0.75) + 
    geom_vline(aes(xintercept = modelmetricsdf.avg.F1), 
               color = "red", 
               linetype = "dashed", 
               linewidth = 2) +
    labs(
      x = NULL,
      y = NULL
      ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 22, face = "bold", margin = ggplot2::margin(t = 10)),
      axis.title.y = element_text(size = 22, face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text = element_text(size = 20),  
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5)
    ) 
  
  
  print(metrics_plot)
}



distribution.graph <- function(clinrandom, clindata, fpanelrandom, fpaneldata, borutarandom, borutadata, logrrandom, logrdata, pos.outcome){
  
  random.distribution <- function(randomdata, pos.outcome){
    output <- data.frame(run = 1:25000,
                         F1 = NA)
    
    for (i in 1:25000) {
      test_conmat <- confusionMatrix(randomdata[[i]][[1]][[3]],
                                     positive = pos.outcome)
      
      output$F1[i] <- test_conmat[[4]][7]
      
      
    }
    output[is.na(output)] <- 0
    output
  }
  
  
  performance.distribution <- function(data, random, pos.outcome){
    
    temp_data2 <- data.frame(run = 1:25,
                             F1 = NA)
    
    for (i in 1:25) {
      temp_data <- confusionMatrix(data[[i]][[1]][[3]],
                                   positive = pos.outcome)
      
      
      temp_data2$F1[i] <- temp_data[[4]][7]
      
    }
    output <- data.frame(run = 1:25,
                         F1 = NA)
    
    for (i in 1:25) {
      output$F1[i] <- p.calc(temp_data2$F1[i], random)
      
    }
    
    output
  }
  
  #clinical set
  clinran <- random.distribution(clinrandom, pos.outcome)
  
  clin <- as.data.frame((clinran$F1 - min(clinran$F1)) /
                          (max(clinran$F1) - min(clinran$F1))) 
  
  colnames(clin)[1] <- "Range"
  clin$Model <- "Clinical Null"
  
  clinmodel <- performance.distribution(clindata, clinran$F1, pos.outcome)
  
  
  colnames(clinmodel)[2] <- "Range"
  clinmodel$Model <- "Clinical"
  
  
  
  
  
  #full panel set
  fpanelran <- random.distribution(fpanelrandom, pos.outcome)
  
  fpanel <- as.data.frame((fpanelran$F1 - min(fpanelran$F1)) /
                            (max(fpanelran$F1) - min(fpanelran$F1))) 
  
  colnames(fpanel)[1] <- "Range"
  fpanel$Model <- "Full Panel Null"
  
  fpanelmodel <- performance.distribution(fpaneldata, fpanelran$F1, pos.outcome)
  
  
  colnames(fpanelmodel)[2] <- "Range"
  fpanelmodel$Model <- "Full Panel"
  
  
  
  
  
  #boruta set
  borutaran <- random.distribution(borutarandom, pos.outcome)
  
  boruta <- as.data.frame((borutaran$F1 - min(borutaran$F1)) /
                            (max(borutaran$F1) - min(borutaran$F1))) 
  
  colnames(boruta)[1] <- "Range"
  boruta$Model <- "Boruta Null"
  
  borutamodel <- performance.distribution(borutadata, borutaran$F1, pos.outcome)
  
  
  colnames(borutamodel)[2] <- "Range"
  borutamodel$Model <- "Boruta"
  
  
  
  
  
  
  #logr set
  logrran <- random.distribution(logrrandom, pos.outcome)
  
  logr <- as.data.frame((logrran$F1 - min(logrran$F1)) /
                          (max(logrran$F1) - min(logrran$F1))) 
  
  colnames(logr)[1] <- "Range"
  logr$Model <- "LogR Null"
  
  logrmodel <- performance.distribution(logrdata, logrran$F1, pos.outcome)
  
  
  colnames(logrmodel)[2] <- "Range"
  logrmodel$Model <- "LogR"
  
  
  
  
  
  
  model_list <- list(clinmodel, fpanelmodel, borutamodel, logrmodel)
  
  #merge all data frames in list
  model_list <- model_list %>% reduce(full_join, by= c('Model', 'Range'))
  
  
  custom_order <- c("Clinical",
                    "Full Panel",
                    "Boruta",
                    "LogR")
  model_list$Model <- factor(model_list$Model, levels = custom_order)
  
  model_list <- model_list %>% replace(is.na(.), 0)
  
  ggplot(model_list, aes(x = Model, y = Range, fill = Model)) +
    geom_boxplot() +
    scale_y_continuous(limits = c(0, 1)) +
    # Assign colors in the exact order the groups will appear on the X-axis
    scale_fill_manual(values = c("blue", "red", "blue", "red")) +
    theme_minimal() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold", size = 22),
          axis.text.y = element_text(face = "bold", size = 18)
    ) 
  
  
}




double.graph <- function(leftrandom, leftdata, rightrandom, rightdata, pos.outcome){
  
  random.distribution <- function(randomdata, pos.outcome){
    output <- data.frame(run = 1:25000,
                         F1 = NA)
    
    for (i in 1:25000) {
      test_conmat <- confusionMatrix(randomdata[[i]][[1]][[3]],
                                     positive = pos.outcome)
      
      output$F1[i] <- test_conmat[[4]][7]
      
      
    }
    output[is.na(output)] <- 0
    output
  }
  
  
  performance.distribution <- function(data, random, pos.outcome){
    
    temp_data2 <- data.frame(run = 1:25,
                             F1 = NA)
    
    for (i in 1:25) {
      temp_data <- confusionMatrix(data[[i]][[1]][[3]],
                                   positive = pos.outcome)
      
      
      temp_data2$F1[i] <- temp_data[[4]][7]
      
    }
    output <- data.frame(run = 1:25,
                         F1 = NA)
    
    for (i in 1:25) {
      output$F1[i] <- p.calc(temp_data2$F1[i], random)
      
    }
    
    output
  }
  
  
  #Left graph set
  leftran <- random.distribution(leftrandom, pos.outcome)
  
  left <- as.data.frame((leftran$F1 - min(leftran$F1)) /
                          (max(leftran$F1) - min(leftran$F1))) 
  
  colnames(left)[1] <- "Range"
  left$Model <- "Clinical Null"
  
  leftmodel <- performance.distribution(leftdata, leftran$F1, pos.outcome)
  
  
  colnames(leftmodel)[2] <- "Range"
  leftmodel$Model <- "Clinical"
  
  
  
  
  
  
  
  #Right graph set
  rightran <- random.distribution(rightrandom, pos.outcome)
  
  right <- as.data.frame((rightran$F1 - min(rightran$F1)) /
                           (max(rightran$F1) - min(rightran$F1))) 
  
  colnames(right)[1] <- "Range"
  right$Model <- "Full Panel Null"
  
  rightmodel <- performance.distribution(rightdata, rightran$F1, pos.outcome)
  
  
  colnames(rightmodel)[2] <- "Range"
  rightmodel$Model <- "Full Panel"
  
  
  
  model_list <- list(left, leftmodel, right, rightmodel)
  
  #merge all data frames in list
  model_list <- model_list %>% reduce(full_join, by= c('Model', 'Range'))
  
  
  model_order <- c("Clinical Null",
                   "Clinical",
                   "Full Panel Null",
                   "Full Panel")
  
  model_list <- model_list %>% replace(is.na(.), 0)
  
  model_list$Model <- factor(model_list$Model, levels = model_order)
  
  
  
  ggplot(model_list, aes(x = Model, y = Range, fill = Model)) +
    geom_boxplot() +
    scale_fill_manual(values = c("skyblue", "red", "skyblue", "red")) +
    theme_minimal() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold", size = 26),
          axis.text.y = element_text(face = "bold", size = 24)
    )
  
}





single.graph <- function(leftrandom, rightdata, randomname, dataname, pos.outcome){
  
  random.distribution <- function(randomdata, pos.outcome){
    output <- data.frame(run = 1:25000,
                         F1 = NA)
    
    for (i in 1:25000) {
      test_conmat <- confusionMatrix(randomdata[[i]][[1]][[3]],
                                     positive = pos.outcome)
      
      output$F1[i] <- test_conmat[[4]][7]
      
      
    }
    output[is.na(output)] <- 0
    output
  }
  
  
  performance.distribution <- function(data, random, pos.outcome){
    
    temp_data2 <- data.frame(run = 1:25,
                             F1 = NA)
    
    for (i in 1:25) {
      temp_data <- confusionMatrix(data[[i]][[1]][[3]],
                                   positive = pos.outcome)
      
      
      temp_data2$F1[i] <- temp_data[[4]][7]
      
    }
    output <- data.frame(run = 1:25,
                         F1 = NA)
    
    for (i in 1:25) {
      output$F1[i] <- p.calc(temp_data2$F1[i], random)
      
    }
    
    output
  }
  
  
  #creating data objects
  leftran <- random.distribution(leftrandom, pos.outcome)
  
  left <- as.data.frame((leftran$F1 - min(leftran$F1)) /
                          (max(leftran$F1) - min(leftran$F1))) 
  
  colnames(left)[1] <- "Range"
  left$Model <- randomname
  
  
  rightmodel <- performance.distribution(rightdata, leftran$F1, pos.outcome)
  
  
  colnames(rightmodel)[2] <- "Range"
  rightmodel$Model <- dataname
  
  
  
  
  
  model_list <- list(left, rightmodel)
  
  #merge all data frames in list
  model_list <- model_list %>% reduce(full_join, by= c('Model', 'Range'))
  
  
  model_order <- c(randomname,
                   dataname)
  
  model_list <- model_list %>% replace(is.na(.), 0)
  
  model_list$Model <- factor(model_list$Model, levels = model_order)
  
  
  
  single.plot <- ggplot(model_list, aes(x = Model, y = Range, fill = Model)) +
    geom_boxplot() +
    scale_fill_manual(values = c("skyblue", "red")) +
    theme_minimal() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold", size = 26),
          axis.text.y = element_text(face = "bold", size = 24)
    )
  print(single.plot)
}









