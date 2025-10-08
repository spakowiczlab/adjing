




p.calc <- function(data, random){
  avg <- (data)
  p <- ecdf(random)
  p(avg)
  
}




svm.metrics <- function(dataobject, output, pos.outcome){
  test_conmat2 <- data.frame(run = 1:25000,
                             F1 = NA)
  
  x = 1
  output <- data.frame(run = 1:1000,
                       Avg.F1 = NA)
  
  for (i in 1:1000) {
    y = x+24
    
    output$Avg.F1[i] <- mean(dataobject$F1[x:y])
    
    x = x+25
  }
  output
}

svm.p.calc <- function(randomdata, pos.outcome, output, dataobject){
  #function that generates averages of metrics for randomized run
  random_met <- svm.metrics(randomdata,
                            random_met,
                            pos.outcome)
  
  #storing experimental run additional metrics   
  temp_conmat2 <- data.frame(run = 1:25,
                             F1 = NA)
  
  #storing percentiles    
  output <- data.frame(
                       F1 = NA)
  
  #converting NAs to 0 so we don't skew the averages    
  # temp_conmat2[is.na(temp_conmat2)] <- 0
  
  #percentile calculation for various metrics  
  
  avg.F1 <- mean(temp_conmat2$F1)
  output$F1[1] <- p.calc(avg.F1,
                         random_met$Avg.F1)
  
  data_wide <- full_join(random_met, temp_conmat2)
  data_wide <- select(data_wide, -run)
  data_graph <- data_wide %>% 
    gather(key="RunVersion", value="Val") %>%
    ggplot( aes(x=RunVersion, y=Val, fill=RunVersion)) +
    geom_violin()
  
  print(data_graph)
  
  output
}

