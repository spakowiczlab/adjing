

outcome.creation <- function(cohort, pos.outcome){
  filter(pt.status, 
         Cohort == cohort, 
         Reason.for.discontinuation == "Therapy complete" | 
           Reason.for.discontinuation == pos.outcome)   
}


p.calc <- function(data, random){
  avg <- (data)
  p <- ecdf(random)
  p(avg)
  
}




