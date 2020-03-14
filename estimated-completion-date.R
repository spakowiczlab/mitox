# estimate-completion-date.R
# 
# This script estimates the completion date of a clinical trial from a few 
# parameters that can be captured as the trial progresses, namely the 
# number of patients currently accrued and the current rate of accrual.
# It was written to support the mitox (2019C0151) study but the function 
# inputs should be sufficiently general to accommodate many studies. 


library(lubridate)

estimatedCompletion <- function(todays.date, 
                                weekly.consent.rate, 
                                accrual.target, 
                                current.accrual,
                                study.duration) {
  
  today <- as.Date(todays.date, format = "%F")
  
  patients.to.go <- accrual.target - current.accrual
  
  weeks.to.full.accrual <- ceiling(patients.to.go / weekly.consent.rate)
  
  weeks.to.completion <- weeks.to.full.accrual + study.duration
  
  date.full.accrual <- today + lubridate::weeks(weeks.to.full.accrual)
  
  date.study.completion <- today + lubridate::weeks(weeks.to.completion)
  
  cat(paste0("Date of full accrual = ", date.full.accrual),
      paste0("Date of study completion = ", date.study.completion))
}

estimatedCompletion(todays.date = "2020-03-13",
                    weekly.consent.rate = 1.25,
                    accrual.target = 63,
                    current.accrual = 5,
                    study.duration = 16)

