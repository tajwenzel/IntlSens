#age.group.limits <- c(1,5,15,25,45,65)
#coverage dat is arranged on these intervals, risk groups (regular risk, increased risk)
#dates=as.Date(cov.eff.in[[1]]$V1, origin = "1970-01-01")
#coveffin<-cov.eff.in
#season<-1


function (riskratio, scenario, coveffin, season, new.cov) 
{
  
  cov <- coveffin[[season]][, 2:22]
  non <- rep(0,length(cov[[1]]))
  eff.pull <- coveffin[[season]][1, 23:43]
  dates <- as.Date(coveffin[[season]]$V1, origin = "1970-01-01")
  calendar <- matrix(rep(0), nrow = length(dates), ncol = 22)
  if (scenario == 1) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- non
    #calendar[t, c(15)] <- non
    #calendar[t, c(16)] <- non
    #calendar[t, c(17)] <- non
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- non
    calendar[t, c(4)] <- non
    calendar[t, c(5)] <- non
    calendar[t, c(6)] <- non
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  
  
  
  
  
  
  #PRESCHOOL
  if (scenario == 2) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]
    #calendar[t, c(2)] <- cov[[9]]
    #calendar[t, c(3)] <- cov[[9]]
    #calendar[t, c(4)] <- cov[[10]]
    #calendar[t, c(5)] <- cov[[10]]
    #calendar[t, c(6)] <- cov[[10]]
    #calendar[t, c(7)] <- cov[[11]]
    #calendar[t, c(8)] <- cov[[12]]
    #calendar[t, c(9)] <- cov[[13]]
    #calendar[t, c(10)] <- cov[[14]]
    #calendar[t, c(11)] <- cov[[14]]
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- new.cov[[season]]#toddlers
    #calendar[t, c(15)] <- non
    #calendar[t, c(16)] <- non
    #calendar[t, c(17)] <- non
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- new.cov#toddlers
    calendar[t, c(4)] <- non
    calendar[t, c(5)] <- non
    calendar[t, c(6)] <- non
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  #PRIMARY
  if (scenario == 3) {
    t = 1:length(dates)
    
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- non
    #calendar[t, c(15)] <- new.cov[[season]]
    #calendar[t, c(16)] <- non
    #calendar[t, c(17)] <- non
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]][[season]][1:length(t)]
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- non
    calendar[t, c(4)] <- new.cov#primary
    calendar[t, c(5)] <- non
    calendar[t, c(6)] <- non
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  #SECONDARY SCHOOL
  if (scenario == 4) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- non
    #calendar[t, c(15)] <- non
    #calendar[t, c(16)] <- new.cov[[season]]
    #calendar[t, c(17)] <- new.cov[[season]]
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- non
    calendar[t, c(4)] <- non
    calendar[t, c(5)] <- new.cov
    calendar[t, c(6)] <- new.cov
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  #PRESCHOOL+PRIMARY
  if (scenario == 5) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- new.cov[[season]]
    #calendar[t, c(15)] <- new.cov[[season]]
    #calendar[t, c(16)] <- non
    #calendar[t, c(17)] <- non
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- new.cov
    calendar[t, c(4)] <- new.cov
    calendar[t, c(5)] <- non
    calendar[t, c(6)] <- non
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
    
   
  }
  #PRESCHOOL+PRIMARY+SECONDARY
  if (scenario == 6) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- new.cov[[season]]
    #calendar[t, c(15)] <- new.cov[[season]]
    #calendar[t, c(16)] <- new.cov[[season]]
    #calendar[t, c(17)] <- new.cov[[season]]
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- new.cov
    calendar[t, c(4)] <- new.cov
    calendar[t, c(5)] <- new.cov
    calendar[t, c(6)] <- new.cov
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  #PRESCHOOL+SECONDARY
  if (scenario == 7) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- new.cov[[season]]
    #calendar[t, c(15)] <- non
    #calendar[t, c(16)] <- new.cov[[season]]
    #calendar[t, c(17)] <- new.cov[[season]]
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
  
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- new.cov
    calendar[t, c(4)] <- non
    calendar[t, c(5)] <- new.cov
    calendar[t, c(6)] <- new.cov
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  #PRIMARY+SECONDARY
  if (scenario == 8) {
    t = 1:length(dates)
    #calendar[t, c(1)] <- cov[[8]]#0-1
    #calendar[t, c(2)] <- cov[[9]]#1-2
    #calendar[t, c(3)] <- cov[[9]]#2-4
    #calendar[t, c(4)] <- cov[[10]]#5-11
    #calendar[t, c(5)] <- cov[[10]] #12-14
    #calendar[t, c(6)] <- cov[[10]] #15-16
    #calendar[t, c(7)] <- cov[[11]] #17-24
    #calendar[t, c(8)] <- cov[[12]] #25-44
    #calendar[t, c(9)] <- cov[[13]]#45-64
    #calendar[t, c(10)] <- cov[[14]]#65-74
    #calendar[t, c(11)] <- cov[[14]] #75+
    #calendar[t, c(12)] <- non
    #calendar[t, c(13)] <- non
    #calendar[t, c(14)] <- non
    #calendar[t, c(15)] <- new.cov[[season]]
    #calendar[t, c(16)] <- new.cov[[season]]
    #calendar[t, c(17)] <- new.cov[[season]]
    #calendar[t, c(18)] <- non
    #calendar[t, c(19)] <- non
    #calendar[t, c(20)] <- non
    #calendar[t, c(21)] <- cov[[7]]
    #calendar[t, c(22)] <- cov[[7]]
    
    calendar[t, c(12)] <- cov[[8]]#0-1
    calendar[t, c(13)] <- cov[[9]]#1-2
    calendar[t, c(14)] <- cov[[9]]#2-4
    calendar[t, c(15)] <- cov[[10]]#5-11
    calendar[t, c(16)] <- cov[[10]] #12-14
    calendar[t, c(17)] <- cov[[10]] #15-16
    calendar[t, c(18)] <- cov[[11]] #17-24
    calendar[t, c(19)] <- cov[[12]] #25-44
    calendar[t, c(20)] <- cov[[13]]#45-64
    calendar[t, c(21)] <- cov[[14]]#65-74
    calendar[t, c(22)] <- cov[[14]] #75+
    calendar[t, c(1)] <- non
    calendar[t, c(2)] <- non
    calendar[t, c(3)] <- non
    calendar[t, c(4)] <- new.cov
    calendar[t, c(5)] <- new.cov
    calendar[t, c(6)] <- new.cov
    calendar[t, c(7)] <- non
    calendar[t, c(8)] <- non
    calendar[t, c(9)] <- non
    calendar[t, c(10)] <- cov[[7]]
    calendar[t, c(11)] <- cov[[7]]
  }
  
  #if(length(unique(as.numeric(calendar[1,])))==1)
  if(dates[1]<=as.Date(paste0(format(as.Date(dates[1], format="%d/%m/%Y"),"%Y"),'-09','-01')))
  {
   xx<-matrix(calendar[dates<=as.Date(paste0(format(as.Date(dates[1], format="%d/%m/%Y"),"%Y"),'-09','-01')),], ncol=22)
   xx[xx>0]<-0
   calendar[1:dim(xx)[1],]<-xx 
   
   vaccine3 <- as_vaccination_calendar(efficacy = eff.pull[1:7] ,
                                      dates = dates, coverage = calendar, no_risk_groups = 2, 
                                      no_age_groups = 11)
  }else{
    
    dates<-c(as.Date(paste0(format(as.Date(dates[1], format="%d/%m/%Y"),"%Y"),'-09','-01')), dates)
    calendar<-rbind(rep(0,22), calendar)
    vaccine3 <- as_vaccination_calendar(efficacy = eff.pull[1:7] ,
                                        dates = dates, coverage = calendar, no_risk_groups = 2, 
                                        no_age_groups = 11)
    
  }
  
  
  eff.out <- as.vector(unlist(vaccine3$efficacy), mode = "numeric")
  v.output <- list(efficacy = eff.out, dates = vaccine3$dates, 
                   calendar = vaccine3$calendar)
  return(v.output)
}
