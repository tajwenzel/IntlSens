
#Function for calculating number of vaccine doses per season per strategy

alt.program.coverage<-function(new.coverage,strategy, strain)
{
  
  #tx_start=as.Date(c("1995-10-01","1996-10-01","1997-10-01","1998-10-01","1999-10-01","2000-10-01",
   #                  "2001-10-01","2002-10-01","2003-10-01","2004-10-01","2005-10-01","2006-10-01",
    #                 "2007-10-01","2008-10-01","2009-09-28","2010-09-13","2011-09-19","2012-08-20",
     #                "2013-08-26"),origin="1970-01-01")
  
  if(strain==1) {cov.effH1N1<-coverageH1};
  if(strain==2) {cov.effH3N2<-coverageH3};
  if(strain==3) {cov.effB<-coverageB};
  
  #getting actual coverage, multiply days between "rates" remember.
  cov.eff<-get(paste0('cov.eff', strain.name[strain]))
  cov.eff.new<-NULL
  intervention.cov<-NULL
  for(years in 1:length(cov.eff))
  {
    log.cov<-NULL
    log.cov<-predict(new.coverage, newdata=data.frame(cov.eff[[years]]$V1+(15584-cov.eff[[years]]$V1[1])))
    
    tt<-vstrategy(riskratio = risk.ratios.ce, scenario=strategy, coveffin = cov.eff, season = years, new.cov = log.cov)
    days<-as.vector(diff(as.matrix(tt$dates))) #difference in days between each row, differs between years
    start<-as.numeric(tt$dates[1])
    overl<-days*tt$calendar[1:length(days),1:22]
    #Now that coverage at each moment in time calcualted; cumulative coverage should be calculated and arrange everything into tables. 
    overl[overl>1]<-1;
    #cov.eff.new<-NULL;
    #filler<-c(rep(0,22))
    #cov.eff.new[[years]]<-matrix(rep(0,dim(cov.eff[[years]])[1]*23),nrow=dim(cov.eff[[years]])[1]);
    #cov.eff.new[[years]][1,]<-rbind(c(start,filler));
    
    #rearrange
    #cov.eff.new[[years]][1:dim(cov.eff[[years]])[1],2:23]<-as.matrix(overl[1:length(days),]);
    #cov.eff.new[[years]]<-cumsum(as.data.frame(cov.eff.new[[years]][,2:23]));
    
    cov.eff.new[[years]]<-cbind.data.frame(tt$dates, rbind(rep(0, 22), cumsum(as.data.frame(overl[,1:22]))));
    #get last row
    intervention.cov[[years]]<-cov.eff.new[[years]][dim(cov.eff.new[[years]])[1],];
  }
  
  alt.program<<-matrix(unlist(intervention.cov), nrow=length(cov.eff), byrow=T)
  #save(intervention.cov,file=paste0('intervention.coverage',scenario,unames[s])) 
}

