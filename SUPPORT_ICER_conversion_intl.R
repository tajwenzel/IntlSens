CONVERT_INCIDENCE_TO_HEALTH_OUTCOME<<-function(incidence1,incidence2,outcome,strain)
{
  #restructure<-function(single.scenario) 
  #{
  # age_group_CROMER <- c("0 - 6 m", "6 m - 4 y", "5 - 14 y", 
  #"15 - 44 y", "45 - 64 y", "65+ y")
  #risk_group <- c("low risk", "high risk")
  #tmp <- single.scenario
  #regroup.age <- data.frame(matrix(rep(0, 12 * 1000), nrow = 1000, 
  # byrow = TRUE))
  #regroup.age[, 1] <- tmp[, 11]/2
  #regroup.age[, 2] <- tmp[, 11]/2 +rowSums(tmp[, 12:13])
  #regroup.age[, 3] <- rowSums(tmp[, 14:15])
  #regroup.age[, 4] <- rowSums(tmp[, 16:19])
  #regroup.age[, 5] <- tmp[, 20]
  #regroup.age[, 6] <- rowSums(tmp[, c(21, 22)])
  #regroup.age[, 7] <- tmp[, 1]/2
  #regroup.age[, 8] <- tmp[, 1]/2 + rowSums(tmp[, 2:3])
  #regroup.age[, 9] <- rowSums(tmp[,4:5])
  #regroup.age[, 10] <- rowSums(tmp[, 6:9])
  #regroup.age[, 11] <- tmp[, 10]
  #regroup.age[, 12] <- rowSums(tmp[, c(11, 12)])
  #regroup.age <- as.data.frame(regroup.age[, 1:12])
  #colnames(regroup.age) <- paste(rep(age_group_CROMER, 2), 
  #                              rep(risk_group, each = length(age_group_CROMER)))
  #return(regroup.age)
  #}
  
  #incidence1<-tableH3[[1]]
  #incidence2<-tableH3[[2]]
  
  
  #checking list length
  inc.length<-length(incidence1)
  missing<-length(incidence1[1][sapply(incidence1,is.null)]) #number of missing
  tab.diff<<-inc.length-missing
  tab.temp1<-incidence1[1:tab.diff]
  tab.temp2<-incidence2[1:tab.diff]
  
  #strain<-c('H1N1','H3N2','B')
  setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
  if(strain==1) {load(death.risk.tables[2]);load(hosp.risk.tables[2]);load(GP.risk.tables[2])}
  if(strain==2) {load(death.risk.tables[3]);load(hosp.risk.tables[3]);load(GP.risk.tables[3])}
  if(strain==3) {load(death.risk.tables[1]);load(hosp.risk.tables[1]);load(GP.risk.tables[1])}
  
  #Computes the number of symptomatic/febrile cases from incidence table
  if(outcome=="cases")
  {
    risk.sample<-rtriangle(n.samples,a=0.309,b=0.513,c=0.396) 
    #generate the percentile of febrile cases from triangular distribution
  } else { 
    #other outcomes are calculated using a risk table saved in RSave
    risk.tab <- apportion(loadRData(paste0(paste('tab_risk', outcome, strain.name[strain],sep="_"),'.R')))
    if(dim(incidence1[[1]])[1]==2500)
    {risk.sample<-rbind(risk.tab[sample.int(n.samples/2.5, n.samples/2.5),], risk.tab[sample.int(n.samples/2.5, n.samples/2.5),], risk.tab[sample.int(n.samples/5, n.samples/5),])}
    if(dim(incidence1[[1]])[1]==3000)
    {risk.sample<-rbind(risk.tab[sample.int(n.samples/3, n.samples/3),], risk.tab[sample.int(n.samples/3, n.samples/3),], risk.tab[sample.int(n.samples/3, n.samples/3),])
    }else{
      risk.sample<-rbind(risk.tab[sample.int(n.samples/2, n.samples/2),], risk.tab[sample.int(n.samples/2, n.samples/2),])}
    risk.sample.ages<-apportion(risk.sample)
  }
  
  risk1<-lapply(1:tab.diff, FUN = function(i) tab.temp1[[i]][1:n.samples,]*risk.sample)
  risk2<-lapply(1:tab.diff, FUN = function(i) tab.temp2[[i]][1:n.samples,]*risk.sample)
  return(list(risk1, risk2))
}

CONVERT_INCIDENCE_TO_QALY<<-function(incidence1,incidence2,strain,discount)
{
  #QALY loss from febrile cases
  QALY.loss.AJ<-read.csv(paste('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/',"DATA/AJ_QALY_list.csv",sep="/"))
  QALY.loss.cases<-sample(QALY.loss.AJ$v,n.samples,replace=TRUE)
  
  #QALY loss from hospitalisations, using distributions from Marc's paper
  QALY.loss.hosp<-rnorm(n.samples,mean=0.018,sd=0.0018); rownames(QALY.loss.hosp)<-NULL
  
  #Utilize conversion function. Calculate the different health outcomes
  table.cases<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="cases",strain); rownames(table.cases)<-NULL
  table.hosp<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="hosp",strain); rownames(table.hosp)<-NULL
  table.death<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="death",strain); rownames(table.death)<-NULL
  table.GP<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="GP",strain=strain); rownames(table.GP)<-NULL
  
  if(discount==3.5) death.QALY.loss<-apportion(c(22.61,23.00,22.32,18.84,12.71,6.52,24.07,24.05,23.18,20.33,14.62,7.28))
  if(discount==1.5) death.QALY.loss<-apportion(c(37.06,37.61,35.19,27.16,16.09,7.41,40.77,40.62,37.84,30.40,19.10,8.40))
  if(discount==0) death.QALY.loss<-apportion(c(60.93,61.62,54.96,38.23,19.75,8.24,69.96,69.45,61.64,44.47,24.14,9.46))
  
  QALY.loss1<-NULL
  #For Status QUo vaccine program
  QALY.cases1<-lapply(1:tab.diff, FUN = function(i) sweep(table.cases[[1]][[i]],MARGIN=1,QALY.loss.cases,'*'))
  QALY.hosp1<-lapply(1:tab.diff, FUN = function(i) sweep(table.hosp[[1]][[i]],MARGIN=1,QALY.loss.hosp,'*'))
  QALY.death1<-lapply(1:tab.diff, FUN = function(i) sweep(table.death[[1]][[i]],MARGIN=2,death.QALY.loss,'*'))
  #QALY.GP1<-lapply(1:tab.diff, FUN = function(i) sweep(table.GP[[1]][[i]],MARGIN=2,death.QALY.loss,'*'))
 
  
  
  QALY.loss2<-NULL
  #For Status QUo vaccine program
  QALY.cases2<-lapply(1:tab.diff, FUN = function(i) sweep(table.cases[[2]][[i]],MARGIN=1,QALY.loss.cases,'*'))
  QALY.hosp2<-lapply(1:tab.diff, FUN = function(i) sweep(table.hosp[[2]][[i]],MARGIN=1,QALY.loss.hosp,'*'))
  QALY.death2<-lapply(1:tab.diff, FUN = function(i) sweep(table.death[[2]][[i]],MARGIN=2,death.QALY.loss,'*'))
  
  QALY.loss1<-Map('+',Map('+',QALY.cases1,QALY.hosp1),QALY.death1)
  QALY.loss2<-Map('+',Map('+',QALY.cases2,QALY.hosp2),QALY.death2)
  
  return(list(QALY.loss1,QALY.loss2,QALY.cases1,QALY.cases2,QALY.hosp1,QALY.hosp2,QALY.death1,QALY.death2))
}

QALYdifferences <<- function(incidence1,incidence2,strain,discount)
{
  temp.sample.QALY <- CONVERT_INCIDENCE_TO_QALY(incidence1,incidence2,strain,discount);
  return(Map('-',temp.sample.QALY[[1]],temp.sample.QALY[[2]])) #should return number of QALY's gained in new program
}

GPHospAvert<<- function(incidence1,incidence2,strain){
  
  table.GP<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="GP",strain=strain)
  table.hosp<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="hosp",strain=strain)
  table.death<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="death",strain=strain)
  table.cases<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="cases",strain=strain)
  
  #GP cases averted by new program
  Diff.GP <-mapply('-',table.GP[[1]],table.GP[[2]],SIMPLIFY=FALSE);
  
  #hospital cases averted by new program
  Diff.Hosp <- mapply('-',table.hosp[[1]],table.hosp[[2]],SIMPLIFY=FALSE);
  
  #formatted in a list with 2 pieces, first is status quo and second is intervention
  avg.GP1<-table.GP[[1]];
  avg.GP2<-table.GP[[2]];
  avg.hosp1<-table.hosp[[1]];
  avg.hosp2<-table.hosp[[2]];
  avg.death1<-table.death[[1]];
  avg.death2<-table.death[[2]];
  avg.cases1<-table.cases[[1]];
  avg.cases2<-table.cases[[2]];
  return(list(Diff.GP,Diff.Hosp,avg.GP1,avg.GP2,avg.hosp1,avg.hosp2, avg.death1,avg.death2, avg.cases1, avg.cases2))
}

#CONVERT_INCIDENCE_TO_QALY(tableB[[1]],tableB[[3]],'B',3.5)
# 6.1 Generate QALY gain between two interventions
#DataSet1 = QALY loss without additional vac.
#DataSet2 = QALY loss with addtional vac.
#strain = "H1N1"/"H3N2"/"B"

#GPHospAvert(tableH3[1],tableH3[2],2)
#strain<-2
#incidence1<-tableH3[[1]]
#incidence2<-tableH3[[2]]

hilo.structure<-function(keep) # Function to restructure model output into data age groups
{
  if(is.vector(keep)==TRUE) {keep<-t(as.matrix(keep))}; #must be matrix format to work
  #Convert age groups and risk groups
  #age.group.limits<<-c(1,2,5,12,15,17,25,45,65,75) #upper limits
  c1 <- keep[,1]
  c2 <- rowSums(keep[,2:3])
  c3 <- rowSums(keep[,4:5])
  c4 <- rowSums(keep[,6:8])
  c5 <- keep[,9]
  c6 <- rowSums(keep[,10:11])
  c7 <- keep[,12]
  c8 <- rowSums(keep[,13:14])
  c9 <- rowSums(keep[,15:16])
  c10 <- rowSums(keep[,17:19])
  c11 <- keep[,20]
  c12 <-rowSums(keep[,21:22])
  regrouped <- cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
  return(regrouped)
}




intl.AgeStuc.CEA.outcomes<-function(Dataset1, Dataset2, Dataset3,Dataset4, strain,i.country, dcount, program1, program2, program3,program4, i.cov, version)
{
  #names
  fname<-c('BE', 'DE', 'FI', 'GB', 'IT', 'LU', 'NL', 'PL', 'PE', 'FR', 'ZI')
  strain.name<-c('H1N1','H3N2','B')
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
  ##cut dataset to shortest if they are different lengths
  recut<-c(dim(Dataset1[[1]])[1], dim(Dataset2[[1]])[1], dim(Dataset3[[1]])[1],dim(Dataset4[[1]])[1] )
  n.samples<<-recut[which.min(recut)]
  tab.diff<-length(Dataset1)
  
  
  ##make sure the right years are counted
  if(tab.diff>14)
  {
   Dataset1<-Dataset1[1:14];    Dataset3<-Dataset3[1:14];
   Dataset2<-Dataset2[1:14];    Dataset4<-Dataset4[1:14]; 
   tab.diff<-14
  }

  #saving names
  ppname<-c("StatusQuo","Preschool","Primary","Secondary","Preschool+Primary School","Preschool+Primary+Secondary","Preschool+Secondary","Primary+Secondary")
  
  
  #####Vaccination cost and dosage calculation
  population<-matrix(popv[1:22], nrow=1, byrow=TRUE)
  status.quo.program<-alt.program.coverage(new.coverage = i.cov, strategy = program1, strain = strain)#status quo
  new.program2<-alt.program.coverage(new.coverage = i.cov,strategy = program2, strain = strain)
  new.program3<-alt.program.coverage(new.coverage = i.cov,strategy = program3, strain = strain)
  new.program4<-alt.program.coverage(new.coverage = i.cov,strategy = program4, strain = strain)
  
  
  #returns number of vaccine doses per season
  vaccines.doses.new4<-t(t(new.program4[,2:23])*as.vector(popv[1:22]))
  vaccines.doses.new2<-t(t(new.program2[,2:23])*as.vector(popv[1:22]))
  vaccines.doses.new3<-t(t(new.program3[,2:23])*as.vector(popv[1:22]))
  vaccines.doses.statquo<-t(t(status.quo.program[,2:23])*as.vector(popv[1:22]))
  
  ################ Mean and SD function 
  estfunc<-function(x) {
    p1<-mean(x); p2<-quantile(x,0.025); p3<-quantile(x,0.975); p4<-sd(x)/sqrt(length(x))
    return(data.frame(p1,p2,p3,p4))
  }
  ########## 7.1.2 #calculates cost per dose########
  vaccine.cost.statquo<-vaccine.cost2<-vaccine.cost3<-vaccine.cost4<-list()
  for(tty in 1:dim(vaccines.doses.statquo)[1])
  {
    rcost.LAIV.school<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples,a=17,b=25,c=20.14), nrow=n.samples) #cost of vaccine per dose per season
    #L7.24+L2.80+L7.64 (reimbursement, service, dose cost) avg nurse 2.80 pounds/100 doses at rate of 36lbs/hr and 20 doses per hr works 7.5 hours a day, 2.5 travel and set up time and lunch
    rcost.GP<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples,a=17,b=25,c=19.66) , nrow=n.samples)
    #L7.24+L2.25+L7.64 (reimbursement, service, dose cost)
    rcost.pharmacy<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples, a=14, b=22, c=17.29), nrow=n.samples)
    #Katie paper BMJ pharmacy
    
    prop65up.pharm<-0.065 
    prop2to4.pharm<-0
    propadults.LR.pharm<-0.055
    propadults.HRpharm<-0.056
    
    #vaccine cost statquo 
    vaccine.cost.statquo[[tty]] <-cbind(t(vaccines.doses.statquo[tty,1:3]*t(rcost.GP[,1:3])), 
                                        t(vaccines.doses.statquo[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                        t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.statquo[tty,7:9]+
                                    t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.statquo[tty,7:9]), 
                                        t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.statquo[tty,10:11]+
                                        t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.statquo[tty,10:11]),
                    t(vaccines.doses.statquo[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+t(vaccines.doses.statquo[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
    
    #vaccine cost new 
    vaccine.cost2[[tty]] <-cbind(t(vaccines.doses.new2[tty,1:3]*t(rcost.GP[,1:3])), 
                                t(vaccines.doses.new2[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.new2[tty,7:9]+
                                    t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.new2[tty,7:9]), 
                                t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.new2[tty,10:11]+
                                    t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.new2[tty,10:11]),
                                t(vaccines.doses.new2[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+
                                  t(vaccines.doses.new2[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
    
    vaccine.cost3[[tty]] <-cbind(t(vaccines.doses.new3[tty,1:3]*t(rcost.GP[,1:3])), 
                                 t(vaccines.doses.new3[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                 t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.new3[tty,7:9]+
                                     t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.new3[tty,7:9]), 
                                 t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.new3[tty,10:11]+
                                     t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.new3[tty,10:11]),
                                 t(vaccines.doses.new3[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+
                                   t(vaccines.doses.new3[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
    
    vaccine.cost4[[tty]] <-cbind(t(vaccines.doses.new4[tty,1:3]*t(rcost.GP[,1:3])), 
                                 t(vaccines.doses.new4[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                 t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.new4[tty,7:9]+
                                     t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.new4[tty,7:9]), 
                                 t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.new4[tty,10:11]+
                                     t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.new4[tty,10:11]),
                                 t(vaccines.doses.new4[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+
                                   t(vaccines.doses.new4[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
  }
  
  
  ################ Begin calculation of outcomes-----
  
  
  #incidence of symptomaic and asymptomatic cases
  cumi1<-Reduce('+', Dataset1)/length(Dataset1)
  cumi2<-Reduce('+', Dataset2)/length(Dataset2)
  cumi3<-Reduce('+', Dataset3)/length(Dataset3)
  cumi4<-Reduce('+', Dataset4)/length(Dataset4)
  
  
  ##sum.outcomes outputs as a list from GPHospAvert. This is list has 10 pieces organized as list(Diff.GP,Diff.Hosp,avg.GP1,avg.GP2,avg.hosp1,avg.hosp2,avg.death1,avg.death2, avg.cases1, avg.cases2)
  
  sum.outcomes <- GPHospAvert(Dataset1,Dataset2, strain=strain);
  sum.outcomes2 <- GPHospAvert(Dataset1,Dataset3, strain=strain);
  sum.outcomes3 <- GPHospAvert(Dataset1,Dataset4, strain=strain);
  
  
  #Simulation output of dimension number simulations * number of seasons
  
  #febrile cases, not infections
  #all.cases1<-sum.outcomes[[9]]
  #all.cases2<-sum.outcomes[[10]]
  #all.cases3<-sum.outcomes2[[10]]
  #all.cases4<-sum.outcomes3[[10]]
  
  
  #calc.cases<-list(all.cases1, all.cases2, all.cases3, all.cases4) #febrile cases
  
  
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #save(calc.cases, file=paste0(fname[i.country],'calc.cases.strat',strain.name[strain],round(last(last(i.cov)),digits = 2),'p', program1, program2, program3, program4,'d', dcount)) 
  
  
  ####This series finds the mean across all simulations for each season and condenses it in to a matrix where the rows are s and the columns are age groups
  
  #avg.cases1<-Reduce('+', all.cases1)/length(all.cases1) #per n.sample simulation
  #avg.cases2<-Reduce('+', all.cases2)/length(all.cases2)
  #avg.cases3<-Reduce('+', all.cases3)/length(all.cases3)
  #avg.cases4<-Reduce('+', all.cases4)/length(all.cases4)
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
  
  aged.cases1<-Reduce('+',sum.outcomes[[9]])/tab.diff
  aged.cases2<-Reduce('+',sum.outcomes[[10]])/tab.diff
  aged.cases3<-Reduce('+',sum.outcomes2[[10]])/tab.diff
  aged.cases4<-Reduce('+',sum.outcomes3[[10]])/tab.diff
  aged.cases<-list(aged.cases1, aged.cases2, aged.cases3, aged.cases4)
  save(aged.cases, file=paste0(fname[i.country],'agestrat.cases',strain.name[strain],'all.programs',version))
  
  
  #number of GP consults
  aged.GP1<-Reduce('+',sum.outcomes[[3]])/tab.diff
  aged.GP2<-Reduce('+',sum.outcomes[[4]])/tab.diff
  aged.GP3<-Reduce('+',sum.outcomes2[[4]])/tab.diff
  aged.GP4<-Reduce('+',sum.outcomes3[[4]])/tab.diff
  aged.GP<-list(aged.GP1, aged.GP2, aged.GP3, aged.GP4)
  save(aged.GP, file=paste0(fname[i.country],'agestrat.GP',strain.name[strain],'all.programs',version))
  
  
  
  #Number of hospitalizations 
  
  aged.hosp1<-Reduce('+',sum.outcomes[[5]])/tab.diff
  aged.hosp2<-Reduce('+',sum.outcomes[[6]])/tab.diff
  aged.hosp3<-Reduce('+',sum.outcomes2[[6]])/tab.diff
  aged.hosp4<-Reduce('+',sum.outcomes3[[6]])/tab.diff
  aged.hosp<-list(aged.hosp1, aged.hosp2, aged.hosp3, aged.hosp4)
  save(aged.hosp, file=paste0(fname[i.country],'agestrat.hosp',strain.name[strain],'all.programs',version))
  
  
  
  #number of deaths
  aged.death1<-Reduce('+',sum.outcomes[[7]])/tab.diff
  aged.death2<-Reduce('+',sum.outcomes[[8]])/tab.diff
  aged.death3<-Reduce('+',sum.outcomes2[[8]])/tab.diff
  aged.death4<-Reduce('+',sum.outcomes3[[8]])/tab.diff
  aged.death<-list(aged.death1, aged.death2, aged.death3, aged.death4)
  save(aged.death, file=paste0(fname[i.country],'agestrat.death',strain.name[strain],'all.programs',version))
  
  
  
  ###############################################################################
  ####### CALCULATE COST from medical GP and Hospitalized cases averted ########
  ##############################################################################
  
  #cost of GP usage in new program
  gp.cost.dist<-matrix(rnorm(n.samples, mean = 39, sd = 8.6),ncol=1)*1.046 #sample from cost distribution
  gp.cost.mat1<-matrix(rep(gp.cost.dist,dim(sum.outcomes[[3]][[1]])[2]), ncol=dim(sum.outcomes[[4]][[1]])[2])
  gp.cost.mat2<-matrix(rep(gp.cost.dist,dim(sum.outcomes[[4]][[1]])[2]), ncol=dim(sum.outcomes[[4]][[1]])[2])
  gp.cost.mat3<-gp.cost.mat4<-matrix(rep(gp.cost.dist,dim(sum.outcomes2[[4]][[1]])[2]), ncol=dim(sum.outcomes2[[4]][[1]])[2])
  gp.cost.complete1<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[3]][[i]]*gp.cost.mat1))/tab.diff
  gp.cost.complete2<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[4]][[i]]*gp.cost.mat2))/tab.diff
  gp.cost.complete3<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[4]][[i]]*gp.cost.mat3))/tab.diff
  gp.cost.complete4<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[4]][[i]]*gp.cost.mat4))/tab.diff
  
  ##########################Hospitalization costs
  
  #cost of hospitalization for new program
  hosp.cost.dist<-matrix(rnorm(n.samples, mean = 911, sd = 215),ncol=1)*1.085
  hosp.cost.mat2<-hosp.cost.mat3<-hosp.cost.mat1<-hosp.cost.mat4<-matrix(rep(hosp.cost.dist,dim(sum.outcomes[[5]][[1]])[2]), ncol=dim(sum.outcomes[[5]][[1]])[2])
  
  hosp.cost.complete2<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[6]][[i]]*hosp.cost.mat2))/tab.diff
  hosp.cost.complete3<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[6]][[i]]*hosp.cost.mat3))/tab.diff
  hosp.cost.complete1<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[5]][[i]]*hosp.cost.mat1))/tab.diff
  hosp.cost.complete4<-Reduce('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[6]][[i]]*hosp.cost.mat4))/tab.diff
  
  
  #healthcare utilization costs
  avg.hc.cost1<-Reduce('+',(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[5]][[i]]*hosp.cost.mat1),lapply(1:tab.diff, FUN = function(i) sum.outcomes[[3]][[i]]*gp.cost.mat1))))/tab.diff
  avg.hc.cost2<-Reduce('+',(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[6]][[i]]*hosp.cost.mat2),lapply(1:tab.diff, FUN = function(i) sum.outcomes[[4]][[i]]*gp.cost.mat2))))/tab.diff
  avg.hc.cost3<-Reduce('+',(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[6]][[i]]*hosp.cost.mat3),lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[4]][[i]]*gp.cost.mat3))))/tab.diff
  avg.hc.cost4<-Reduce('+',(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[6]][[i]]*hosp.cost.mat4),lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[4]][[i]]*gp.cost.mat4))))/tab.diff
  
  
  hc.cost1<-Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[5]][[i]]*hosp.cost.mat1),lapply(1:tab.diff, FUN = function(i) sum.outcomes[[3]][[i]]*gp.cost.mat1))
  hc.cost2<-(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[6]][[i]]*hosp.cost.mat2),lapply(1:tab.diff, FUN = function(i) sum.outcomes[[4]][[i]]*gp.cost.mat2)))
  hc.cost3<-(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[6]][[i]]*hosp.cost.mat3),lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[4]][[i]]*gp.cost.mat3)))
  hc.cost4<-(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[6]][[i]]*hosp.cost.mat4),lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[4]][[i]]*gp.cost.mat4)))
  
  tot.hccost1<-Reduce('+', hc.cost1)/tab.diff
  tot.hccost2<-Reduce('+', hc.cost2)/tab.diff
  tot.hccost3<-Reduce('+', hc.cost3)/tab.diff
  tot.hccost4<-Reduce('+', hc.cost4)/tab.diff
  
  net.cost12<-lapply(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes[[2]][[i]]*hosp.cost.mat2),lapply(1:tab.diff, FUN = function(i) sum.outcomes[[1]][[i]]*gp.cost.mat2)), rowSums)
  net.cost13<-lapply(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[2]][[i]]*hosp.cost.mat3),lapply(1:tab.diff, FUN = function(i) sum.outcomes2[[1]][[i]]*gp.cost.mat3)), rowSums)
  net.cost14<-lapply(Map('+',lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[2]][[i]]*hosp.cost.mat4),lapply(1:tab.diff, FUN = function(i) sum.outcomes3[[1]][[i]]*gp.cost.mat4)), rowSums)
  
  
  QalyGain12 <- QALYdifferences(Dataset1, Dataset2, strain=strain,discount=dcount)
  QalyGain13 <- QALYdifferences(Dataset1, Dataset3, strain=strain,discount=dcount)
  QalyGain14 <- QALYdifferences(Dataset1, Dataset4, strain=strain,discount=dcount)
  
  
  ###########################################################################
  ## QALY LOSSES CALCULATIONS
  ##########################################################################
  
  QALY.loss.AJ<-read.csv(paste('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/',"DATA/AJ_QALY_list.csv",sep="/"))
  QALY.loss.cases<-sample(QALY.loss.AJ$v,n.samples,replace=TRUE)
  QALY.loss.hosp<-rnorm(n.samples,mean=0.018,sd=0.0018); rownames(QALY.loss.hosp)<-NULL
  
  #QALY DISCOUNT FOR INFLATION 
  if(dcount==3.5) death.QALY.loss<-apportion(c(22.61,23.00,22.32,18.84,12.71,6.52,24.07,24.05,23.18,20.33,14.62,7.28))
  if(dcount==1.5) death.QALY.loss<-apportion(c(37.06,37.61,35.19,27.16,16.09,7.41,40.77,40.62,37.84,30.40,19.10,8.40))
  if(dcount==0) death.QALY.loss<-apportion(c(60.93,61.62,54.96,38.23,19.75,8.24,69.96,69.45,61.64,44.47,24.14,9.46))
  
  #order return(list(QALY.loss1,QALY.loss2,QALY.cases1,QALY.cases2,QALY.hosp1,QALY.hosp2,QALY.death1,QALY.death2))
  #where label '2' denotes experimental observation and '1' the status quo
  QALY.loss<-CONVERT_INCIDENCE_TO_QALY(Dataset1,Dataset2,strain=strain,dcount)
  QALY.l3<-CONVERT_INCIDENCE_TO_QALY(Dataset1,Dataset3,strain=strain,dcount)
  QALY.l4<-CONVERT_INCIDENCE_TO_QALY(Dataset1,Dataset4,strain=strain,dcount)
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #QALY LOSS TOTAL
  QALY.loss1<-QALY.loss[[1]]#status quo
  QALY.loss2<-QALY.loss[[2]] #new program
  QALY.loss3<-QALY.l3[[2]]
  QALY.loss4<-QALY.l4[[2]]
  
  tot.QALY.loss<-list(QALY.loss1, QALY.loss2, QALY.loss3, QALY.loss4)
  save(tot.QALY.loss, file=paste0(fname[i.country],'tot.qaly.loss',strain.name[strain],'all.programs','d', dcount))
  
  sum.cost1<-sum.cost2<-sum.cost3<-sum.cost4<-list()
  #for net benefit, NOT MARGINAL
  
    sum.cost1 <- Map('+',vaccine.cost.statquo[1:14], hc.cost1) #healthcare costs + vaccine costs
    sum.cost2 <- Map('+',vaccine.cost2[1:14], hc.cost2)
    sum.cost3 <- Map('+',vaccine.cost3[1:14], hc.cost3)
    sum.cost4 <- Map('+',vaccine.cost4[1:14], hc.cost4)
  
  tot.costs.loss<-list(sum.cost1, sum.cost2, sum.cost3, sum.cost4)
  
  save(tot.costs.loss, file=paste0(fname[i.country],'tot.cost.loss',strain.name[strain],'all.programs','d', dcount))

  
  #QALY LOSS DUE TO INFLUENZA INFECTION
  
  QALY.cases1<-(Reduce('+', QALY.loss[[3]])/tab.diff+Reduce('+', QALY.l3[[3]])/tab.diff+Reduce('+', QALY.l4[[3]])/tab.diff)/3
  QALY.cases2<-Reduce('+', QALY.loss[[4]])/tab.diff
  QALY.cases3<-Reduce('+', QALY.l3[[4]])/tab.diff
  QALY.cases4<-Reduce('+', QALY.l4[[4]])/tab.diff
  
  #QALY LOSS DUE TO HOSPITALIZATION
  QALY.hosp1<-(Reduce('+', QALY.loss[[5]])/tab.diff+Reduce('+', QALY.l3[[5]])/tab.diff+Reduce('+', QALY.l4[[5]])/tab.diff)/3
  QALY.hosp2<-Reduce('+', QALY.loss[[6]])/tab.diff
  QALY.hosp3<-Reduce('+', QALY.l3[[6]])/tab.diff
  QALY.hosp4<-Reduce('+', QALY.l4[[6]])/tab.diff
  
  #QALY LOSS DUE TO Death
  QALY.death1<-(Reduce('+', QALY.loss[[7]])/tab.diff+Reduce('+', QALY.l3[[7]])/tab.diff+Reduce('+', QALY.l4[[7]])/tab.diff)/3
  QALY.death2<-Reduce('+', QALY.loss[[8]])/tab.diff
  QALY.death3<-Reduce('+', QALY.l3[[8]])/tab.diff
  QALY.death4<-Reduce('+', QALY.l4[[8]])/tab.diff
  
  #QALY loss total averaged
  QALY.tot1<-(Reduce('+', QALY.loss[[1]])/tab.diff+Reduce('+', QALY.l3[[1]])/tab.diff+Reduce('+', QALY.l4[[1]])/tab.diff)/3
  QALY.tot2<-Reduce('+', QALY.loss[[2]])/tab.diff
  QALY.tot3<-Reduce('+', QALY.l3[[2]])/tab.diff
  QALY.tot4<-Reduce('+', QALY.l4[[2]])/tab.diff
  
  #QALY LOSS DUE TO HOSPITALIZATION and DISEASE
  
  #where dataset1 is always the status quo
  
  net.vaccine.cost2<-lapply(Map('-',vaccine.cost2[1:14], vaccine.cost.statquo[1:14]), rowSums)
  net.vaccine.cost3<-lapply(Map('-',vaccine.cost3[1:14], vaccine.cost.statquo[1:14]), rowSums)
  net.vaccine.cost4<-lapply(Map('-',vaccine.cost4[1:14], vaccine.cost.statquo[1:14]), rowSums)
  
  net.cost12<-lapply(net.cost12, unname)
  net.cost13<-lapply(net.cost13, unname)
  net.cost14<-lapply(net.cost14, unname)
  
  
  net.cost122 <- Map('-',net.vaccine.cost2[1:14],net.cost12)
  net.cost132 <- Map('-',net.vaccine.cost3[1:14],net.cost13)
  net.cost142 <- Map('-',net.vaccine.cost4[1:14],net.cost14)
  
  net.qaly.gain12<-lapply(QalyGain12, rowSums)
  net.qaly.gain13<-lapply(QalyGain13, rowSums)
  net.qaly.gain14<-lapply(QalyGain14, rowSums)
  
  ICER12<-unlist(mapply('/',net.cost122, net.qaly.gain12, SIMPLIFY=FALSE))
  ICER13<-unlist(mapply('/',net.cost132, net.qaly.gain13, SIMPLIFY=FALSE))
  ICER14<-unlist(mapply('/',net.cost142, net.qaly.gain14, SIMPLIFY=FALSE))
  
  
  ICER.set2<-list(net.qaly.gain12,net.cost122,ICER12) #final set
  ICER.set3<-list(net.qaly.gain13,net.cost132,ICER13) #final set
  ICER.set4<-list(net.qaly.gain14,net.cost142,ICER14) #final set
  
  names(ICER.set2)<-c('net QALY','net cost', 'icer')
  names(ICER.set4)<-names(ICER.set3)<-names(ICER.set2)
  
  inv.names<-c('Status Quo','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
        save(ICER.set2, file=paste0(fname[i.country],'ICER.intl',strain.name[strain],'p', program1, program2,'d', dcount))
        save(ICER.set3, file=paste0(fname[i.country],'ICER.intl',strain.name[strain],'p', program1, program3,'d', dcount)) 
        save(ICER.set4, file=paste0(fname[i.country],'ICER.intl',strain.name[strain],'p', program1, program4,'d', dcount)) 
 
  
  ############################
  
  
  reorg<-function(d1,d2,d3,d4, i.country)
  {
    
    inv1<-rowSums(d1)
    inv2<-rowSums(d2)
    inv3<-rowSums(d3)
    inv4<-rowSums(d4)
    
    rem<-data.frame(rbind(
    unname(cbind(inv1,rep(fname[i.country], denom), rep(inv.names2[program1], denom))),
    unname(cbind(inv2,rep(fname[i.country], denom), rep(inv.names2[program2], denom))),
    unname(cbind(inv3,rep(fname[i.country], denom), rep(inv.names2[program3], denom))),
    unname(cbind(inv4,rep(fname[i.country], denom), rep(inv.names2[program4], denom)))
    ))
    
    return(rem)
  }
  
  denom<-length(rowSums(Dataset1[[1]]))
  
  GP.costs.l<-reorg(gp.cost.complete1, gp.cost.complete2, gp.cost.complete3, gp.cost.complete4, i.country)
  Hosp.costs.l<-reorg(hosp.cost.complete1, hosp.cost.complete2, hosp.cost.complete3, hosp.cost.complete4, i.country)
  sum.costs.l<-reorg(Reduce('+',sum.cost1)/tab.diff, Reduce('+',sum.cost2)/tab.diff, Reduce('+',sum.cost3)/tab.diff, Reduce('+',sum.cost4)/tab.diff, i.country)
  aged.cases.l<-reorg(aged.cases1, aged.cases2, aged.cases3, aged.cases4, i.country)
  QALY.tot.l<-reorg(QALY.tot1, QALY.tot2, QALY.tot3, QALY.tot4, i.country)
  QALY.death.l<-reorg(QALY.death1, QALY.death2, QALY.death3, QALY.death4, i.country)
  cumi.l<-reorg(cumi1, cumi2, cumi3, cumi4, i.country)
  
  
 add.dose<-data.frame(rbind(
   cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset2)[[i]])/sum(vaccines.doses.statquo[i,]-vaccines.doses.new2[i,])))/tab.diff, rep(fname[i.country], denom), rep(inv.names2[program2], denom)),
  cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset3)[[i]])/sum(vaccines.doses.statquo[i,]-vaccines.doses.new3[i,])))/tab.diff, rep(fname[i.country], denom), rep(inv.names2[program3], denom)),
 cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset4)[[i]])/sum(vaccines.doses.statquo[i,]-vaccines.doses.new4[i,])))/tab.diff, rep(fname[i.country], denom), rep(inv.names2[program4], denom))
 ))
  
  
 avt.dose<-data.frame(rbind(cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset2)[[i]])/sum(vaccines.doses.new2[i,])))/tab.diff,rep(fname[i.country], denom), rep(inv.names2[program2], denom)),
 cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset3)[[i]])/sum(vaccines.doses.new3[i,])))/tab.diff,rep(fname[i.country], denom), rep(inv.names2[program3], denom)),
 cbind(Reduce('+',lapply(1:tab.diff, function(i)  rowSums(Map('-',Dataset1, Dataset4)[[i]])/sum(vaccines.doses.new4[i,])))/tab.diff, rep(fname[i.country], denom), rep(inv.names2[program4], denom))))
  
 avt.dose[[i.country]]<<-avt.dose;
 add.dose[[i.country]]<<-add.dose;
 aged.cases.l[[i.country]]<<-aged.cases.l;
 cumi.l[[i.country]]<<-cumi.l;
 QALY.death.l[[i.country]]<<-QALY.death.l;
 QALY.tot.l[[i.country]]<<-QALY.tot.l;
 GP.costs.l[[i.country]]<<-GP.costs.l;
 Hosp.costs.l[[i.country]]<<-Hosp.costs.l;
 sum.costs.l[[i.country]]<<-sum.costs.l;
 
 csumn<-csums<-c()
  
  
  csums<-cbind(rowSums(aged.cases1), rowSums(aged.cases2), rowSums(aged.cases3), rowSums(aged.cases4),
             rowSums(aged.death1), rowSums(aged.death2), rowSums(aged.death3), rowSums(aged.death4),
             rowSums(aged.GP1), rowSums(aged.GP2), rowSums(aged.GP3), rowSums(aged.GP4),
             rowSums(aged.hosp1), rowSums(aged.hosp2), rowSums(aged.hosp3), rowSums(aged.hosp4))

  #Table of counts

csumn<-matrix(c(paste0(fname[i.country]),paste0(strain.name[strain]), as.numeric(colMeans(csums))), nrow=1)
colnames(csumn)<-c('country', 'strain', 'cases1', 'cases2', 'cases3','cases4',
                   'deaths1', 'deaths2', 'deaths3','deaths4', 
                   'GP consults1', 'GP consults2', 'GP consults3','GP consults4',
                   'hospitalized1', 'hospitalized2', 'hospitalized3', 'hospitalized4')

ifelse(i.country==1,avg.out<<-csumn, avg.out<<-rbind(avg.out, csumn))


#########################

finaltab<-c()
#table of costs
for(oo in 1:3)
{
  ifelse(oo==1,finaltab<-cbind(fname[i.country],strain.name[strain], 
                               estfunc(rowSums(gp.cost.complete1))[oo], estfunc(rowSums(gp.cost.complete2))[oo], estfunc(rowSums(gp.cost.complete3))[oo],estfunc(rowSums(gp.cost.complete4))[oo],
                               estfunc(rowSums(hosp.cost.complete1))[oo], estfunc(rowSums(hosp.cost.complete2))[oo], estfunc(rowSums(hosp.cost.complete3))[oo], estfunc(rowSums(hosp.cost.complete4))[oo],
                               estfunc(rowSums(QALY.tot1))[oo], estfunc(rowSums(QALY.tot2))[oo], estfunc(rowSums(QALY.tot3))[oo], estfunc(rowSums(QALY.tot4))[oo],
                               estfunc(rowSums(QALY.cases1))[oo], estfunc(rowSums(QALY.cases2))[oo], estfunc(rowSums(QALY.cases3))[oo], estfunc(rowSums(QALY.cases4))[oo],
                               estfunc(rowSums(QALY.hosp1))[oo], estfunc(rowSums(QALY.hosp2))[oo], estfunc(rowSums(QALY.hosp3))[oo], estfunc(rowSums(QALY.hosp4))[oo],
                               estfunc(rowSums(QALY.death1))[oo], estfunc(rowSums(QALY.death2))[oo], estfunc(rowSums(QALY.death3))[oo], estfunc(rowSums(QALY.death4))[oo]),
                                                                                                                                                                                                                                                                  finaltab<-t(unname(cbind( 
                                estfunc(rowSums(gp.cost.complete1))[oo], estfunc(rowSums(gp.cost.complete2))[oo], estfunc(rowSums(gp.cost.complete3))[oo],estfunc(rowSums(gp.cost.complete4))[oo],
                                estfunc(rowSums(hosp.cost.complete1))[oo], estfunc(rowSums(hosp.cost.complete2))[oo], estfunc(rowSums(hosp.cost.complete3))[oo], estfunc(rowSums(hosp.cost.complete4))[oo],
                                estfunc(rowSums(QALY.tot1))[oo], estfunc(rowSums(QALY.tot2))[oo], estfunc(rowSums(QALY.tot3))[oo], estfunc(rowSums(QALY.tot4))[oo],
                                estfunc(rowSums(QALY.cases1))[oo], estfunc(rowSums(QALY.cases2))[oo], estfunc(rowSums(QALY.cases3))[oo], estfunc(rowSums(QALY.cases4))[oo],
                                estfunc(rowSums(QALY.hosp1))[oo], estfunc(rowSums(QALY.hosp2))[oo], estfunc(rowSums(QALY.hosp3))[oo], estfunc(rowSums(QALY.hosp4))[oo],
                                 estfunc(rowSums(QALY.death1))[oo], estfunc(rowSums(QALY.death2))[oo], estfunc(rowSums(QALY.death3))[oo], estfunc(rowSums(QALY.death4))[oo]))))

  if(oo==1){
    colnames(finaltab)<-c('country','strain',
                          'GP costs1', 'GP costs2', 'GP costs3', 'GP costs4', 'hospital costs1',
                          'hospital costs2','hospital costs3','hospital costs4',
                          'QALYs lost1', 'QALYs lost2','QALYs lost3','QALYs lost4',
                          'QALYs lost case1', 'QALYs lost case2','QALYs lost case3','QALYs lost case4',
                          'QALYs lost hosp1', 'QALYs lost hosp2','QALYs lost hosp3','QALYs lost hosp4',
                          'QALYs lost death1', 'QALYs lost death2','QALYs lost death3','QALYs lost death4')
    coststab<-melt(finaltab)
  }else{coststab<-cbind(coststab,finaltab[,1])}
}

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
colnames(coststab)<-c('country','strain','variable','value','CIL','CIU')

#coststab<<-coststab
save(coststab, file=paste0(fname[i.country],'CEAcosts',strain.name[strain],'p', program1, program2, program3, program4,'d', dcount)) 
}












