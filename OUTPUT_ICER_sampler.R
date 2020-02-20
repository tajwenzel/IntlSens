infection_odin_intl <- function(population, initial_infected, new.vcalendar, contact_matrix, susceptibility, transmissibility, interval) 
{
  # Extract the date used from the vaccine calendar
  begin_date <- as.Date(paste0(format(new.vcalendar$dates[1], "%Y"),"-09-01"))
  t <- as.numeric(seq(begin_date, begin_date + 7*52, interval))
  
  no_groups <- 22
  no_risk_groups <- no_groups/nrow(contact_matrix)
  no_age_groups <- no_groups/no_risk_groups
  
  calendar1 <- new.vcalendar$calendar[c(nrow(new.vcalendar$calendar),1:nrow(new.vcalendar$calendar)),]
  dates <- as.numeric(c(t[1], new.vcalendar$dates))
  
  
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover 
  # risk groups
  new_cij <- matrix(rep(0,no_groups*no_groups), nrow = no_groups)
  for (k in 1:no_risk_groups) {
    for (l in 1:no_risk_groups) {
      lk <- (k - 1)*no_age_groups + 1
      ll <- (l - 1)*no_age_groups + 1
      new_cij[lk:(lk + no_age_groups - 1), ll:(ll + no_age_groups - 1)] <- contact_matrix
    }
  }
  # Set the parameter values
  mod <- gen_seeiir_ag_vacc_nw(no_groups = 22, cij = new_cij, trans = transmissibility,
                               pop = population[1:no_groups],
                               I0 = initial_infected[1:no_groups],
                               susc = rep(susceptibility,no_risk_groups),
                               alpha = c(rep(new.vcalendar$efficacy[1],9),
                                         rep(new.vcalendar$efficacy[7],2),
                                         rep(new.vcalendar$efficacy[8],9),
                                         rep(new.vcalendar$efficacy[14],2)),
                               dates = dates,
                               calendar = calendar1[,1:no_groups],
                               gamma1 = 2/0.8, gamma2 = 2/1.8)
  
  output <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  
  susceptibles<-data.frame(output[,2:23])
  if(colnames(susceptibles)[1]!='S[1]') stop
  
  cumulative.cases <- output[, ((ncol(output) - 2*no_groups + 1):(ncol(output) - 1*no_groups + 1))]
  if(colnames(cumulative.cases)[1]!='cumI[1]') stop
  
  cumulative.cases.novacc <- output[, ((ncol(output) - 1*no_groups + 1):(ncol(output)))]
  if(colnames(cumulative.cases.novacc)[1]!='cumIn[1]') stop

  #Returning the differences in cumulative infections from one week to the other
  incident <- data.frame(cumulative.cases[2:(nrow(cumulative.cases)), ] - cumulative.cases[1:(nrow(cumulative.cases) - 1), ])
  incident.novacc <- data.frame(cumulative.cases.novacc[2:(nrow(cumulative.cases.novacc)), ] - cumulative.cases.novacc[1:(nrow(cumulative.cases.novacc) - 1), ])
  
                                          
  #mutate(Time = as.Date(t[1:nrow(incident)], origin = "1970-01-01"), incident)
  #Cleanup and add Time column
  IDeffects<-incident.novacc/susceptibles[2:dim(susceptibles)[1],]

  return(list( mutate(incident, Time = as.Date(t[1:nrow(incident)], origin = "1970-01-01")),
               mutate(incident.novacc, Time = as.Date(t[1:nrow(incident.novacc)], origin = "1970-01-01")),
               mutate(IDeffects, Time = as.Date(t[1:nrow(IDeffects)], origin = "1970-01-01"))
  ))
}

gen_seeiir_ag_vacc_nw <- odin::odin({
  # Number of groups
  no_groups <- user()
  
  # INITIAL CONDITIONS
  # Population size by age/risk group
  pop[] <- user()
  # Initial infection by age/risk group
  I0[] <- user()
  
  # MODEL PARAMETERS
  # Susceptibility
  susc[] <- user()
  
  # Transmissibility
  trans <- user()
  
  # Latent periods
  gamma1 <- user()
  gamma2 <- user()
  
  # Vaccine related variables 
  dates[] <- user()
  calendar[,] <- user()
  
  # efficacy
  alpha[] <- user()
  
  # Contact matrix
  cij[,] <- user()
  
  # Force of infection
  lambda[] <- trans * susc[i] * (sum(sij[i,]))
  
  # Vaccination. The rate is a step function that changes at each date according
  # to the passed calendar
  vI[] <- interpolate(dates, calendar,"constant")
  # Vaccination is given as a fraction vaccination, here we scale it to 
  # a rate
  sumN[] <- if (vI[i]>0) (S[i]+E1[i]+E2[i]+I1[i]+I2[i]+R[i]) else 0
  v[] <- if (sumN[i]>0) vI[i]*pop[i]/sumN[i] else 0
  
  # Transmission matrix
  sij[,] <- cij[i,j] * (I1[j] + I2[j] + I1v[j] + I2v[j])
  
  # Newly infected
  newInf[] <- lambda[i] * S[i]
  newInfv[] <- lambda[i] * Sv[i]
  
  # THE DERIVATIVES OF THE SEEIIR MODEL
  # Derivatives of the not vaccinated group
  deriv(S[]) <- -newInf[i] - v[i] * S[i]
  deriv(E1[]) <- newInf[i] - gamma1 * E1[i] - v[i] * E1[i]
  deriv(E2[]) <- gamma1 * (E1[i] - E2[i]) - v[i] * E2[i]
  deriv(I1[]) <- gamma1 * E2[i]  - gamma2 * I1[i] - v[i] * I1[i]
  deriv(I2[]) <- gamma2 * (I1[i] - I2[i]) - v[i] * I2[i]
  deriv(R[]) <- gamma2 * I2[i] - v[i] * R[i]
  
  # Derivatives vaccination group
  deriv(Sv[]) <- -newInfv[i] + v[i] * (1-alpha[i]) * S[i]
  deriv(E1v[]) <- newInfv[i] - gamma1 * E1v[i] + v[i] * E1[i]
  deriv(E2v[]) <- gamma1 * (E1v[i] - E2v[i]) + v[i] * E2[i]
  deriv(I1v[]) <- gamma1 * E2v[i]  - gamma2 * I1v[i] + v[i] * I1[i]
  deriv(I2v[]) <- gamma2 * (I1v[i] - I2v[i]) + v[i] * I2[i]
  deriv(Rv[]) <- gamma2 * I2v[i] + v[i] * (R[i] + alpha[i] * S[i])
  
  # Tracking the cumulative amount of infections over time for output of incidence
  deriv(cumI[]) <- newInf[i] + newInfv[i]
  deriv(cumIn[])<-newInf[i]
  
  # Initial value of the variables
  initial(S[1:no_groups]) <- pop[i] - I0[i]
  initial(E1[1:no_groups]) <- 0
  initial(E2[1:no_groups]) <- 0
  initial(I1[1:no_groups]) <- I0[i]
  initial(I2[1:no_groups]) <- 0
  initial(R[1:no_groups]) <- 0
  
  initial(Sv[1:no_groups]) <- 0
  initial(E1v[1:no_groups]) <- 0
  initial(E2v[1:no_groups]) <- 0
  initial(I1v[1:no_groups]) <- 0
  initial(I2v[1:no_groups]) <- 0
  initial(Rv[1:no_groups]) <- 0
  initial(cumI[1:no_groups]) <- 0
  initial(cumIn[1:no_groups])<-0
  
  # Set dimension of all variables/parameters
  dim(dates) <- user()
  dim(calendar) <- user()
  
  dim(pop) <- no_groups
  dim(I0) <- no_groups
  dim(susc) <- no_groups
  dim(lambda) <- no_groups
  dim(v) <- no_groups
  dim(vI) <- no_groups
  dim(sumN) <- no_groups  
  dim(alpha) <- no_groups
  dim(cij) <- c(no_groups, no_groups)
  dim(sij) <- c(no_groups, no_groups)
  
  dim(S) <- no_groups
  dim(E1) <- no_groups
  dim(E2) <- no_groups
  dim(I1) <- no_groups
  dim(I2) <- no_groups
  dim(R) <- no_groups
  dim(Sv) <- no_groups
  dim(E1v) <- no_groups
  dim(E2v) <- no_groups
  dim(I1v) <- no_groups
  dim(I2v) <- no_groups
  dim(Rv) <- no_groups
  dim(cumI) <- no_groups
  dim(newInf) <- no_groups
  dim(newInfv) <- no_groups
  dim(cumIn) <-no_groups
}, verbose = F)



vir.structure<-function(keep) # Function to restructure model output into data age groups
{
  if(is.vector(keep)==TRUE) {keep<-t(as.matrix(keep))}; #must be matrix format to work
  #Convert age groups and risk groups
  #age.group.limits<<-c(1,2,5,12,15,16,25,45,65,75) #upper limits
  c1 <- sum(keep[,c(1,2,3,12,13,14)])
  c2 <- sum(keep[,c(4,5,15,16)])
  c3 <- sum(keep[,c(6,7,8,17,18,19)])
  c4 <- sum(keep[,c(9,20)])
  c5 <- sum(keep[,c(10,11,21,22)])
  regrouped <- c(c1,c2,c3,c4,c5)
  return(regrouped)
}

hilo.structure<-function(keep) # Function to restructure model output into data age groups
{
  if(is.vector(keep)==TRUE) {keep<-t(as.matrix(keep))}; #must be matrix format to work
  #Convert age groups and risk groups
  #age.group.limits<<-c(1,2,5,12,15,17,25,45,65,75) #upper limits
  c1 <- sum(keep[,c(1)])
  c2 <- sum(keep[,c(2,3)])
  c3 <- sum(keep[,c(4,5)])
  c4 <- sum(keep[,c(6,7,8)])
  c5 <- sum(keep[,c(9)])
  c6 <- sum(keep[,c(10,11)])
  c7 <- sum(keep[,c(12)])
  c8 <- sum(keep[,c(13,14)])
  c9 <- sum(keep[,c(15,16)])
  c10 <- sum(keep[,c(17,18,19)])
  c11 <- sum(keep[,c(20)])
  c12 <-sum(keep[,c(21,22)])
  regrouped <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
  return(regrouped)
  
}





singleSENS.samp<-function(strain, i.country, program, num.samp, new.cov, cov.name, version)
{
  #Reconstruct contract matricies from posterior contact data \\
  
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
    
    ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PostContactSample',version, strain.name[strain])))
    flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain])))
    if(length(ct.samp)==0 | length(flu.samp)==0)
    {
    stop(paste(fname[i.country], 'Post sample or PostContactSample not found'))
     }
    load(ct.samp);
    load(flu.samp); 
   
    thousand.samp<-function(dd) 
    { 
      fits<-post.sample.i[dd,] #posterior sample row
      initial.risk<-(rep(10^fits[9], length(age.group.sizes)));
      poly.data<-ct.function(rand.contact.ids[[dd]])
      initial.infected <- stratify_by_risk(initial.risk, risk.ratios.ce) #increased risk group is second
      s.class<-c(fits[6], #0-<1
                 fits[6], #1-2
                 fits[6], #2-4
                 fits[6], #5-11
                 fits[6], #12-14
                 fits[7], #15-<16
                 fits[7], #16-24
                 fits[7], #25-44 
                 fits[7], #45-64
                 fits[8], #65-74
                 fits[8]); #75+) 
      #susceptibility length must match number of age groups
      #duration.i<-1.8;
      #ODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODE#ODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODE 
      
      cases<-infection_odin_intl(population=popv, initial_infected = initial.infected, new.vcalendar = vcalendar2, contact_matrix = poly.data, susceptibility = s.class, transmissibility = fits[5], interval=7)
      #Cases contains whol epi curve
      annual.incidence<-as.vector(colSums(cases[[1]][,1:22]))
      ifelse(dd==1, total.cases<<-unname(annual.incidence), total.cases<<-rbind(total.cases,unname(annual.incidence)))
      unvacc.incidence<-as.vector(colSums(cases[[2]][,1:22]))
      ifelse(dd==1, unvacc.cases<<-unname(unvacc.incidence), unvacc.cases<<-rbind(unvacc.cases,unname(unvacc.incidence)))
      ID.measure<-as.vector(colSums(cases[[3]][,1:22]))
      ifelse(dd==1, ID.calc<<-unname(ID.measure), ID.calc<<-rbind(ID.calc,unname(ID.measure)))
    }
    
  ###------------------Calculate incidence-------------------------
    
  cov.eff.in<-cov.eff.data[[strain]]
  setDT(ili.counts$total.monitored, keep.rownames = TRUE)
  total.monitored<-ili.counts$total.monitored
  ct.function<<-function(x) {contact_matrix(as.matrix(polymod[x,]),age_sizes[,1], age.group.limits)}
  
  polymod<<-NULL
  polymod<-as.data.frame(fread(file.path('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata',paste(allpolymod[i.country]))))
  polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
  
  ####Start with clean variables so clear all output variables names
  total.samp<-num.samp
  total.keep<<-vector('list',14)
  total.incidence<<-c()
  obsv.keep<<-vector('list',14)
  obsv.incidence<<-vector('list',total.samp)
  unvacc.keep<<-vector('list',14)
  unvacc.incidence<<-c()
  
  year.length<-length(total.keep)
  for(i.season in 1:year.length)
  {
    
    #load GB mcmcbatch, and other country ct.table
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
    
    post.sample.i<-post.sample[[i.season]] 
    colnames(post.sample.i) <- c("eps1", "eps2", "eps3", "psi", "q",
                                 "susc1", "susc2", "susc3", "I0")
    rand.contact.ids<-rand.contact.samp[[i.season]]
    
    #----------------Coverage
    date.labels<-format(as.Date(cov.eff.in[[i.season]]$V1, origin="1970-01-01"), "%Y") 
    
    #----------------Vaccine Strategy
    #set new coverage for the coverage program with new.cov
    log.cov<-predict(new.cov, newdata=data.frame(cov.eff.in[[i.season]]$V1+(15584-cov.eff.in[[i.season]]$V1[1])))
    #----------------Vaccine Strategy                                       
    #set new coverage for the coverage program with new.cov                 
    
    #----------------Vaccine Strategy                                           
    #set new coverage for the coverage program with new.cov                     
    vcalendar2<-vstrategy2(risk.ratios.ce,program,cov.eff.in,i.season, log.cov)
    #create vaccine    
    
    #vcalendar1<-vstrategy(risk.ratios.ce,program,cov.eff.in,i.season, new.cov) #create vaccine    
    
    
    #we want to save as a list of 14 seasons with 1000 R0 samples each season
    sapply(1:dim(post.sample.i)[1],FUN=thousand.samp) ##initialize 1000sample function
    total.keep[[i.season]]<-unname(total.cases); 
    unvacc.keep[[i.season]]<-unname(unvacc.cases); 
   
  } 
  #end season loop
  
  p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving
  
  cov.name<-c(0.3, 0.55, 0.7)
  ###Save output to file
  kname<-paste0(fname[i.country],'SensStrat',p.name[program],cov.name,strain.name[strain], version)  
  save(total.keep,file=kname) 
  jname<-paste0(fname[i.country],'unvacc.incidenceIR',p.name[program],cov.name,strain.name[strain], version)  
  save(unvacc.keep,file=jname) 
  lname<-paste0(fname[i.country],'IDeffects',p.name[program],cov.name,strain.name[strain], version)  
  save(unvacc.keep,file=lname) 
  
}


#singleSENS.samp(strain=strainpull, i.country = country, program=1,num.samp = n.samples, new.cov=mL55,cov.name=2,version='v1');

