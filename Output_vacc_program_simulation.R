
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



infection_odin_R0 <- function(population, initial_infected, new.vcalendar, contact_matrix, susceptibility, transmissibility, interval, tstep, duration) 
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
                                         rep(new.vcalendar$efficacy[14],2))
                               ,
                               dates = dates,
                               calendar = calendar1[,1:no_groups],
                               gamma1 = 2/0.8, gamma2 = 2/1.8)
  
  y1 <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  
  S<-y1[,2:23]
  Infec.unvacc<-y1[, (ncol(y1) - no_groups + 1):ncol(y1)]
  
  y<-y1[, ((ncol(y1) - no_groups + 1)-22):(ncol(y1) - no_groups)]
  # Returning the differences in cumulative infections from one week to the other
  y2 <- data.frame(y[2:(nrow(y)), ] - y[1:(nrow(y) - 1), ])
  Infec.unvacc2<-data.frame(Infec.unvacc[2:nrow(Infec.unvacc),]-Infec.unvacc[1:(nrow(Infec.unvacc)-1),])
  
  s.prop<-c(1-risk.ratios.ce[1,],risk.ratios.ce[1,])
  
  # matrix for high and low risk groups
  beta.ll<-(new_cij[1:11,1:11]*S[tstep,1:11])*transmissibility
  beta.hl<-(new_cij[12:22,1:11]*S[tstep,1:11])*transmissibility
  beta.lh<-(new_cij[1:11,12:22]*S[tstep,12:22])*transmissibility
  beta.hh<-(new_cij[12:22,12:22]*S[tstep,12:22])*transmissibility
  raw.beta.mat<-rbind(cbind(new_cij[1:11,1:11], new_cij[1:11,12:22]), cbind(new_cij[12:22,1:11], new_cij[12:22,12:22]))
  beta.mat<-rbind(cbind(beta.ll, beta.lh),cbind(beta.hl, beta.hh))
  
  #beta.overall<-new_cij[1:11,1:11]*(susc[tstep,1:11]+susc[tstep,12:22])*transmissibility
  
  #beta.mat[beta.mat < 0] <- 0
  lam<-as.numeric(eigen(beta.mat)$vectors[,which.max(eigen(beta.mat)$values)])
  Ihl<-lam/(sum(lam))
  
  R0<-sum(((beta.ll*s.prop[1:11]+beta.hl*s.prop[12:22])/duration)*Ihl[1:11]+((beta.hh*s.prop[12:22]+beta.lh*s.prop[1:11])/duration)*Ihl[12:22])
  
  
  ###-------------
  assortcoeff <- function(m) {
    tr <- sum(diag(m))
    sumsq <- sum (rowSums(m)*colSums(m))
    (tr - sumsq) / (1 - sumsq)
  }
  
  assort<-assortcoeff(beta.mat)
  wrap<-c(tstep, R0, assort)
  return(wrap)
}
   

##################
vacc.program.simulation<-function(program, i.country, strain.choice, Reff, version)
{
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));

    #if(i.country==4)
     # {ct.samp<-list.files(pattern=glob2rx(paste0('UKPostContactSample',version, strain.name[strain.choice])))
      #flu.samp<-list.files(pattern=glob2rx(paste0('UKPostSample',version, strain.name[strain.choice])))
      #}else{
      ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PostContactSample',version, strain.name[strain.choice])))
    flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain.choice])))
    #}
    
    if(length(ct.samp)==0 | length(flu.samp)==0)
    {
      stop(print(paste(fname[i.country], 'base version used')))
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
      duration.i<-2/1.8;
      #ODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODE#ODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODE 
      #ODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODE#ODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODEDEODEODEODEODEODEODEODEODEODEODEODEODEODEODEODE 
      
      R0<-infection_odin_R0(population=popv, initial_infected = initial.infected, contact_matrix = poly.data, susceptibility = s.class, transmissibility = fits[5], new.vcalendar = vcalendar2, interval=7, duration=duration.i,tstep=tstep2)
      
      
      ifelse(dd==1, total.R0<<-R0[1:2], total.R0<<-rbind(total.R0,R0[1:2]))
      ifelse(dd==1, total.assort<<-R0[3], total.assort<<-c(total.assort,R0[3]))
      ifelse(dd==1, total.transmit<<-fits[5], total.transmit<<-c(total.transmit,fits[5]))
    }
    
    cov.eff.in<-cov.eff.data[[strain.choice]]
    setDT(ili.counts$total.monitored, keep.rownames = TRUE)
    total.monitored<-ili.counts$total.monitored
    ct.function<<-function(x) {contact_matrix(as.matrix(polymod[x,]),age_sizes[,1], age.group.limits)}
    
    polymod<<-NULL
    polymod<-as.data.frame(fread(file.path('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata',paste(allpolymod[i.country]))))
    polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
    
    ####Start with clean variables so clear all output variables names
    
    total.keep<<-vector('list',14)
    total.akeep<<-vector('list',14)
    total.tr<<-vector('list',14)
    total.R0<<-c()
    total.assort<<-c()
    total.samp<-num.samp
   
    
    year.length<-length(total.keep)
    
    for(i.season in 1:year.length)
    {
      
      #load GB mcmcbatch, and other country ct.table
      setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
      #ifelse(Reff==T,tstep2<-round(runif(1,10,30)), tstep2<-round(runif(1,1,4)))
      ifelse(Reff==T,tstep2<-12, tstep2<-2)
      post.sample.i<-post.sample[[i.season]] 
      colnames(post.sample.i) <- c("eps1", "eps2", "eps3", "psi", "q",
                                   "susc1", "susc2", "susc3", "I0")
      rand.contact.ids<-rand.contact.samp[[i.season]]
      
      #----------------Coverage
      date.labels<-format(as.Date(cov.eff.in[[i.season]]$V1, origin="1970-01-01"), "%Y") 
      
      #----------------Vaccine Strategy
      #set new coverage for the coverage program with new.cov
      log.cov<-predict(new.cov, newdata=data.frame(cov.eff.in[[i.season]]$V1+(15584-cov.eff.in[[i.season]]$V1[1])))
                    
      vcalendar2<-vstrategy2(risk.ratios.ce,program,cov.eff.in,i.season, log.cov)
      
      #we want to save as a list of 14 seasons with 1000 R0 samples each season
      sapply(1:dim(post.sample.i)[1],FUN=thousand.samp) ##initialize 1000sample function
      total.keep[[i.season]]<-total.R0; 
      total.akeep[[i.season]]<-total.assort; 
      total.tr[[i.season]]<-total.transmit
      
    } 
    #end season loop
  
  p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving #intervention names for saving
  ###Save output to file
  
  if(tstep2<5)
  {
    kname<-paste0(fname[i.country],'R0',p.name[program],strain.name[strain.choice])  
    save(total.keep,file=kname) 
    sname<-paste0(fname[i.country],'Assort',p.name[program],strain.name[strain.choice])  
    save(total.akeep,file=sname)
    jname<-paste0(fname[i.country],'transmissibility',p.name[program],strain.name[strain.choice])  
    save(total.transmit,file=jname)
  }else{
    kname<-paste0(fname[i.country],'RE',p.name[program],strain.name[strain.choice])  
    save(total.keep,file=kname) 
    sname<-paste0(fname[i.country],'Assort',p.name[program],strain.name[strain.choice])  
    save(total.akeep,file=sname)
    jname<-paste0(fname[i.country],'transmissibility',p.name[program],strain.name[strain.choice])  
    save(total.transmit,file=jname)
  }
} #end program loop








###############################
#EPI CURVE
####################################
#parameters from MCMC
#post.sample



