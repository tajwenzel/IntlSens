
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

vir.structure.rows<-function(keep) # Function to restructure model output into data age groups
{
  if(is.vector(keep)==TRUE) {keep<-t(as.matrix(keep))}; #must be matrix format to work
  #Convert age groups and risk groups
  #age.group.limits<<-c(1,2,5,12,15,17,25,45,65,75) #upper limits
  c1 <- rowSums(keep[,c(1,2,3,12,13,14)])
  c2 <- rowSums(keep[,c(4,5,15,16)])
  c3 <- rowSums(keep[,c(6,7,8,17,18,19)])
  c4 <- rowSums(keep[,c(9,20)])
  c5 <- rowSums(keep[,c(10,11,21,22)])
  regrouped <- cbind(c1,c2,c3,c4,c5)
  return(regrouped)
}


hilo.structure<-function(keep) # Function to restructure model output into data age groups
{
  if(is.vector(keep)==TRUE) {keep<-t(as.matrix(keep))}; #must be matrix format to work
  #Convert age groups and risk groups
  #age.group.limits<<-c(1,2,5,12,15,16,25,45,65,75) #upper limits
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

library('odin')

gen_seeiir_ag_vacc <- odin::odin({
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
  vI[] <- interpolate(dates, calendar, "constant")
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
  initial(cumI[1:no_groups]) <- 0
  initial(cumIn[1:no_groups])<-0
  
  initial(Sv[1:no_groups]) <- 0
  initial(E1v[1:no_groups]) <- 0
  initial(E2v[1:no_groups]) <- 0
  initial(I1v[1:no_groups]) <- 0
  initial(I2v[1:no_groups]) <- 0
  initial(Rv[1:no_groups]) <- 0
  
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

infection_odin_R0 <- function(population, initial_infected, vaccine_calendar, contact_matrix, susceptibility, transmissibility, interval, tstep) 
{
  # Extract the date used from the vaccine calendar
  begin_date <- as.Date(paste0(format(vaccine_calendar$dates[1], "%Y"),"-09-01"))
  t <- as.numeric(seq(begin_date, begin_date + 7*52, interval))
  
  no_groups <- 22
  no_risk_groups <- no_groups/nrow(contact_matrix)
  no_age_groups <- no_groups/no_risk_groups
  
  calendar <- vcalendar1$calendar[c(nrow(vcalendar1$calendar),1:nrow(vcalendar1$calendar)),]
  dates <- as.numeric(c(t[1], vcalendar1$dates))
  
  
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
  
  
  #odes <- infectionODEs(popv, initial.infected, vaccine_calendar=vcalendar,
  #contact_matrix=poly.data[[dd]], 
  #input as a contact matrix, and MUST match iteration number from mcmc
  #susceptibility<-s.class,
  #susceptibility length must match number of age groups
  #transmissibility <- p.transmit,
  #infection_delays=c(2/0.8,2/1.8), interval=7) 
  
  
  duration.i<-1.8;
  # Set the parameter values
  mod <- gen_seeiir_ag_vacc(no_groups = 22, cij = new_cij, trans = transmissibility,
                            pop = population[1:no_groups],
                            I0 = initial_infected[1:no_groups],
                            susc = rep(susceptibility,no_risk_groups),
                            alpha = 
                              c(rep(vaccine_calendar$efficacy[1],9),rep(vaccine_calendar$efficacy[7],2),rep(vaccine_calendar$efficacy[8],9),rep(vaccine_calendar$efficacy[14],2)),
                            dates = dates,
                            calendar = calendar[,1:no_groups],
                            gamma1 = 2/0.8, gamma2 = 2/duration.i)
  
  y <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  
  susc<-y[,2:23]
  #check cases
  cases <- y[, (ncol(y) - no_groups + 1):ncol(y)]
  #Returning the differences in cumulative infections from one week to the other
  I <- data.frame(cases[2:(nrow(cases)), ] - cases[1:(nrow(cases) - 1), ])
  s.prop<-c(1-risk.ratios.ce[1,],risk.ratios.ce[1,])
  
  # matrix
  beta.ll<-(new_cij[1:11,1:11]*susc[tstep,1:11])*transmissibility
  beta.hl<-(new_cij[12:22,1:11]*susc[tstep,1:11])*transmissibility
  beta.lh<-(new_cij[1:11,12:22]*susc[tstep,12:22])*transmissibility
  beta.hh<-(new_cij[12:22,12:22]*susc[tstep,12:22])*transmissibility
  
  #eigen(beta.ll)$vectors[,which.max(eigen(beta.ll)$values)]
  
  beta.mat<-rbind(cbind(beta.ll, beta.lh),cbind(beta.lh, beta.hh))
  lam<-eigen(beta.mat)$vectors[,which.max(eigen(beta.mat)$values)]
  Ihl<-lam/(sum(lam))
  
  R0<-sum(((beta.ll*s.prop[1:11]+beta.hl*s.prop[12:22])/(2/duration.i))*Ihl[1:11]+((beta.hh*s.prop[1:11]+beta.lh*s.prop[12:22])/(2/duration.i))*Ihl[12:22])
  return(R0)
}
