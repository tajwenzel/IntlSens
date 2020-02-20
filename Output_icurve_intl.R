library(fluEvidenceSynthesis)
library(plyr)
library(pander)
library(data.table)
library(triangle)
library(car)
library(proto)
library(ggplot2)
library(MASS)
library(colorRamps)
library(RColorBrewer)
library(tableone)
library(dplyr)
library(tidyr)
library(broom)
library(triangle)
library(mcmc)
library(reshape2)
library(cowplot)
library(ggpubr)
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

rm(list = ls())
#load workspace
#get sampling functions for the International data
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_odin.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_ICER_conversion_22ages.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/SUPPORT_R0_loaders.R')
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
vstrategy<-dget('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_cov_strategy_uk.R')


##load general parameters-------------------
load('Intl_general_parameters.RData')
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/ili.counts.rda')
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/virological.rda')
strain.name<-c('H1N1','H3N2','B')

strainpull<-strain<-strain.choice<-1
country<-4
i.country<-country
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))

#fiftyfive<-coverageB[[15]]$V12
#fiftyfive<-c(rep(0,1),fiftyfive[2:18]+0.012); 
#fiftyfive.short<-fiftyfive[seq(1,dim(coverageB[[15]])[1],4)]
#fiftyfive.cov<-list(fiftyfive.short, fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short)

#new.ili<-list()
#new.ili$ili<-ili.counts$ili[ order(as.numeric(row.names(ili.counts$ili))), ]
#new.ili$total.monitored<-ili.counts$total.monitored[ order(as.numeric(ili.counts$total.monitored$rn)), ]
#####load incidence sampler for international data



incidence.samp.intl<-function(strain.choice, i.country, program, num.samp, new.cov, version)
{
  #Reconstruct contract matricies from posterior contact data \\
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PostContactSample',version, strain.name[strain.choice])))
  flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain.choice])))
  
  if(length(ct.samp)==0 | length(flu.samp)==0)
  {
    stop(paste(fname[i.country], 'base version used'))
  }
  
  load(ct.samp);
  load(flu.samp); 
    
    #tname<-paste0('trace',i.season)
    #post.name<-paste0('post', i.season)
    #assign(tname, ggplot(data=melt(set)) + facet_wrap( ~ Var2, ncol=3, scales="free" ) + geom_line(aes(x=Var1, y=value)), envir = .GlobalEnv)
    #assign(post.name, ggplot(data=melt(set)) + facet_wrap( ~ Var2, ncol=3, scales="free" ) + geom_density(aes(x = value)) + geom_rug(aes(x = value, y = 0), position = position_jitter(height = 0))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10)), envir = .GlobalEnv)
    #----------------Contact matrix
    #generate new contact matrix
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
      
      cases<-infection_odin(population=popv, initial_infected = initial.infected, vaccine_calendar = vcalendar1, contact_matrix = poly.data, susceptibility = s.class, transmissibility = fits[5], interval=7)
      #Cases contains whole epi curve
      total.incidence<-colSums(cases[,1:22])
      ifelse(dd==1, total.cases<<-unname(total.incidence), total.cases<<-rbind(total.cases,unname(total.incidence)));
      #unvacc.incidence<-as.vector(colSums(cases[[2]][,1:22]))
      #ifelse(dd==1, unvacc.cases<<-unname(unvacc.incidence), unvacc.cases<<-rbind(unvacc.cases,unname(unvacc.incidence)))
        #rearrange ode output into 5 age groups given provided in virological/ili data and add dates
        obsv<-vir.structure.rows(cases)
        #cases represents all cases in the population. Cases needs to be re-proportioned into the amount of cases that would be observed for ILI; which is ILI total monitored over populations divided into the the 5 age classes. 
        
        cases.monitored<-obsv[,1:5]*(ltab[(((i.season-1)*52)+1):(i.season*52),])
        
        #for observed cases; defined as number of people in age group who would be recorded for ILI and who would have been positive for strain if tested. Based on Baguelin Supplement, 2013 Figure S6 caption. Samples from m+ = Binomial (Z-number of ili cases in monitored population, epsilon- age-group specific ascertainment rate). This process is represented below.
        obsv.odes<-c()
        obsv.cases<-c()
        
        for (ii in 1:dim(cases)[1]){ 
          #loop through rows to sample based on age group and date. Ascertainment rate is the same through time. 
          obsv.odes<-c(rbinom(1,round(cases.monitored[ii,1]), fits[1]),    #epsilons are always applied on 5 age groups
                       rbinom(1,round(cases.monitored[ii,2]), fits[1]),
                       rbinom(1,round(cases.monitored[ii,3]), fits[2]),
                       rbinom(1,round(cases.monitored[ii,4]), fits[2]),
                       rbinom(1,round(cases.monitored[ii,5]), fits[3]))
          # gets number of odes that would have been observed through ILI monitoring
          ifelse(ii==1, obsv.cases<-cbind(cases$Time[ii],matrix(obsv.odes, nrow=1)), obsv.cases<-rbind(obsv.cases, cbind(cases$Time[ii],matrix(obsv.odes, nrow=1))))
        } 
        
        #row bind all subsequent dates within loop after ii=1
        obsv.incidence[[dd]]<<-obsv.cases;
        
        
        #return(list(total.cases, obsv.incidence))
      }
  
  #########clean variables
  total.samp<-num.samp
  total.keep<<-vector('list',14)
  total.incidence<<-c()
  obsv.keep<<-vector('list',14)
  obsv.incidence<<-vector('list', num.samp)
  unvacc.keep<<-vector('list',14)
  unvacc.incidence<<-c()
  
  #define ILI data 
  short<-stratify_by_age(age_sizes$V1, c(5,15,45,65))
  setDT(ili.counts$total.monitored, keep.rownames = TRUE)
  ll<-apply(ili.counts$total.monitored[,2:6],1, function(x) x/as.data.frame(t(short)))
  ltab<-do.call(rbind.data.frame, ll)
  
  #set coverage
  cov.eff.in<-cov.eff.data[[strain.choice]]
  ct.function<<-function(x) {contact_matrix(as.matrix(polymod[x,]),age_sizes[,1], age.group.limits)}
  
  #set contact rates
  polymod<<-NULL
  polymod<-as.data.frame(fread(file.path('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata',paste(allpolymod[i.country]))))
  polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
  
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
      vcalendar1<<-vstrategy(risk.ratios.ce,program,cov.eff.in,i.season, new.cov) #create vaccine    
      
      #we want to save as a list of 14 seasons with 1000 R0 samples each season
      #sapply(1:dim(post.sample.i)[1],FUN=thousand.samp) ##initialize 1000sample function
      rrg<-sapply(1:num.samp,FUN=thousand.samp) ##initialize 1000sample function
      total.keep[[i.season]]<-unname(total.cases); 
      obsv.keep[[i.season]]<-obsv.incidence;
      
    } 
    #end season loop
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving
  ###Save output to file
  kname<-paste0(i.country,'sample',p.name[program],strain.name[strain.choice])  
  save(total.keep,file=kname) 
  if(program==1)
  {hname<-paste0(i.country,'obsv',p.name[program],strain.name[strain.choice])  
  save(obsv.keep,file=hname)}
  
} #end strain loop



###########Invoke function


n.samples<-100
country<-4
for(country in 1:(length(fname)))
{
  if(country==4) next
  #end vaccine program loop
  incidence.samp.intl(strain=strainpull, i.country = country, program=1,num.samp = n.samples, new.cov=fiftyfive.cov, version=3);
  incidence.samp.intl(strain=strainpull, i.country = country, program=4,num.samp = n.samples, new.cov=fiftyfive.cov, version=3);
  incidence.samp.intl(strain=strainpull, i.country = country, program=2,num.samp = n.samples, new.cov=fiftyfive.cov, version=3);
  incidence.samp.intl(strain=strainpull, i.country = country, program=3,num.samp = n.samples, new.cov=fiftyfive.cov, version=3);
}



########################################################################################
####### Support functions ########
####################################################################################################

i.outcomes.loader<-function(strain, i.country)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
  #Rearrange tables to correct order of interventions------------------
 
  
  load(paste0(i.country,'sample', p.name[1],strain.name[strain])); avg.incidence1<-total.keep;
  load(paste0(i.country,'sample', p.name[2],strain.name[strain])); avg.incidence2<-total.keep;
  load(paste0(i.country,'sample', p.name[3],strain.name[strain])); avg.incidence3<-total.keep;
  load(paste0(i.country,'sample', p.name[4],strain.name[strain])); avg.incidence4<-total.keep;
 
  
load(paste0(i.country,'obsv', p.name[1], strain.name[strain]));  avg.incidence0.1<-obsv.keep;
  #load in obsv cases, and assign new name
  o.name<-paste0('o.table',strain.name[strain]);
  
  assign(o.name, avg.incidence0.1, envir = .GlobalEnv)
  
  var.name<-paste0('table',strain.name[strain])
  
  assign(var.name, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4), envir = .GlobalEnv)
}

hypergeo.conf<-function (gamma, x, N, n)
{
  if(n==0|N==0){return(0)} 
  else {
    alpha <- 1 - gamma
    D0 <- ceiling(N * qbeta(gamma, x + 1, n - x))
    p1 <- phyper(x, D0, N, n)
    p2 <- phyper(x, D0 - 1, N, n)
    
    if(p1 <= alpha & p2 > alpha) {return(D0)}
    if(p1 > alpha) {D0 <- D0 + 1; return(D0)}
    if(p2 <= alpha){D0 <- D0 - 1; return(D0)}
  }
}


######################## MAIN GRAPH FUNCTION/LOOP FOR COMBINED YEARLY GRAPHS------




##########----------
version<-'v1'
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain.choice])))
  
  if(length(flu.samp)==0)
  {
    stop(print(paste('missing posterior sample for', fname[i.country], version, strain.name[strain.choice])))
  }
  
  load(flu.samp); 
  
  for(i.season in 1:14)
  {
    cparms.name<-paste0(fname[i.country],'pars',i.season)
    post.sample.i<-post.sample[[i.season]]
    colnames(post.sample.i) <- c("eps1", "eps2", "eps3", "psi", "q",
                                 "susc1", "susc2", "susc3", "I0")
    
    sample.m<-melt(post.sample.i)
    
    assign(cparms.name, ggplot(sample.m, aes(x=value)) + 
             geom_histogram(aes(y=..density..), colour="black", fill="white")+
             geom_density(alpha=.2, fill="#FF6666") +facet_wrap(~Var2, ncol=2, nrow=5, scales='free')+labs(title=paste0('Posterior Parameters'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10)))
    
    #post.pars <- ddply(sample.m, 'X2', summarise,
     #                 mean1 = mean(value),
      #                median1=median(value),
       #               CRL= quantile(value, 0.025),
        #              CRU= quantile(value, 0.025)
         #             )
    
    #post.pars<-cbind(post.pars, rep(i.season, 9), rep(fname[i.country], 9))
    
    post.pars<-cbind(sample.m, rep(i.season, dim(sample.m)[1]), rep(fname[i.country], dim(sample.m)[1]))
    
    ifelse(i.country==1 & i.season==1, post.pars.out<-post.pars, post.pars.out<-rbind(post.pars.out, post.pars))
    
  }
}   



#colnames(post.pars.out)<-c('parameter','mean','median','CRL','CRU','season','country')
#post.pars.out$country<-factor(post.pars.out$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
#col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

####susceptibility
colnames(post.pars.out)<-c('num','parameter','value','season','country')
post.pars.out$country<-factor(post.pars.out$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
susc.pars<-post.pars.out[post.pars.out$parameter=='susc1'|post.pars.out$parameter=='susc2'|post.pars.out$parameter=='susc3',]
sus.p<-ggplot(susc.pars, aes(x=country, y=value, fill=country))+geom_violin(scale='count', draw_quantiles = 0.5)+facet_wrap(~parameter, nrow=3)+geom_hline(yintercept=0)+labs(x='Contact Matrix', y='Parameter Value')
ggsave(paste0("suscparameters",strain.name[strain],'allcountries', version, ".png"),device='png',plot=sus.p, width=30,height=18,units='cm')

#####epsilons
eps.pars<-post.pars.out[post.pars.out$parameter=='eps1'|post.pars.out$parameter=='eps2'|post.pars.out$parameter=='eps3',]
eps.p<-ggplot(eps.pars, aes(x=country, y=value, fill=country))+geom_violin(scale='width', draw_quantiles = 0.5)+facet_wrap(~parameter, nrow=3, scales = 'free')+geom_hline(yintercept=0)+labs(x='Contact Matrix', y='Parameter Value')
ggsave(paste0("epsparameters",strain.name[strain],'allcountries',version, ".png"),device='png',plot=eps.p, width=30,height=18,units='cm')


####"I0"    "psi"   "q" 
rem.pars<-post.pars.out[post.pars.out$parameter=='I0'|post.pars.out$parameter=='psi'|post.pars.out$parameter=='q',]
rem.p<-ggplot(rem.pars, aes(x=country, y=value, fill=country))+geom_violin(scale='width', draw_quantiles = 0.5)+facet_wrap(~parameter, nrow=3, scales = 'free')+geom_hline(yintercept=0)+labs(x='Contact Matrix', y='Parameter Value')
ggsave(paste0("remparameters",strain.name[strain],'allcountries',version, ".png"),device='png',plot=rem.p, width=30,height=18,units='cm')



for(i.country in 1:length(fname))
{

i.outcomes.loader(1,i.country)

  work.table<-get(paste0('o.table',strain.name[strain]))
  
  #reduce sampled lists to season list of means and standard error
  ktabB<-stabB<-list()
  for(oo in 1:length(work.table))
  { 
    #get means
    m.pot<-Reduce('+', work.table[[oo]])/length(o.tableH1N1[[1]])
    ktabB[[oo]]<-m.pot
    
    #standard deviation
    m.se<-sqrt(Reduce('+',lapply(lapply(work.table[[oo]], function(x) x-m.pot), function(x) x^2))/length(work.table[[1]]))
    stabB[[oo]]<-m.se[,2:6]
  }
  
  
  #staring graph data.frames
  #date for data
  begin_date <- as.Date("1995-09-01")
  interval<-7
  t <- as.Date(as.numeric(seq(begin_date, begin_date+365*14, interval)), origin='1970-01-01')
  
  #calculate ILI (mij+) | strain 
  names(virological$pos.by.strain)<-c('B','H3N2','H1N1')
  #virological is not ordered 1->3 so need to reassign
  if(strain==1) {vir.strain<-virological$pos.by.strain$H1N1}
  if(strain==2) {vir.strain<-virological$pos.by.strain$H3N2}
  if(strain==3) {vir.strain<-virological$pos.by.strain$B}
  
  #row rearrange
  short<-stratify_by_age(age_sizes$V1, c(5,15,45,65))
  ll<-apply(ili.counts$total.monitored[,2:6],1, function(x) x/as.data.frame(t(short)))
  ltab<-do.call(rbind.data.frame, ll)
  
  #calculate
  rlB<-(vir.strain[,1:5]/virological$no.samples[,1:5])*ili.counts$ili[,1:5]
  rlB[is.na(rlB)] <- 0 #remove NAs
  t<-t[1:dim(rlB)[1]]
  rlb.g<-data.frame(cbind(t,rlB)) #add dates
  
  colnames(rlb.g)<-c('date',"[0,5)",   "[5,15)" , "[15,45)", "[45,65)", "[65,+)" )
  net.set<-melt(rlb.g, id.vars='date')
  
  rlb.ci<-matrix(rep(0,dim(rlb.g)[1]*(dim(rlb.g)[2]-1)), ncol=(dim(rlb.g)[2]-1)) #allocate memory for confidence intervals
  
  t2 <- as.Date(as.numeric(seq(as.Date("1995-09-01"), (begin_date+365*14), interval)), origin='1970-01-01') #edwin starts model at 09/01, although should be week 35 of 1995 which is 08/28
  t2<-t2[1:(14*52)]
  for(j in 1:5)
  {
    for(i in 1:length(t2))
    {
      rlb.ci[i,j]<-hypergeo.conf(gamma = .95, x = as.numeric(vir.strain[i,j]), as.numeric(ili.counts$ili[i,j]), as.numeric(virological$no.samples[i,j]))
    }
  }
  
  #MODEL OUTCOMES----
  ###start model outcomes rerrange
  #call tables and rearrange list
  
  ktab<-do.call(rbind.data.frame, ktabB)
  stab<-do.call(rbind.data.frame, stabB)
  #add standard error to dataframe
  xx<-ktab; colnames(xx)<-c('date',"[0,5)",   "[5,15)" , "[15,45)", "[45,65)", "[65,+)") 
  upper<-cbind(ktab[,1],ktab[,2:6]+1.96*stab); colnames(upper)<-c('date',"[0,5)",   "[5,15)" , "[15,45)", "[45,65)", "[65,+)") 
  lower<-cbind(ktab[,1],ktab[,2:6]-1.96*stab); colnames(lower)<-c('date', "[0,5)",   "[5,15)" , "[15,45)", "[45,65)", "[65,+)"); 
  
  #remove negatives
  upper[upper<0]<-0 
  lower[lower<0]<-0
  
  #melt data frame of model mean outcomes  
  net.set2<-cbind(melt(xx, id.vars=c('date')), melt(upper, id.var=c('date'))$value, melt(lower, id.var=c('date'))$value)
  colnames(net.set2)<-c('date', 'variable','value', 'upper','lower'); 
  net.set2[,1]<-as.Date(net.set2[,1], origin = '1970-01-01')
  #clean data frame and add dates and confidence intervals
  rci<-data.frame(t2,rlb.ci[1:length(t2),])
  colnames(rci)<-c('date', "[0,5)",   "[5,15)" , "[15,45)", "[45,65)", "[65,+)")
  rrci<-melt(rci, id.vars='date')

  
  
  ####loop for graphing age-structured incidence for each season---
  #each is assigned a name 'p' and number of the season
  
  #for(i.season in 1:14){
   # ptemp <- ggplot(data=melt(rlb.g[(((i.season-1)*52)+1):(i.season*52),], id.vars='date'))+facet_wrap( ~variable, ncol=1, scales="free" )+geom_errorbar(data= melt(rci[(((i.season-1)*52)+1):(i.season*52),], id.vars='date'), aes( ymin=0, ymax=value, x=date), width=.1, alpha=0.5)+geom_point(aes(x=date, y=value), size=1)+ylab('Number of positive ILI cases')+theme(legend.position="none")+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10))
    
  #  net.set2l<-cbind(melt(xx[(((i.season-1)*52)+1):(i.season*52),], id.vars=c('date')), melt(upper[(((i.season-1)*52)+1):(i.season*52),], id.var='date')$value, melt(lower[(((i.season-1)*52)+1):(i.season*52),], id.var='date')$value)
    #colnames(net.set2l)<-c('date', 'variable','value', 'upper','lower'); 
    #net.set2l[,1]<-as.Date(net.set2l[,1], origin = '1970-01-01')
    
    
   # p2<-ptemp+geom_line(data=net.set2l, aes(x=date, y=value, color='red'))+geom_ribbon(data=net.set2l,aes(ymax=upper,ymin=lower, x=date, fill='band'), alpha=0.3)+scale_fill_manual("",values="red")+theme(legend.position="none")+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10))
    
    #dname<-paste0('g', i.season)
    #assign(dname, p2)
  #}
  
  
  ###Incidence curve for all years  
  #all years incidence curve
  
  k1<-ggplot(data=net.set)+facet_wrap(~variable, ncol=1, scales="free" )+geom_point(aes(x=as.Date(date), y=as.numeric(value)), size=1)+ geom_errorbar(data= rrci, aes( ymin=0, ymax=as.numeric(value), x=date), width=.1, alpha=0.5)+theme(legend.position="none")
  tt2<-k1+geom_line(data=net.set2, aes(x=date, y=value, colour='red'))+geom_ribbon(data=net.set2,aes(ymax=upper,ymin=lower, x=date, fill='band', colour='red'), alpha=0.5)+labs(title= paste(strain.name[strain], fname[i.country],p.name[program], 'Incidence Curve'))+theme(legend.position="none")
  
  ggsave(paste0("ICurve",strain.name[strain],fname[i.country], ".png"),device='png',plot=tt2, width=43,height=20,units='cm')

}
  ###########PULL DATA FOR CONTACT MATRICIES----------------------
  
  
for(i.country in 1:length(fname))
{

  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
    load(paste0('post.means',fname[i.country], 'H1N1')); post.meansH1N1<-post.means
    load(paste0('prior.means',fname[i.country], 'H1N1')); prior.meansH1N1<-prior.means
  
  for(i.season in 1:14)
  {
    #pull means for each season 
    post.mat.mean<-post.means[[i.season]]*age.group.sizes
    prior.mat.mean<-prior.means[,,i.season]*age.group.sizes
    
    
    #Code for 7 compartment plot
    col.post<-cbind(post.mat.mean[,1],(post.mat.mean[,2]+post.mat.mean[,3]),(post.mat.mean[,4]+post.mat.mean[,5]),(post.mat.mean[,6]+post.mat.mean[,7]), post.mat.mean[,8],post.mat.mean[,9], (post.mat.mean[,10]+post.mat.mean[,11]))
    post.mat.s<-cbind(col.post[1,],(col.post[2,]+col.post[3,]),(col.post[4,]+col.post[5,]),(col.post[6,]+col.post[7,]), col.post[8,],col.post[9,], (col.post[10,]+col.post[11,]))
    
    post.m.ggplot <- melt(post.mat.s)
    #prior
    
    col.prior<-cbind(prior.mat.mean[,1],(prior.mat.mean[,2]+prior.mat.mean[,3]),(prior.mat.mean[,4]+prior.mat.mean[,5]),(prior.mat.mean[,6]+prior.mat.mean[,7]), prior.mat.mean[,8],prior.mat.mean[,9], (prior.mat.mean[,10]+prior.mat.mean[,11]))
    prior.mat.s<-cbind(col.prior[1,],(col.prior[2,]+col.prior[3,]),(col.prior[4,]+col.prior[5,]),(col.prior[6,]+col.prior[7,]), col.prior[8,],col.prior[9,], (col.prior[10,]+col.prior[11,]))
    pre.m.ggplot<- melt(prior.mat.s)
    
    b.range<-round(c(quantile(c(post.m.ggplot$value, pre.m.ggplot$value), c(.001,.125,.25,.375,.5,.725,.8,.9, 1))),2)
    #section data into discrete intervals
    
    combo<-c(post.m.ggplot$value, pre.m.ggplot$value)
    if(any(b.range==0))
       {
    combo2<-combo[combo >= b.range[max(which(b.range==0))+1]]
    c.range<-round(quantile(combo2,c(0.001,.125,.25,.375,.5,.725,.8,.9, 1)),2)
    c.range2<- na.omit(c(0, 1e-5, c.range, round(c.range[11])+1))
    }else{
      c.range2<-b.range;
      c.range2[1]<-0;
      c.range2<-na.omit(c.range2)
    }
    discrete1<-cut(post.m.ggplot$value, breaks=c.range2, right=FALSE, ordered_result = FALSE)
    post.discrete<-cbind(post.m.ggplot, discrete1)
    discrete2<-cut(pre.m.ggplot$value, breaks=c.range2, right=FALSE, ordered_result = FALSE)
    prior.discrete<-cbind(pre.m.ggplot, discrete2)
    #dev.off()
    #dev.off()
    cprior.name<-paste0(fname[i.country], 'cpr',i.season)
    cpost.name<-paste0(fname[i.country],'cpo',i.season)
    
    b.label<-c('0',
               paste0('[0','-',c.range2[3],')'),
               paste0('[',c.range2[3],'-',c.range2[4],')'),
               paste0('[',c.range2[4],'-',c.range2[5],')'),
               paste0('[',c.range2[5],'-',c.range2[6],')'),
               paste0('[',c.range2[6],'-',c.range2[7],')'),
               paste0('[',c.range2[7],'-',c.range2[8],')'),
               paste0('[',c.range2[8],'-',c.range2[9],')'),
               paste0('[',c.range2[10],'-', c.range2[11],')'),
               paste0('[',c.range2[11],'+',')'))
    
    b.label<-b.label[1:length(c.range2)];
    Assort.label<-round(A.cdata$mean2[A.cdata$c.repeat==fname[i.country]&A.cdata$Season==(i.season+1994)],3)
    
    assign(cpost.name, ggplot(data=post.discrete, aes(x=Var1,y=Var2,fill=factor(discrete1)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Posterior Contact Matrix\n (Discrete Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=b.range, discrete=TRUE)+scale_fill_viridis(option='inferno',labels=b.label[is.na(match(sort(seq(1,9,1)), sort(match(unique(post.discrete$discrete1), levels(post.discrete$discrete1)))))==F], name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10))+
      annotate("text", x = 8, y = 8, label = paste0('\u03C1 =', Assort.label)))
    
    assign(cprior.name, ggplot(data=prior.discrete, aes(x=Var1,y=Var2,fill=factor(discrete2)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Prior Contact Matrix\n (Discrete Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=b.range, discrete=TRUE)+scale_fill_viridis(option='inferno',labels=b.label[is.na(match(sort(seq(1,9,1)), sort(match(unique(prior.discrete$discrete2), levels(prior.discrete$discrete2)))))==F], name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=10),axis.title=element_text(size=10)))
                                                                                           
  } 

   # meh<-lapply(1:14, function(yy) {ggarrange(get(paste0('g',yy)), ggarrange(get( paste0(fname[i.country],'cpr',yy)), get(paste0(fname[i.country],'cpo',yy)),  nrow = 2, labels=c(fname[i.country])))})
    #xx<-marrangeGrob(meh,
                     #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                     #labels=c('A'),
     #                nrow = 1, ncol=1)
    #ggsave(paste0("icurve.mtx.parms",strain.name[strain], fname[i.country],".pdf"),device='pdf',plot=xx, width=54,height=24,units='cm')
    
    
    ##########################################
    ###########GRAPH OF ALL OUTPUT#####################------------------  
    ########################################## 
    all3<-lapply(1:14, function(yy) {ggarrange(get(paste0('g',yy)), ggarrange(get(paste0(fname[i.country],'cpr',yy)), get(paste0(fname[i.country],'cpo',yy)),  nrow = 2, labels=c(fname[i.country])), get(paste0(fname[i.country], 'pars', yy)), ncol=3, nrow=1)})
    graphs.allpost<-marrangeGrob(all3,
                                 #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                                 #labels=c('A'),
                                 nrow = 1, ncol=1)
    ggsave(paste0("icurve.mtx.allparms",strain.name[strain], fname[i.country],".pdf"),device='pdf',plot=graphs.allpost, width=60,height=28,units='cm')
    
    
}





##########################################
###########GRAPH OF ALL OUTPUT#####################------------------  
########################################## 
all3<-lapply(1:14, function(yy) {ggarrange(get(paste0('g',yy)), ggarrange(get(paste0(fname[i.country],'cpr',yy)), get(paste0(fname[i.country],'cpo',yy)),  nrow = 2, labels=c(fname[i.country])), get(paste0(fname[i.country], 'pars', yy)), ncol=3, nrow=1)})
graphs.allpost<-marrangeGrob(all3,
                 #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                 #labels=c('A'),
                 nrow = 1, ncol=1)
ggsave(paste0("icurve.mtx.allparms",strain.name[strain], fname[i.country],".pdf"),device='pdf',plot=graphs.allpost, width=60,height=28,units='cm')

    
    
    #####################------------------  
#Check MCMC output via Gelman Statistic
    
library(coda)
mcmc.diagnotic<-function(i.country, strain.choice, mark,pull, version)
{
  
  
  polymod<-fread(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Inputdata', paste(allpolymod[i.country])))
  polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
  
  #Reconstruct contract matricies from posterior contact data \\
  
  #if(i.country==4)
  #{setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
   # flu.tablesB<-list.files(pattern=glob2rx('GBflu.u5*B'));
    #flu.tablesH1<-list.files(pattern=glob2rx('GBflu.u5*H1'));
    #flu.tablesH3<-list.files(pattern=glob2rx('GBflu.u5*H3'));
    
    #ct.tablesB<-list.files(pattern=glob2rx('GBct.u5*B'))
    #ct.tablesH1<-list.files(pattern=glob2rx('GBct.u5*H1'))
    #ct.tablesH3<-list.files(pattern=glob2rx('GBct.u5*H3'))
  #}else{
        setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
        
        if(strain.choice==1){ct.tables<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version, '*H1', sep='')))}
        if(strain.choice==2){ct.tables<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.', version, '*H3', sep='')))}
        if(strain.choice==3){ct.tables<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.', version, '*B', sep='')))}
        
        if(strain.choice==3){flu.tables<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.', version, '*B',sep='')));}
        if(strain.choice==1){flu.tables<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version, '*H1',sep='')));}
        if(strain.choice==2){flu.tables<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version, '*H3',sep='')));}
      #}
  
  if(strain.choice==1){year.length<-length(ct.tables)}
  if(strain.choice==2){year.length<-length(ct.tables)}
  if(strain.choice==3){year.length<-length(ct.tables)}
  
  gelman.output<-list()
  last.pull<<-list()
  for(hh in 1:year.length)
  {
    
    load(flu.tables[[hh]]); 
    chain1<-as.mcmc(mcmc.result$batch[(mark:dim(mcmc.result$batch)[1]),])
    last.pull[[hh]]<<- mcmc.result$batch[(dim(mcmc.result$batch)[1]),]
    #load(flu.tablesB); chain3<- as.mcmc(mcmc.result$batch[(mark:dim(post.sample)[1]),])
    
    if(pull==F)
    {
      load(flu.tables.b[[hh]]); chain2<- as.mcmc(mcmc.result$batch[(mark:dim(mcmc.result$batch)[1]),])
      xx<-mcmc.list(chain1, chain2)
      f<-gelman.diag(xx, confidence = 0.95, transform=T, autoburnin=F, multivariate = T)
      
      gelman.output[[hh]]<- f
      gelman.plot(xx)
    }
    
  }
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata','IntlSens'));
  save(last.pull, file=paste0(fname[i.country],strain.name[strain.choice], 'last.pull'))
}

version<-'v1'
for(i.country in 1:length(fname))
{if(i.country==4) next
  mcmc.diagnotic(i.country,strain.choice=2,50000, pull=T, version='v1')}


mcmc.diagnotic(i.country=1,strain.choice=2,50000, pull=T, version='v1')

chain1.1<-cbind(seq(1,dim(chain1)[1],1), chain1, rep('chain1',dim(chain1)[1]))
colnames(chain1.1) <- c('iter',"eps1", "eps2", "eps3", "psi", "q",
                        "susc1", "susc2", "susc3", "I0", 'chain')
#chain2.1<-cbind(seq(1,dim(chain2)[1],1), chain2, rep('chain2',dim(chain2)[1]))
#colnames(chain2.1) <- c('iter',"eps1", "eps2", "eps3", "psi", "q",
 #                       "susc1", "susc2", "susc3", "I0", 'chain')
#chain3.1<-cbind(seq(1,dim(chain3)[1],1), chain3, rep('chain2',dim(chain3)[1]))
#colnames(chain3.1) <- c('iter',"eps1", "eps2", "eps3", "psi", "q",
#                        "susc1", "susc2", "susc3", "I0", 'chain')


#a.c<-data.frame(rbind(chain1.1, chain2.1, chain3.1))
a.cm<-melt(a.c, id.vars = c('iter','chain'))
aa1<-ggplot(data=a.cm) + facet_wrap( ~ variable, ncol=3, scales="free" ) + geom_line(aes(x=iter, y=value, group=chain, color=chain))
ggsave(paste0("trace",strain.name[strain], hh,".tiff"),device='tiff',plot=aa1, width=54,height=24,units='cm')

apply(chain1, 2, median)
if(pull==F)
{
apply(chain2, 2, median)
apply(chain3, 2, median)
}
}



#####################------------------  




library(devtools)
install_github("bobverity/RgeoProfile")
library(RgeoProfile)
library(RCurl)
library(gdata)

u = 'http://maps.googleapis.com/maps/api/staticmap?center=51.517478,-0.103755&zoom=12&size=640x640&scale=2&maptype=hybrid&language=en-EN&sensor=false'
o = getURLContent(u, verbose = TRUE, useragent = getOption("HTTPUserAgent")) 


ss2<- as.matrix(ss[complete.cases(ss)==T,])
d <- geoData(as.numeric(ss2[,3]), as.numeric(ss2[,4]))

# set model and MCMC parameters
p <- geoParams(data = d, sigma_mean = 9,sigma_var= 6,sigma_squared_shape = 2, chains = 3, 
               burnin = 1e5, samples = 1e5)

#------------------------------------------------------------------
# run model
#------------------------------------------------------------------
# run MCMC
breaks_lon <- seq(46.7,47.3,l=101)
breaks_lat <- seq(-122.6,-123,l=101)
m <- geoMCMC(data = d, params = p, lambda=NULL)




TH_mask <- north_london_mask[which(north_london_mask$NAME == "Tower Hamlets"),]
prob_masked <- geoMask(probSurface = m$posteriorSurface, params = p, mask = TH_mask,
                       operation = "outside", scaleValue = 1e-9)
gp_masked <- geoProfile(prob_masked$prob)

# plot new surface
mapMask <- geoPlotMap(params = p, data = d, source = s, surface = gp_masked, 
                      breakPercent = seq(0,25,l = 11))
mapMask

# hs of masked surface
hs_mask <- geoReportHitscores(params = p, source = s, surface = gp_masked)
hs_mask
setwd('/Users/Natasha/Downloads/', sheet=1)
ss<-read.csv('Cat map points.csv')
library(rgdal)
s1 <- geoShapefile('/Users/Natasha/Downloads/Thurston_ZipCodes')
