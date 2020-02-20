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
library(matrixStats)
library(ggforce)
library(ks)
library(pryr)
library(ggthemes) # Load
library(odin)


rm(list = ls())


#We run realization of epidemiological model gleaned from the MCMC for a status quo estimate, then rerun the model under each of the three interventions. Next we calclate the total cost of strategy and the total qalys lost. After all strategies have had simulations order stategies by total cost or total qaly loss. Calculate the icer for each of ordered strategies and delete dominated strategies. From that framework you get Icers and willingness-to-pay.


#install 
##################################################
####### FILE LOADING FUNCTIONS FOR GRAPHS ########
#################################################
#loaders for all files
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/SUPPORT_ICER_loaders.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/SUPPORT_ICER_conversion_intl.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/OUTPUT_ICER_sampler.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/OUTPUT_post.sampler_intl.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/Output_vaccine_doses_INTL.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_ICER_conversion_22ages.R')
vstrategy<-dget('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_cov_strategy_uk.R')


##coverage rates
load('cov.function55.RData')
load('cov.function70.RData')
load('cov.function30.RData')


load('Intl_general_parameters.RData') ##grab saved workspace
##load general parameters-------------------

load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/ili.counts.rda')
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/virological.rda')

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
death.risk.tables<-list.files(pattern=glob2rx('tab_risk_death_*')); #loads in the order strain B,H1,H3
hosp.risk.tables<-list.files(pattern=glob2rx('tab_risk_hosp*'));
GP.risk.tables<-list.files(pattern=glob2rx('tab_risk_GP*'))
strain.name<-c('H1N1','H3N2','B') ##put all labels in


####INTL FUNCTION Start=================================================
strain<-strainpull<-1 ##set strain
i.country<-country<-4  #set country
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))


##### loaders


cov.sens.loader<-function(strain, i.country, end.cov, version)
{
  
  #if(i.country==4) {setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB/coverage/0.55')}else{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
  #Rearrange tables to correct order of interventions------------------
  sq.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStratStatusQuo', end.cov, strain.name[strain], version)))
  sens.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStrat*','*chool', end.cov, strain.name[strain], version)))
  
  if(length(sens.tables)==0){
    stop('missing SensStrat intervention files');
  }
  if(length(sq.tables)==0)
  {stop('missing SensStrat SQ files')}  
  
 
  load(sens.tables[1]); avg.incidence1<-total.keep; #preschool
  load(sens.tables[2]); avg.incidence2<-total.keep; #primary
  load(sens.tables[3]); avg.incidence3<-total.keep; #secondary
  load(sq.tables[1]); avg.incidence4<-total.keep; #status quo
  
  
  var.name<-paste0('intl.table',strain.name[strain])
  
  assign(var.name, list(avg.incidence4,avg.incidence1,avg.incidence2, avg.incidence3), envir = .GlobalEnv)
  #save(tableH1,file=paste0('AnnualIncH1'))
}



#############################################################################################
####### 5 Load samples to make graphs ########
############################################################################################
inv.names2<-c('Status Quo','Preschool', 'Primary', 'Secondary')
strainpull<-1
#sens.loader(strain=strainpull,i.country=country);
cov.sens.loader(strain=strainpull,i.country=country, end.cov=0.55, version = 'v1')
#############################################################################################
####### 6 QALY Conversion Equations ########
############################################################################################
strain<-strainpull
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
if(strain==1) {load(death.risk.tables[2]);load(hosp.risk.tables[2]);load(GP.risk.tables[2])}
if(strain==2) {load(death.risk.tables[3]);load(hosp.risk.tables[3]);load(GP.risk.tables[3])}
if(strain==3) {load(death.risk.tables[1]);load(hosp.risk.tables[1]);load(GP.risk.tables[1])}
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))

#load all conversion functions


##################################
####### 7 CALCULATE ICERs ########
#################################

####### 7.1 Single strain ICER calculation ########
Dataset1<-NULL
Dataset2<-NULL
Dataset3<-NULL
Dataset4<-NULL


if(strain==1)
{ Dataset1<-intl.tableH1N1[[1]]
  Dataset2<-intl.tableH1N1[[2]]
  Dataset3<-intl.tableH1N1[[3]]
  Dataset4<-intl.tableH1N1[[4]]}
if(strain==2)
  {Dataset1<-intl.tableH3N2[[1]]
  Dataset2<-intl.tableH3N2[[2]]
  Dataset3<-intl.tableH3N2[[3]]
  Dataset4<-intl.tableH3N2[[4]]}
if(strain==3)
  {Dataset1<-intl.tableB[[1]]
  Dataset2<-intl.tableB[[2]]
  Dataset3<-intl.tableB[[3]]
  Dataset4<-intl.tableB[[4]]}


strainpull<-2
strain<-strainpull
dcount.l<-c(3.5)

for(strainpull in 1:length(strain.name))
{
  cov.sens.loader(strain=strainpull,i.country=country, end.cov=0.55, version='v1')
  
for( dd in 1:length(dcount.l))
{
intl.CEA.allstrains(Dataset1 = Dataset1, Dataset2 = Dataset2, Dataset3=Dataset3, Dataset4=Dataset4, strain = strainpull, i.country=country, dcount = dcount.l[dd], program1=1, program2=2, program3 = 3, program4=4, i.cov=fiftyfive.cov)
}
}


GP.costs.l<-Hosp.costs.l<-sum.costs.l<-aged.cases.l<-QALY.tot.l<-QALY.death.l<-cumi.l<-add.dose<-avt.dose<-list()

for(country in 1:length(fname))
{

    cov.sens.loader(strain=strainpull,i.country=country, end.cov=0.55, version='v1')
  
  Dataset1<-NULL
  Dataset2<-NULL
  Dataset3<-NULL
  Dataset4<-NULL
  
  
  if(strainpull==1)
  { Dataset1<-intl.tableH1N1[[1]]
  Dataset2<-intl.tableH1N1[[2]]
  Dataset3<-intl.tableH1N1[[3]]
  Dataset4<-intl.tableH1N1[[4]]}
  if(strainpull==2)
  {Dataset1<-intl.tableH3N2[[1]]
  Dataset2<-intl.tableH3N2[[2]]
  Dataset3<-intl.tableH3N2[[3]]
  Dataset4<-intl.tableH3N2[[4]]}
  if(strainpull==3)
  {Dataset1<-intl.tableB[[1]]
  Dataset2<-intl.tableB[[2]]
  Dataset3<-intl.tableB[[3]]
  Dataset4<-intl.tableB[[4]]}
  
  

strain<-strainpull
dcount.l<-c(3.5)
dd<-1
  intl.AgeStuc.CEA.outcomes(Dataset1 = Dataset1, Dataset2 = Dataset2, Dataset3=Dataset3, Dataset4=Dataset4, strain = strainpull,i.country=country, dcount = dcount.l[dd], program1=1, program2=2, program3 = 3, program4=4, i.cov=mL55, version = 'v1')

}



######################################################################################################
#T-tests and KS tests
########################################################################################################
strain.choice<-1 ##we only examine 'H1N1' in this analysis
version<-'v1'


ks.posterior<-function(i.country, version)
{
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata','GB'));
flu.samp.gb<-list.files(pattern=glob2rx(paste0('GBPostSample',version, strain.name[strain.choice])))

load(flu.samp.gb); pars.dist.gb<-post.sample
#pars.dist.gb<-do.call(rbind.data.frame, post.sample.gb)
post.sample<-c()

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain.choice])))
load(flu.samp)

pars.dist<-post.sample
ks.yes<-list()
#pars.dist<-do.call(rbind.data.frame, post.sample)

#ks.test(post.sample.gb[[3]][,1], post.sample[[3]][,1])
for(i.season in 1:14)
{
    eps1<-ks.test(pars.dist.gb[[i.season]][,1], pars.dist[[i.season]][,1])
    eps2<-ks.test(pars.dist.gb[[i.season]][,2], pars.dist[[i.season]][,2])
    eps3<-ks.test(pars.dist.gb[[i.season]][,3], pars.dist[[i.season]][,3])
    psi<-ks.test(pars.dist.gb[[i.season]][,4], pars.dist[[i.season]][,4])
    q<-ks.test(pars.dist.gb[[i.season]][,5], pars.dist[[i.season]][,5])
    sus1<-ks.test(pars.dist.gb[[i.season]][,6], pars.dist[[i.season]][,6])
    sus2<-ks.test(pars.dist.gb[[i.season]][,7], pars.dist[[i.season]][,7])
    sus3<-ks.test(pars.dist.gb[[i.season]][,8], pars.dist[[i.season]][,8])
    I0<-ks.test(pars.dist.gb[[i.season]][,9], pars.dist[[i.season]][,9])
    
    ks.yes[[i.season]]<-data.frame(cbind( c("eps1", "eps2", "eps3", "psi", "q","susc1", "susc2", "susc3", "I0")
                             , rep('GB', 9), rep(fname[i.country], 9), c(eps1$p.value, eps2$p.value, eps3$p.value, psi$p.value, q$p.value, sus1$p.value, sus2$p.values, sus3$p.value, I0$p.value)))
    
   # ifelse((i.season==1), ks.comp[[i.season]]<-ks.yes, ks.comp<-rbind(ks.comp, ks.yes))
  }
  
  return(ks.yes)
}



##### 
par.c<-list()
par.names<-c("eps1", "eps2", "eps3", "psi", "q","susc1", "susc2", "susc3", "I0")
for(i.country in 1:(length(fname)-1))
{
  for(par in 1:length(par.names))
  {
LUpars<-do.call(rbind.data.frame, ks.posterior(i.country, 'v1'))
hhy<-sum(as.numeric(as.character(LUpars[LUpars$X1==par.names[par],]$X4))>0.05)/14

par.r<-c(fname[i.country], par.names[par], as.numeric(hhy))

ifelse(par==1, par.r2<-par.r, par.r2<-rbind(par.r2, par.r))
  }
  
  par.c[[i.country]]<-par.r2
}

ks.posterior(1, 'v1')
ks.posterior(3, 'v1')
ks.posterior(2, 'v1')
ks.posterior(5, 'v1')
ks.posterior(6, 'v1')
ks.posterior(7, 'v1')
ks.posterior(8, 'v1')
ks.posterior(9, 'v1')
ks.posterior(10, 'v1')


t.test.posterior<-function(i.country, version)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata','GB'));
  flu.samp.gb<-list.files(pattern=glob2rx(paste0('GBPostSample',version, strain.name[strain.choice])))
  
  load(flu.samp.gb); pars.dist.gb<-post.sample
  #pars.dist.gb<-do.call(rbind.data.frame, post.sample.gb)
  post.sample<-c()
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain.choice])))
  load(flu.samp)
  
  pars.dist<-post.sample
  ks.yes<-list()
  #pars.dist<-do.call(rbind.data.frame, post.sample)
  
  #ks.test(post.sample.gb[[3]][,1], post.sample[[3]][,1])
  for(i.season in 1:14)
  {
    eps1<-t.test(pars.dist.gb[[i.season]][,1], pars.dist[[i.season]][,1])
    eps2<-t.test(pars.dist.gb[[i.season]][,2], pars.dist[[i.season]][,2])
    eps3<-t.test(pars.dist.gb[[i.season]][,3], pars.dist[[i.season]][,3])
    psi<-t.test(pars.dist.gb[[i.season]][,4], pars.dist[[i.season]][,4])
    q<-t.test(pars.dist.gb[[i.season]][,5], pars.dist[[i.season]][,5])
    sus1<-t.test(pars.dist.gb[[i.season]][,6], pars.dist[[i.season]][,6])
    sus2<-t.test(pars.dist.gb[[i.season]][,7], pars.dist[[i.season]][,7])
    sus3<-t.test(pars.dist.gb[[i.season]][,8], pars.dist[[i.season]][,8])
    I0<-t.test(pars.dist.gb[[i.season]][,9], pars.dist[[i.season]][,9])
    
    ks.yes[[i.season]]<-data.frame(cbind( c("eps1", "eps2", "eps3", "psi", "q","susc1", "susc2", "susc3", "I0")
                                          , rep('GB', 9), rep(fname[i.country], 9), c(eps1$p.value, eps2$p.value, eps3$p.value, psi$p.value, q$p.value, sus1$p.value, sus2$p.values, sus3$p.value, I0$p.value)))
    
    # ifelse((i.season==1), ks.comp[[i.season]]<-ks.yes, ks.comp<-rbind(ks.comp, ks.yes))
  }
  
  return(ks.yes)
}



###loop to average all results
par.t<-list()
par.names<-c("eps1", "eps2", "eps3", "psi", "q","susc1", "susc2", "susc3", "I0")
for(i.country in 1:(length(fname)-1))
{
  for(par in 1:length(par.names))
  {
    LUpars<-do.call(rbind.data.frame, t.test.posterior(i.country, 'v1'))
    hhy<-sum(as.numeric(as.character(LUpars[LUpars$X1==par.names[par],]$X4))>0.05)/14
    
    par.r<-c(fname[i.country], par.names[par], as.numeric(hhy))
    
    ifelse(par==1, par.r2<-par.r, par.r2<-rbind(par.r2, par.r))
  }
  
  par.t[[i.country]]<-par.r2
}


t.test.posterior(1, 'v1')
t.test.posterior(3, 'v1')








########################################################################################
#Start of graphs
########################################################################################
post.pars.out<-c()

par.grapher<-function(i.country, strain.choice, version)
{
  post.pars<-c()
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
    
    post.pars<-cbind(sample.m, rep(i.season, dim(sample.m)[1]), rep(fname[i.country], dim(sample.m)[1]))
    
    ifelse(i.country==1 & i.season==1, post.pars.out<<-post.pars, post.pars.out<<-rbind(post.pars.out, post.pars))
    
  }
}   

for(i.country in 1:(length(fname)-1))
{par.grapher(i.country=i.country, strain.choice = 1, version = 'v1')}



#####################################################################################################################

library(randomcoloR)
n <- 10
palette <- distinctColorPalette(n)
palette<-c("#D6CCCC", "#8CDDBF", "#D4956D", "#BD93D8", "#7353E0", "#79ADD8", "#8AE36F", "#E0DA6A", "#DE5ADD", "#DD6392")
setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")


#susceptibility
colnames(post.pars.out)<-c('num','parameter','value','season','country')
post.pars.out$country<-factor(post.pars.out$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
susc.pars<-post.pars.out[post.pars.out$parameter=='susc1'|post.pars.out$parameter=='susc2'|post.pars.out$parameter=='susc3',]

susc.ffm<-ddply(susc.pars, .(parameter, country), summarize, ms=mean(value), LB=quantile(value, 0.025), UB=quantile(value, 0.975))
susc.ffm<-cbind(susc.ffm,rep('Susceptibility',dim(susc.ffm)[1]))
colnames(susc.ffm)<-c('parameter', 'country','ms', 'LB', 'UB', 'Susceptibility')
sus.p<-ggplot(susc.ffm, aes(x=country, y=ms, group=parameter, color=country))+theme_bw()+scale_color_manual(values=palette)+geom_pointrange(aes(ymax=UB, ymin=LB, shape=parameter), inherit.aes = T, position=position_dodge(0.5))+geom_hline(yintercept=0)+labs(x='Contact Matrix', y='Parameter Value')+geom_errorbar(aes(ymin = LB, ymax = UB),position=position_dodge(0.5), width = 0.2)+guides(color=FALSE)+theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position='right', legend.title=element_blank(),legend.text=element_text(size=14))+facet_grid(.~Susceptibility)+scale_shape_discrete(labels=c('susc1'=expression(sigma['6m-14 yrs old']), 'susc2'=expression(sigma['15-64 yrs old']), 'susc3'=expression(sigma['65+ yrs old   '])))



##############################epsilons
eps.pars<-post.pars.out[post.pars.out$parameter=='eps1'|post.pars.out$parameter=='eps2'|post.pars.out$parameter=='eps3',]

eps.ffm<-ddply(eps.pars, .(parameter, country), summarize, ms=mean(value), LB=quantile(value, 0.025), UB=quantile(value, 0.975))

eps.ffm<-cbind(eps.ffm,rep('Case Ascertainment',dim(eps.ffm)[1]))
colnames(eps.ffm)<-c('parameter', 'country','ms', 'LB', 'UB', 'Case Ascertainment')

eps.p<-ggplot(eps.ffm, aes(x=country, y=ms, group=parameter, color=country))+theme_bw()+scale_color_manual(values=palette)+geom_pointrange(aes(ymax=UB, ymin=LB, shape=parameter), inherit.aes = T, position=position_dodge(0.5))+geom_hline(yintercept=0)+labs(x='Contact Matrix', y='Parameter Value')+geom_errorbar(aes(ymin = LB, ymax = UB),position=position_dodge(0.5), width = 0.2)+guides(color=FALSE)+scale_shape_discrete(labels=c('eps1'=expression(epsilon['6m-15 yrs old']), 'eps2'=expression(epsilon['16-64 yrs old']), 'eps3'=expression(epsilon['65+ yrs old   '])))+theme(axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position='right',legend.title=element_blank(), legend.text=element_text(size=14))+facet_grid(.~`Case Ascertainment`)




#"I0"    "psi"   "q" 
rem.pars<-post.pars.out[post.pars.out$parameter=='I0'|post.pars.out$parameter=='psi'|post.pars.out$parameter=='q',]
rem.ffm<-ddply(rem.pars, .(parameter, country), summarize, ms=mean(value), LB=quantile(value, 0.025), UB=quantile(value, 0.975))

levels(rem.ffm$parameter)<-c(levels(rem.ffm$parameter), 'Probability of Outside Infection', 'Virus Transmissibility','Initial Infections Coefficient')

rem.ffm$parameter[rem.ffm$parameter=='psi']<-'Probability of Outside Infection'
rem.ffm$parameter[rem.ffm$parameter=='q']<-'Virus Transmissibility'
rem.ffm$parameter[rem.ffm$parameter=='I0']<-'Initial Infections Coefficient'



rem.p<-ggplot(rem.ffm, aes(x=country, y=ms, group=parameter, color=country))+theme_bw()+scale_color_manual(values=palette)+facet_wrap(~parameter, scale='free', ncol=3)

rem.p2<-rem.p+geom_pointrange(aes(ymax=UB, ymin=LB), inherit.aes = T, position=position_dodge(0.5))+geom_hline(yintercept=0)+labs(x='Contact Matrix')+geom_errorbar(aes(ymin = LB, ymax = UB),position=position_dodge(0.5), width = 0.2)+guides(color=FALSE)+theme(axis.title.y=element_blank())


library(ggpubr)
pars.out<-ggarrange(sus.p, eps.p,rem.p2, common.legend=F,legend='right',ncol=1, nrow=3)
pars.out2<-annotate_figure(pars.out,
                         top = text_grob("Posterior Parameters"))

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('paroutputintl','3.5',version,".png"),plot = pars.out2, width=11, height=8.5, units='in',device='png')




#################################################################
#Health Outcomes Graphs
##################################################################




GP.costs.l<-do.call(rbind.data.frame, GP.costs.l)
colnames(GP.costs.l)<-c('value','country','Strategy')
GP.costs.l$Strategy<-factor(GP.costs.l$Strategy, levels = c("Preschool",  "Primary",    "Secondary",  "Status Quo"), labels=c('I1', 'I2', 'I3', 'SQ'))
GP.costs.l$country<-factor(GP.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(GP.costs.l$country)

GP.costs<-ddply(GP.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value)))/1e3,
       LB=quantile(as.numeric(as.character(value)), 0.025)/1e3, UB=quantile(as.numeric(as.character(value)), 0.975)/1e3)


GPcost.graph<-ggplot(GP.costs.l, aes(as.factor(Strategy), (as.numeric(as.character(value)))/1e6, fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_boxplot()+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+facet_wrap(~country, scale='fixed', ncol=3, nrow=4)+theme_bw()+coord_flip()+labs(fill='Intervention\n Strategty', y='Annual Cost of GP Consults (£ Millions)', x='Intervention Strategy')+geom_hline(yintercept = 0, color='dark gray')+theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('GPcosts','3.5',version,".png"),plot = GPcost.graph, width=11, height=8.5, units='in',device='png')

##Total Costs
sum.costs.l<-do.call(rbind.data.frame, sum.costs.l)
colnames(sum.costs.l)<-c('value','country','Strategy')
levels(sum.costs.l$Strategy)<-c('I1: 2-4 years old', 'I2: 5-11 years old', 'I3: 12-16 years old', 'SQ: High Risk & 65+')
sum.costs.l$country<-factor(sum.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(sum.costs.l$country)

sum.costs<-ddply(sum.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))/1e6),
    LB=quantile(as.numeric(as.character(value))/1e6, 0.025), UB=quantile(as.numeric(as.character(value))/1e6, 0.975))


influenza.costs.graph<-ggplot(sum.costs.l, aes(x=as.factor(country), y=((as.numeric(as.character(value)))/1e6), fill=as.factor(country)))+geom_violin(scale='width', draw_quantiles = 0.5 )+facet_wrap(~Strategy, ncol=4, nrow=2)+theme_bw()+labs(fill='Contact Matrix', y='Annual Cost of Influenza Illness (£ Millions)', x='Contact Matrix')+scale_fill_manual(values=palette)+scale_y_continuous(limits = c(125, 220))

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('influenzacosts','3.5',version,".png"),plot = influenza.costs.graph, width=11, height=8.5, units='in',device='png')

#scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))

#=============================QALYS total
QALY.tot.l<-do.call(rbind.data.frame, QALY.tot.l)
colnames(QALY.tot.l)<-c('value','country','Strategy')
levels(QALY.tot.l$Strategy)<-c('I1: 2-4 years old', 'I2: 5-11 years old', 'I3: 12-16 years old', 'SQ: High Risk & 65+')
QALY.tot.l$country<-factor(QALY.tot.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(QALY.tot.l$country)

QALY.tot.sum<-ddply(QALY.tot.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))),
LB=quantile(as.numeric(as.character(value)), 0.025), UB=quantile(as.numeric(as.character(value)), 0.975))



QALYs.tot.graph<-ggplot(QALY.tot.sum, aes(x=as.factor(country), y=(m1), fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=6))+facet_wrap(~Strategy, scale='free', ncol=4, nrow=3)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Hospitalization Costs (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(country), ymin=m1, ymax=UB), width=.2,  position='dodge', inherit.aes = F)

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('hospitalcosts','3.5',version,".png"),plot = hosp.costs.graph, width=11, height=8.5, units='in',device='png')

#despite similar costs, QALY losses differ


##QALY loss from Death
QALY.death.l<-do.call(rbind.data.frame, QALY.death.l)
levels(QALY.death.l$Strategy)
colnames(QALY.death.l)<-c('value','country','Strategy')
QALY.death.l$Strategy<-factor(QALY.death.l$Strategy, levels=c("I1", "I2", "I3", "SQ"))
QALY.death.l$country<-factor(QALY.death.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(QALY.death.l$country)

QALY.death.sum<-ddply(QALY.death.l, .(country, Strategy), summarise, m1=mean(round(as.numeric(as.character(value)))),
  LB=quantile(round(as.numeric(as.character(value))), 0.025), UB=quantile(round(as.numeric(as.character(value))), 0.975))


QALYs.lost.graph<-ggplot(QALY.death.sum, aes(x=as.factor(country), y=(m1/1e6), fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=6))+facet_wrap(~Strategy, scale='free', ncol=4, nrow=3)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Hospitalization Costs (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(country), ymin=m1/1e6, ymax=UB/1e6), width=.2,  position='dodge', inherit.aes = F)

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('hospitalcosts','3.5',version,".png"),plot = hosp.costs.graph, width=11, height=8.5, units='in',device='png')



###############################Symptomatic Cases
aged.cases.l<-do.call(rbind.data.frame, aged.cases.l)
colnames(aged.cases.l)<-c('value','country','Strategy')
levels(aged.cases.l$Strategy)<-c( 'SQ','I1', 'I2', 'I3')
aged.cases.l$country<-factor(aged.cases.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(aged.cases.l$country)

aged.sum<-ddply(aged.cases.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))),
UB=quantile(as.numeric(as.character(value)), 0.975), LB=quantile(as.numeric(as.character(value)), 0.025))


ggplot(aged.cases.l, aes(as.factor(country), (as.numeric(as.character(value)))/100000, fill=as.factor(country)))+scale_fill_manual(values=palette)+geom_violin(draw_quantiles = c(0.5), scale='width')+scale_y_continuous(breaks=scales::pretty_breaks(n=8))+facet_wrap(~Strategy, scale='free', ncol=2, nrow=2)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Symptomatic Cases (Hundred  Thousands)', x='Contact Matrix')


############Hospital Costs
Hosp.costs.l<-do.call(rbind.data.frame, Hosp.costs.l)
colnames(Hosp.costs.l)<-c('value','country','Strategy')


Hosp.costs.l$Strategy<-factor(Hosp.costs.l$Strategy, levels = c("Preschool",  "Primary",    "Secondary",  "Status Quo"), labels=c('I1', 'I2', 'I3', 'SQ'))
Hosp.costs.l$country<-factor(Hosp.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
#col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(Hosp.costs.l$country)
hosp.sum<-ddply(Hosp.costs.l, .(country, Strategy), summarise, m1=mean(round(as.numeric(as.character(value))))/1000,
LB=quantile(round(as.numeric(as.character(value))), 0.025)/1000, UB=quantile(round(as.numeric(as.character(value)))/1000, 0.975))
hosp.costs.graph<-ggplot(hosp.sum, aes(x=as.factor(country), y=(m1/1e6), fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=6))+facet_wrap(~Strategy, scale='free', ncol=4, nrow=3)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Hospitalization Costs (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(country), ymin=m1/1e6, ymax=UB/1e6), width=.2,  position='dodge', inherit.aes = F)

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('hospitalcosts','3.5',version,".png"),plot = hosp.costs.graph, width=11, height=8.5, units='in',device='png')





#####Cumulative Incidence
ks.compare<-function(xx, data.test)
{
  data.l<-data.test
  colnames(data.l)<-c('value','country','Strategy')
  levels(data.l$Strategy)<-c('SQ','S1', 'S2', 'S3')
  data.l$country<-factor(data.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
for(yy in 1:(length(fname)-1))
{
  if(xx==yy) next
inv1<-ks.test(as.numeric(as.character(data.l$value[data.l$country==fname[xx]& data.l$Strategy=='S1'])), as.numeric(as.character(data.l$value[data.l$country==fname[yy]& data.l$Strategy=='S1'])))

inv2<-ks.test(as.numeric(as.character(data.l$value[data.l$country==fname[xx]& data.l$Strategy=='S2'])), as.numeric(as.character(data.l$value[data.l$country==fname[yy]& data.l$Strategy=='S2'])))

inv3<-ks.test(as.numeric(as.character(data.l$value[data.l$country==fname[xx]& data.l$Strategy=='S3'])), as.numeric(as.character(data.l$value[data.l$country==fname[yy]& data.l$Strategy=='S3'])))

invq<-ks.test(as.numeric(as.character(data.l$value[data.l$country==fname[xx]& data.l$Strategy=='SQ'])), as.numeric(as.character(data.l$value[data.l$country==fname[yy]& data.l$Strategy=='SQ'])))

ks.yes<-data.frame(cbind(c('SQ', 'S1','S2','S3'), rep(fname[xx], 4), rep(fname[yy], 4), c(invq$p.value, inv1$p.value, inv2$p.value, inv3$p.value)))

ifelse((yy==1), ks.comp<-ks.yes, ks.comp<-rbind(ks.comp, ks.yes))
}
  
  return(ks.comp)
}

ks.compare(10, cumi.l)


program1<-1
program2<-2
program3<-3
program4<-4

hyy<-list()
for(i.country in 1:length(fname))
{
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
load(file=paste0(fname[i.country],'CEAcosts',strain.name[strain],'p', program1, program2, program3,program4,'d', dcount.l[dd]))
hyy[[i.country]]<-coststab

}



##########################################################################################
#Epi Pulls for Age-stratified Deaths, Cases, GP Appointments, and Hospitalizations
#####################################################################################################################
cols2<-c("#FF61CC","#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF" )


######loader
age.strat.loader<-function(country, vout, version, cov.var, dcount, strategy)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #strategy.summary<-NULL
  age.tables<-list.files(pattern=glob2rx(paste0(fname[country],'agestrat.',vout,'.allstrainp*',cov.f[cov.var],dcount.l[dcount],version)))
  load(age.tables[strategy], envir = .GlobalEnv)
} 

long.inv.names<-c("Preschool"  ,"Primary"  ,  "Secondary", 'Status Quo')

indirect.strat.loader<-function(country, version, cov.var, strain, strategy)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[country])))
  #strategy.summary<-NULL
 if(strategy>1)
   {inv.ID.tables<-list.files(pattern=glob2rx(paste0(fname[country],'unvacc.incidenceIR*chool', cov.var,strain.name[strain],version)))
  load(inv.ID.tables[strategy-1], envir = .GlobalEnv)
  }else{
  inv.SQ.tables<-list.files(pattern=glob2rx(paste0(fname[country],'unvacc.incidenceIRSt*', cov.var,strain.name[strain],version)))
  load(inv.SQ.tables[1], envir = .GlobalEnv)}
  #unvacc.SQ<-unvacc.keep}
  
  
} 



indirect.strat.loader(1,'v1',0.55,1,1)
#########v.out=variable (e.g. death)

epistack.indirect<-function(country, i.version, cov.var, strain)
{
  agg.hosp<-c()
  age.small<-list()
  for(i.program in 1:length(long.inv.names))
  {
    
    risk.age<-list()
    prop.small<-list()
    redis<-c()
      indirect.strat.loader(country=country, version=i.version,cov.var, strain, strategy=i.program)
      assign(paste0('unvacc.keep', i.program), unvacc.keep)
      redis<-get(paste0('unvacc.keep', i.program), envir = .GlobalEnv)
      redis<-Reduce('+', redis,)/ length(redis)
      
      #redis<-aged.cases[[iter]]
      
      age.small[[i.program]]<-cbind(redis[,1]+redis[,2]+redis[,12]+redis[,13], #0-2, 
                               redis[,3]+redis[,14], #2-4,
                               redis[,4]+redis[,15], #5-11
                               redis[,5]+redis[,6]+redis[,16]+redis[,17],
                               redis[,7]+redis[,18],
                               redis[,8]+redis[,19]+redis[9]+redis[,20]+redis[10]+redis[,21],
                               redis[11]+redis[,22]#65+
      )
      colnames(age.small[[i.program]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25-64 yr olds', '65+ yr olds')
      
    }
    
    test<-melt(age.small)
    test$L1<-factor(test$L1)
    test$L1<-factor(test$L1, levels=c('1','2','3','4'),  labels=c('SQ','I1', 'I2', 'I3'))
  
  agg.hosp<-test
  colnames(agg.hosp)<-c('X1','Age','value','Strategy')
  
  hosp.m2<-ddply(agg.hosp[,2:4], .(Age, Strategy), summarise, m1=mean(value), LB=quantile(value, 0.025, na.rm=F), UB=quantile(value, 0.975, na.rm=F))
  #hosp.m2$CM<-factor(hosp.m2$CM, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
  #relevel
  hosp.m2$Age<-factor(hosp.m2$Age, levels=c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25-64 yr olds', '65+ yr olds'))
      
      epistack<-ggplot(data=hosp.m2)+geom_col(aes(x=Strategy,y=m1/1e6, fill=Age))+theme_bw()+scale_fill_brewer()+
        labs(subtitle=paste(fname[country],'Coverage=', paste0(cov.var*100,'%')), fill = 'Age Strata')+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12), axis.text.y = element_text(size=14),axis.title.y = element_text(size=12), legend.position = 'right')

      
      
      return(epistack)
  ##savers
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
  #ggsave(paste0('epistack',country, i.vout,cov.f[cov.var],strain.name[strain],i.version,".png"),plot =epistack , width=11, height=8.5, units='in',device='png')
}  

stack.out<-ggarrange(epistack.indirect(country=1,cov.var = 0.55,strain = 1,i.version = 'v1')+
                       rremove("x.text")+rremove('x.ticks')+rremove('xy.title'),
                    epistack.indirect(country=2,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("xy.text")+rremove('ticks')+rremove('xy.title'),
                    epistack.indirect(country=3,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("x.text")+rremove('x.ticks')+rremove('xy.title'),
                    epistack.indirect(country=4,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("xy.text")+rremove('ticks')+rremove('xy.title'),
                    epistack.indirect(country=5,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("x.text")+rremove('x.ticks')+rremove('xy.title'),
                    epistack.indirect(country=6,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("xy.text")+rremove('ticks')+rremove('xy.title'),
                    epistack.indirect(country=7,cov.var = 0.55,strain = 1,i.version = 'v1')+
                      rremove("x.text")+rremove('x.ticks')+rremove('xy.title'),
                    epistack.indirect(country=8,cov.var = 0.55,strain = 1,i.version = 'v1')+rremove('xy.title'),
                    epistack.indirect(country=9,cov.var = 0.55,strain = 1,i.version = 'v1')+rremove('xy.title'),
                    common.legend=T,legend='bottom',ncol=2, nrow=5)



stack.out2<-annotate_figure(stack.out,
                           top = text_grob('Average Incidence Among Unvaccinated'),
                           left='Average Count (Millions)',
                           bottom='Intervention Strategy')

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('stackIDoutcomeintl',version,".png"),plot = stack.out2, width=8.5, height=11, units='in',device='png')




########################################################################################
####### 8 Calculate Cases Averted Per Dose ########
###################################################################################################




diminish<-function(dds)
{tero<-lapply(1:tab.diff, FUN=function(i) 
  cbind(
    Map('+',Map('+',Map('+', dds[[i]][,1], dds[[i]][,2]),dds[[i]][,12]),dds[[i]][,13]),
    Map('+', dds[[i]][,3], dds[[i]][,14]),
    Map('+', dds[[i]][,4], dds[[i]][,15]),
    Map('+', Map('+', Map('+', dds[[i]][,5], dds[[i]][,6]),  dds[[i]][,16]), dds[[i]][,17]),
    Map('+', dds[[i]][,7], dds[[i]][,18]),
    Map('+', dds[[i]][,8], dds[[i]][,19]),
    Map('+', dds[[i]][,9], dds[[i]][,20]),
    Map('+', Map('+', Map('+',dds[[i]][,10], dds[[i]][,11]), dds[[i]][,21]), dds[[i]][,22])))
return(tero)}

diminish2<-function(dds)
{tero<-lapply(1:tab.diff, FUN=function(i) 
  cbind(
    Map('+',Map('+',Map('+', dds[,1], dds[,2]),dds[,12]),dds[,13]),
    Map('+', dds[,3], dds[,14]),
    Map('+', dds[,4], dds[,15]),
    Map('+', Map('+', Map('+', dds[,5], dds[,6]),  dds[,16]), dds[,17]),
    Map('+', dds[,7], dds[,18]),
    Map('+', dds[,8], dds[,19]),
    Map('+', dds[,9], dds[,20]),
    Map('+', Map('+', Map('+',dds[,10], dds[,11]), dds[,21]), dds[,22])))
return(tero)}

diminish3<-function(dds)
{tero<-cbind.data.frame(
    dds[,1]+dds[,2]+dds[,12]+dds[,13],
    dds[,3]+dds[,14],
    dds[,4]+dds[,15],
    dds[,5]+dds[,6]+dds[,16]+dds[,17],
    dds[,7]+dds[,18],
    dds[,8]+dds[,19],
    dds[,9]+dds[,20],
    dds[,10]+dds[,11]+dds[,21]+dds[,22], stringsAsFactors=F)
return(tero)}

i.cov<-mL55
strain<-2

status.quo.program<-alt.program.coverage(new.coverage = i.cov, strategy = 1, strain = strain)#status quo
new.program2<-alt.program.coverage(new.coverage = i.cov,strategy = 2, strain = strain)
new.program3<-alt.program.coverage(new.coverage = i.cov,strategy = 3, strain = strain)
new.program4<-alt.program.coverage(new.coverage = i.cov,strategy = 4, strain = strain)


#returns number of vaccine doses per season
vaccines.doses.new4<-t(t(new.program4[,2:23])*as.vector(popv[1:22]))
vaccines.doses.new2<-t(t(new.program2[,2:23])*as.vector(popv[1:22]))
vaccines.doses.new3<-t(t(new.program3[,2:23])*as.vector(popv[1:22]))
vaccines.doses.statquo<-t(t(status.quo.program[,2:23])*as.vector(popv[1:22]))


ca.out12<-ca.out13<-ca.out14<-NULL

for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #cov.sens.loader(strain=strainpull,i.country=i.country, end.cov=0.55, version='u2')
  cov.sens.loader(strain=strainpull,i.country=i.country, end.cov=0.55, version='v1')
  
  #cases.avert2<-lapply(1:tab.diff, FUN=function(i) Map('-', calc.cases[[1]], calc.cases[[2]])[[i]]/rowSums(vaccines.doses.new2)[i])
  if(strain==1)
  {inc.avert2<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH1N1[[1]], intl.tableH1N1[[2]])[[i]]/sum(vaccines.doses.new2[i,]))
  inc.avert3<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH1N1[[1]], intl.tableH1N1[[3]])[[i]]/sum(vaccines.doses.new3[i,]))
  inc.avert4<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH1N1[[1]], intl.tableH1N1[[4]])[[i]]/sum(vaccines.doses.new4[i,]))}
  
  if(strain==2)
  {inc.avert2<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH3N2[[1]], intl.tableH3N2[[2]])[[i]]/sum(vaccines.doses.new2[i,]))
  inc.avert3<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH3N2[[1]], intl.tableH3N2[[3]])[[i]]/sum(vaccines.doses.new3[i,]))
  inc.avert4<-lapply(1:tab.diff, FUN=function(i) Map('-', intl.tableH3N2[[1]], intl.tableH3N2[[4]])[[i]]/sum(vaccines.doses.new4[i,]))}
  
  
  ca.m12<-Reduce('+', inc.avert2)/tab.diff #not age structured
  ca.m13<-Reduce('+', inc.avert3)/tab.diff
  ca.m14<-Reduce('+', inc.avert4)/tab.diff
  
  ca.m12<-diminish3(ca.m12)
  ca.m13<-diminish3(ca.m13)
  ca.m14<-diminish3(ca.m14)
  
  ca.f12<-data.frame(rep(fname[i.country], dim(ca.m12)[1]), ca.m12, rowSums(ca.m12), rep(inv.names[[1]], dim(ca.m12)[1]))
  ca.f13<-data.frame(rep(fname[i.country], dim(ca.m13)[1]), ca.m13, rowSums(ca.m13), rep(inv.names[[2]], dim(ca.m13)[1]))
  ca.f14<-data.frame(rep(fname[i.country], dim(ca.m14)[1]), ca.m14, rowSums(ca.m14), rep(inv.names[[3]], dim(ca.m14)[1]))
  
  ifelse(i.country==1, ca.out12<-ca.f12, ca.out12<-rbind(ca.out12, ca.f12))
  ifelse(i.country==1, ca.out13<-ca.f13, ca.out13<-rbind(ca.out13, ca.f13))
  ifelse(i.country==1, ca.out14<-ca.f14, ca.out14<-rbind(ca.out14, ca.f14))
}

colnames(ca.out12)<-c('country','0-1 yrs', '2-4 yrs', '5-11 yrs', '12-16 yrs', '17-24 yrs', '25-44 yrs', '45-64 yrs', '65+ yrs','total', 'program')

colnames(ca.out13)<-colnames(ca.out12)
colnames(ca.out14)<-colnames(ca.out12)

all.interventions<-rbind(rbind(ca.out12, ca.out13), ca.out14)


dds.m<-melt(all.interventions, id.vars = c('country', 'program'))
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
strat.ddsm<-dds.m[!dds.m$variable=='total',]


young.ages<-ggplot(strat.ddsm[strat.ddsm$variable %in% c('0-1 yrs', '2-4 yrs', '5-11 yrs', '12-16 yrs'),], aes(country, as.numeric(as.character(value))*1000, fill=country))+ theme_grey()+
  geom_boxplot()+facet_wrap(variable~program, ncol=3 )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+theme(legend.position = 'none')+ labs(fill='Age Groups')

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
strain<-2
ggsave(paste0('infecaverted.young', strain.name[strain], version,".png"),plot = young.ages, width=11, height=8.5, units='in',device='png')

old.ages<-ggplot(strat.ddsm[strat.ddsm$variable %in% c('17-24 yrs', '25-44 yrs', '45-64 yrs', '65+ yrs'),], aes(country, as.numeric(as.character(value))*1000, fill=country))+ theme_grey()+
  geom_boxplot()+facet_wrap(variable~program, ncol=3 )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+ labs(fill='Age Groups')+theme(legend.position = 'none')
setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
#strain<-1
ggsave(paste0('infecaverted.old', strain.name[strain], version,".png"),plot = old.ages, width=11, height=8.5, units='in',device='png')

total.ages<-ggplot(dds.m[dds.m$variable %in% c('total'),], aes(country, as.numeric(as.character(value))*1000, fill=country))+ theme_grey()+
  geom_boxplot()+facet_wrap(variable~program, ncol=3 )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+ labs(fill='Age Groups')+theme(legend.position = 'none')
setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
#strain<-1
ggsave(paste0('infecaverted.totals', strain.name[strain], version,".png"),plot = total.ages, width=11, height=8.5, units='in',device='png')

averted.totals<-dds.m[dds.m$variable %in% c('total'),]

avert.t<-ddply(averted.totals, .(country,program), summarise, tmu=round(mean(value)*1000), tmuLB95=round((mean(value)-1.95*sd(value))*1000), tmuUB95=round((mean(value)+1.95*sd(value))*1000))
avert.t<-avert.t %>% arrange(program)

infec.avert1<-ggplot(dds.m[dds.m$program=='Preschool',], aes(country, as.numeric(as.character(value))*1000, fill=variable))+ theme_grey()+geom_boxplot()+facet_wrap(~variable, scales ='free_y' )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+ labs(fill='Age Groups', title = 'Strategy I1: Preschool Vaccination')+theme(legend.position = 'none')

infec.avert2<-ggplot(dds.m[dds.m$program=='Primary',], aes(country, as.numeric(as.character(value))*1000, fill=variable))+ theme_grey()+geom_boxplot()+facet_wrap(~variable, scales ='free_y' )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+ labs(fill='Age Groups', title = 'Strategy I2: Primary School Vaccination')+theme(legend.position = 'none')

infec.avert3<-ggplot(dds.m[dds.m$program=='Secondary',], aes(country, as.numeric(as.character(value))*1000, fill=variable))+ theme_grey()+geom_boxplot()+facet_wrap(~variable, scales ='free_y' )+xlab('Contact Matrix')+ylab('Infections Averted per 1000 Vaccine Dose')+ labs(fill='Age Groups', title = 'Strategy I3: Secondary School Vaccination')+theme(legend.position = 'none')

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
strain<-1
ggsave(paste0('infecaverted.preschool', strain.name[strain], version,".png"),plot = infec.avert1, width=11, height=8.5, units='in',device='png')
ggsave(paste0('infecaverted.primary', strain.name[strain], version,".png"),plot = infec.avert2, width=11, height=8.5, units='in',device='png')
ggsave(paste0('infecaverted.sec',strain.name[strain],  version,".png"),plot = infec.avert3, width=11, height=8.5, units='in',device='png')



#####################################
#ICER graphs and calculations
#######################################
library('plyr')
library('ggpubr')

dcount<-3.5
intl.ICER.loader<-function(strain, i.country, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  strategy.summary<-NULL
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.strat', strain.name[strain],'p', '*', 'd',dcount)))
} #outputs as ICER.set2 list, QALY diff is difference from status Quo



inv.names<-c('Strategy I1: 2-4 yr olds', 'Strategy I2: 5-11 yr olds', 'Strategy I3: 12-16 yr olds')
cos.f<-qds.f<-qds.m<-cos.m<-NULL




cep.grapher<-function(i.country, strain, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.intl', strain.name[strain],'p', '*', 'd',dcount)))
  
  cos.f<-qds.f<-qds.m<-cos.m<-NULL
  
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    load(cea.tables[[ff]]) #H1N1
    
    if(ff==1) cer<-ICER.set2 #grab ICERs
    if(ff==2) cer<-ICER.set3 #grab ICERs
    if(ff==3) cer<-ICER.set4 #grab ICERs
    
    #m.qa<-median(qa)
    #m.co<-median(co)
    
    qa1<-unlist(cer$`net QALY`[1:14])
    co1<-unlist(cer$`net cost`[1:14])
    
    #qa1<-Reduce('+', cer$`net QALY`[1:14])/tab.diff
    #co1<-Reduce('+', cer$`net cost`[1:14])/tab.diff
    qds<-data.frame(cbind(rep(fname[i.country], length(qa1)), qa1, rep(inv.names[[ff]], length(qa1))))
    cos<-data.frame(cbind(rep(fname[i.country], length(co1)), co1, rep(inv.names[[ff]], length(co1))))
    
    #qa<-unlist(cer$`net QALY`[1:14])
    #co<-unlist(cer$`net cost`[1:14])
    
    #qds<-data.frame(cbind(rep(fname[i.country], length(qa)), qa, rep(inv.names[[ff]], length(qa))))
    #cos<-data.frame(cbind(rep(fname[i.country], length(co)), co, rep(inv.names[[ff]], length(co))))
    
    colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qds.f<-qds, qds.f<-rbind(qds.f, qds))
    ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    
    meds.f<-cbind(as.numeric(as.character(qds[,2])), as.numeric(as.character(cos[,2])))
    
    xy.points<-data.frame(meds.f)
    names(xy.points)<-c('x', 'y')
    kd <- ks::kde(xy.points, compute.cont=TRUE)
    contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont["10%"])[[1]])
    if((contour_95$y[length(contour_95$y)]!= contour_95$y[1]) | (contour_95$x[length(contour_95$x)]!= contour_95$x[1]))
    {contour_95$x[length(contour_95$x)+1]<-contour_95$x[1]
    contour_95$y[length(contour_95$y)+1]<-contour_95$y[1]}
   
    contour_95 <- data.frame(contour_95, rep(inv.names[[ff]], length(contour_95$x)))
    colnames(contour_95)<-c('level','x','y','program')
    
    p1<-geom_path(data=contour_95, aes(x, y, color=program))
    contour.name<-paste0('cont', fname[i.country], ff)
    assign(contour.name,  p1)
  }
  
  
  colnames(qds.f)<-  c('country', 'net QALYS', 'Intervention')
  colnames(cos.f)<- c('country', 'net costs', 'Intervention')
  
  counter<-c()
  cost.utility<-as.matrix(cbind(qds.f[,1:2], cos.f[,2], qds.f[,3]))
  cost.utility<-data.frame(na.omit(unname(cost.utility)))
  colnames(cost.utility)<-c('country', 'QALYS', 'Costs', 'Intervention')
  
  
  mean.p<-ddply(as.data.frame(cost.utility), .(Intervention), summarize, Rate1=mean(as.numeric(as.character(QALYS))), Rate2=mean(as.numeric(as.character(Costs))), costLB=quantile(as.numeric(as.character(Costs)), 0.025), costUB=quantile(as.numeric(as.character(Costs)), 0.975), qalyLB=quantile(as.numeric(as.character(QALYS)), 0.025), qalyUB=quantile(as.numeric(as.character(QALYS)), 0.975), cost.sd=sd(as.numeric(as.character(Costs))), qal.sd=sd(as.numeric(as.character(QALYS))))
  
  mean.po<-mean.p[order(mean.p$Rate2, decreasing = F),]
  
  not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  
  rudi.icer<-diff(rbind(c(0,0),cbind(mean.po$Rate1, mean.po$Rate2)))*-1
  rudi.icer.ub<-diff(rbind(c(0,0),cbind(mean.po$costUB, mean.po$qalyUB)))*-1
  rudi.icer.lb<-diff(rbind(c(0,0),cbind(mean.po$costLB, mean.po$qalyLB)))*-1
  rudi.icer.sd<-diff(rbind(c(0,0),cbind(mean.po$cost.sd, mean.po$qal.sd)))*-1
  
  div<-rudi.icer[,2]/rudi.icer[,1]
  div.ub<-rudi.icer.ub[,1]/rudi.icer.ub[,2]
  div.lb<-rudi.icer.lb[,1]/rudi.icer.lb[,2]
  div.sd<-rudi.icer.sd[,1]/rudi.icer.sd[,2]
  
  
  flip<-as.matrix(mean.po[which(div >0&div <20000),])
  
  if(dim(flip)[1]<2)
  {
    mean.f<-data.frame(t(as.matrix(c(flip[1:9], 
                                     round(div[which(div >0&div <20000)]), 
                                     round(div.ub[which(div >0&div <20000)]), #quantile
                                     round(div.lb[which(div >0&div <20000)]), #quantile
                                     round(div[which(div >0&div <20000)]-1.96*div.sd[which(div >0&div <20000)]), round(div[which(div >0&div <20000)]+1.96*div.sd[which(div >0&div <20000)]))))) #sd
  }else{
    mean.f<-data.frame(cbind(as.matrix(flip[,1:9]),
                             round(div[which(div >0&div <20000)]), 
                             round(div.ub[which(div >0&div <20000)]), #quantile
                             round(div.lb[which(div >0&div <20000)]), #quantile
                             round(div[which(div >0&div <20000)]-1.96*div.sd[which(div >0&div <20000)]), round(div[which(div >0&div <20000)]+1.96*div.sd[which(div >0&div <20000)]))) #sd
  }  
  
  sq<-matrix(c('SQ: Risk Groups, 65+ yrs old',rep(0,13)), nrow=1)
  
  if(dim(na.omit(mean.f))[2]>0 & dim(na.omit(mean.f))[1]>0)
  { colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  levels(mean.f$Intervention)= c(levels(mean.f$Intervention),'SQ: Risk Groups, 65+ yrs old')
  mean.ff<-data.frame(as.matrix(mean.f[order(mean.f$Rate2, decreasing=F),]))
  colnames(mean.ff)<-colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  }else{
    mean.ff<-NULL
  }
  
  not.ce<-mean.po[which(div <0 | div >20000),]
  counter<-1
  
  if(dim(not.ce)[1]!=0 & is.null(mean.ff)==F)
  {
    counter<-counter+1
    
    ##second icer eval
    rudi.icer2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff[,2])), as.numeric(as.character(mean.ff[,3])))))*-1
    rudi.icer.ub2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costUB)), as.numeric(as.character(mean.ff$qalyUB)))))*-1
    rudi.icer.lb2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costLB)), as.numeric(as.character(mean.ff$qalyLB)))))*-1
    rudi.icer.sd2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$cost.sd)), as.numeric(as.character(mean.ff$qal.sd)))))*-1
    #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
    
    
    div2<-rudi.icer2[,2]/rudi.icer2[,1]
    div.ub2<-rudi.icer.ub2[,1]/rudi.icer.ub2[,2]
    div.lb2<-rudi.icer.lb2[,1]/rudi.icer.lb2[,2]
    div.sd2<-rudi.icer.sd2[,1]/rudi.icer.sd2[,2]
    
    flip2<-as.matrix(mean.f[which(div2 >0&div2 <20000),])
    
    if(dim(flip2)[1]<2)
    {
      
      mean.f2<-data.frame(t(as.matrix(c(flip2[1:9], 
                                        round(div2[which(div2 >0&div2 <20000)]), 
                                        round(div.ub2[which(div2 >0&div2 <20000)]), #quantile
                                        round(div.lb2[which(div2 >0&div2 <20000)]), #quantile
                                        round(div2[which(div2 >0&div2 <20000)]-1.96*div.sd2[which(div2 >0&div2 <20000)]), round(div2[which(div2 >0&div2 <20000)]+1.96*div.sd2[which(div2 >0&div2 <20000)]))))) #sd
    }else{
      mean.f2<-data.frame(cbind(as.matrix(flip2[,1:9]), 
                                round(div2[which(div2 >0&div2 <20000)]), 
                                round(div.ub2[which(div2 >0&div2 <20000)]), #quantile
                                round(div.lb2[which(div2 >0&div2 <20000)]), #quantile
                                round(div2[which(div2 >0&div2 <20000)]-1.96*div.sd2[which(div2 >0&div2 <20000)]), round(div2[which(div2 >0&div2 <20000)]+1.96*div.sd2[which(div2 >0&div2 <20000)]))) #sd
    }
    
    
    
    if(dim(na.omit(mean.f2))[2]==14)
    {
      colnames(mean.f2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      levels(mean.f2$Intervention)= c(levels(mean.f2$Intervention),'SQ: Risk Groups, 65+ yrs old')
      mean.ff2<-data.frame(as.matrix(mean.f2[order(mean.f2$Rate2, decreasing=F),]))
      colnames(mean.ff2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
    }else{
      mean.ff2<-NULL
    }
    
    not.ce2<-rbind(not.ce, mean.f[which(div2 <0 | div2 >20000),1:9])
    
    
    
    
    if(dim(not.ce2)[1]!=3)
    {
      if(any(dim(not.ce2)!=dim(not.ce)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer3<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2[,2])), as.numeric(as.character(mean.ff2[,3]))))*-1
        rudi.icer.ub3<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costUB)), as.numeric(as.character(mean.ff2$qalyUB))))*-1
        rudi.icer.lb3<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costLB)), as.numeric(as.character(mean.ff2$qalyLB))))*-1
        rudi.icer.sd3<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$cost.sd)), as.numeric(as.character(mean.ff2$qal.sd))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div3<-rudi.icer3[,2]/rudi.icer3[,1]
        div.ub3<-rudi.icer.ub3[,1]/rudi.icer.ub3[,2]
        div.lb3<-rudi.icer.lb3[,1]/rudi.icer.lb3[,2]
        div.sd3<-rudi.icer.sd3[,1]/rudi.icer.sd3[,2]
        
        flip3<-as.matrix(mean.f2[which(div3 >0&div3 <20000),])
        
        #if(dimflip3[1]!=0)
        #{
        mean.f3<-data.frame(t(as.matrix(c(flip3[1:9], 
                                          round(div3[which(div3 >0&div3 <20000)]), 
                                          round(div.ub3[which(div3 >0&div3 <20000)]), #quantile
                                          round(div.lb3[which(div3 >0&div3 <20000)]), #quantile
                                          round(div3[which(div3 >0&div3 <20000)]-1.96*div.sd3[which(div3 >0&div3 <20000)]), round(div3[which(div3 >0&div3 <20000)]+1.96*div.sd3[which(div3 >0&div3 <20000)]))))) #sd
        
        
        #}
        
        if(dim(na.omit(mean.f3))[2]==14)
        { colnames(mean.f3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f3$Intervention)= c(levels(mean.f3$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff3<-data.frame(rbind(sq, as.matrix(mean.f3[order(mean.f2$Rate2, decreasing=F),])))
        colnames(mean.ff3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff3<-NULL}
        not.ce3<-rbind(not.ce2, mean.f2[which(div3 <0 | div3 >20000),1:9])
        
      }
    }
  }
  
  
  ######CE evaluation
  
  if(counter>1)
  {
    mean.ff<-get(paste0('mean.ff', counter))
    not.ce<-get(paste0('not.ce', counter))
  }
  
  
  ######CE evaluation
  #not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  #counter<-1
  #rudi.icer<-diff(rbind(rep(0,6),cbind(rev(mean.po[,3]), rev(mean.po[,4]))))
  #mean.f<-mean.po[which(rudi.icer[,2]/rudi.icer[,1] >0&rudi.icer[,2]/rudi.icer[,1] <20000),]
  #levels(mean.f$Intervention)= c("Preschool"  ,"Primary"  ,  "Secondary", 'Status Quo')
  #mean.ff1<-data.frame(rbind(c(fname[i.country],'Status Quo',as.numeric(0), as.numeric(0)), mean.f[order(as.numeric(as.character(mean.f$Rate2, decreasing=F))),]))
  
  #not.ce1<-mean.po[which(rudi.icer[,2]/rudi.icer[,1] <0 | rudi.icer[,2]/rudi.icer[,1] >20000),]
  
  #if(length(not.ce)!=0 & dim(mean.f)[1]!=0)
  #{
  #counter<-counter+1
  ##second icer eval
  #rudi.icer2<-diff(rbind(c(0,0),cbind(rev(mean.f[,3]), rev(mean.f[,4]))))
  #mean.f2<-mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] >0&rudi.icer2[,2]/rudi.icer2[,1] <20000),]
  #levels(mean.f2$Intervention)= c("Preschool"  , "Primary" , "Secondary"  ,  'Status Quo')
  #mean.ff2<-data.frame(rbind(c(fname[i.country],'Status Quo',as.numeric(0),as.numeric(0)), mean.f2[order(mean.f2$Rate2, decreasing=F),]))
  
  # not.ce2<-rbind(not.ce, mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] <0 | rudi.icer2[,2]/rudi.icer2[,1] >20000),])
  
  
  #if(any(dim(not.ce2)!=dim(not.ce))&length(mean.ff2)!=4&dim(mean.f2)[1]!=0)
  #{
  # counter<-counter+1
  ##3rd CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  #rudi.icer3<-diff(rbind(c(0,0),cbind(rev(mean.f2[,3]), rev(mean.f2[,4]))))
  #mean.f3<-mean.f2[which(rudi.icer3[,2]/rudi.icer3[,1] >0&rudi.icer3[,2]/rudi.icer3[,1] <20000),]
  #levels(mean.f3$Intervention)= c("Preschool"  ,"Primary",  "Secondary"  , 'Status Quo')
  #mean.ff3<-data.frame(rbind(c(fname[i.country], 'Status Quo',as.numeric(0),as.numeric(0)), mean.f3[order(mean.f3$Rate2, decreasing=F),]))
  
  #not.ce3<-rbind(not.ce2, mean.f[which(rudi.icer3[,2]/rudi.icer3[,1] <0 | rudi.icer3[,2]/rudi.icer3[,1] >20000),])
  #}
  #}
  
  
  #mean.ff<-get(paste0('mean.ff', counter))
  #not.ce<-get(paste0('not.ce', counter))
  #colnames(mean.ff)<-c('country','Intervention', 'QALYS', 'Costs')
  #colnames(not.ce)<-c('country','Intervention', 'QALYS', 'Costs')
  
  
  
  ####
  cep<-cep1<-NULL
  col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #llg[llg$cost>1e6|llg$cost < (-1e6),]<-NA
  cep<-ggplot(cost.utility, aes(as.numeric(as.character(QALYS)), as.numeric(as.character(Costs)), color=Intervention))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(breaks=scales::pretty_breaks(n=6), limits = c(-20000,31000))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = 20000, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0(fname[i.country])), fill='Intervention')+theme_bw()+scale_colour_manual(values=c(col[2], col[3], col[4], col[1]))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  #+ylab('Incremental Cost in million ₤/year')+xlab('QALY Difference')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=',strain.name[strain]), paste0('coverage=',cov.f[cov.var])), fill='Intervention')
  
  for(fff in 1:length(inv.names))
  {add.on<-get(paste0('cont',fname[i.country], fff))
  cep<-cep+add.on}
  #cep.name<-paste0('cep',strain, cov.var, dcount.l[dcount])
  
  #assign(cep.name, 
  
  #mean.ff<-data.frame(mean.ff)
  
  if(is.null(mean.ff))
  {
    levels(not.ce$Intervention)=c(levels(not.ce$Intervention),'SQ: Risk Groups, 65+ yrs old')
    not.ce<-data.frame(rbind(sq[1:9], as.matrix(not.ce[order(not.ce$Rate2, decreasing=F),])))
    
    cep1<-cep+geom_point(data=not.ce, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=10, color=c(col[1], col[2], col[4], col[3]), show.legend=F)+theme_bw()
  }else{ if(dim(mean.ff)[1]!=4) 
  {
    colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
    cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_path(data=mean.ff, aes(x=as.numeric(Rate1), y=as.numeric(Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=as.numeric(as.character(not.ce$Rate1)), y=as.numeric(as.character(not.ce$Rate2))), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw()
  }else{
    colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
    cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+theme_bw()
  }
  }
  
  
  #cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_line(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=not.ce$Rate1, not.ce$Rate2), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw()
  
  #)
  
  #assign(cep.name, cep1)
  #cep<-cep1<-NULL
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
  #ggsave(paste0('CEP',strain.name[strain],cov.f[cov.var],dcount,".png"),plot = cep, width=11, height=8.5, units='in',device='png')
  return(cep1)
  
}



all323<-all3232<-NULL
all323<-ggarrange(cep.grapher(4,1,3.5)+rremove("xy.title"),
                  cep.grapher(1,1,3.5)+rremove("xy.title"),
                  cep.grapher(2,1,3.5)+rremove("xy.title"), 
                  cep.grapher(3,1,3.5)+rremove("xy.title"), 
                  cep.grapher(10,1,3.5)+rremove("xy.title"), 
                  cep.grapher(5,1,3.5)+rremove("xy.title"), 
                  cep.grapher(6,1,3.5)+rremove("xy.title"), 
                  cep.grapher(7,1,3.5)+rremove("xy.title"), 
                  cep.grapher(8,1,3.5)+rremove("xy.title"), 
                  cep.grapher(9,1,3.5)+rremove("xy.title"), 
                  common.legend=T, legend='bottom',nrow=5, ncol=2)

all3232<-annotate_figure(all323,
                        bottom = text_grob("QALY Differences"),
                        left = text_grob("Cost in million ₤/year", rot=90)
)

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('CEPallcountries','3.5',version,".png"),plot = all3232, width=8.5, height=11, units='in',device='png')

library(scales)

allh3<-allh32<-NULL
allh3<-ggarrange(cep.grapher(4,2,3.5)+rremove("xy.title"),
                  cep.grapher(1,2,3.5)+rremove("xy.title"),
                  cep.grapher(2,2,3.5)+rremove("xy.title"), 
                  cep.grapher(3,2,3.5)+rremove("xy.title"), 
                  cep.grapher(10,2,3.5)+rremove("xy.title"), 
                  cep.grapher(5,2,3.5)+rremove("xy.title"), 
                  cep.grapher(6,2,3.5)+rremove("xy.title"), 
                  cep.grapher(7,2,3.5)+rremove("xy.title"), 
                  cep.grapher(8,2,3.5)+rremove("xy.title"), 
                  cep.grapher(9,2,3.5)+rremove("xy.title"), 
                  common.legend=T, legend='bottom',nrow=5, ncol=2)

allh32<-annotate_figure(allh3,
                         bottom = text_grob("QALY Differences"),
                         left = text_grob("Cost in million ₤/year", rot=90)
)

setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
ggsave(paste0('CEPallcountriesH3N2','3.5',version,".png"),plot = all3232, width=8.5, height=11, units='in',device='png')


############################Infection averted per additional dose
#X1=infections, 
#X2= country
#X3=intervention

add.dose.l<-do.call(rbind,add.dose)
levels(add.dose.l$X3)<-c('S1: 2-4 year olds', 'S2: 5-11 year olds','S3: 12-16 year olds')
add.dose.l$country<-relevel(add.dose.l$country, "GB")
pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/infecavertedperadddose.samescale.pdf', height=8.5, width =11 )
  print(ggplot(add.dose.l, aes(x=X3, y=(-as.numeric(as.character(X1)))*1000, fill=X3)) +
          geom_violin(draw_quantiles = c(0.5))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+ xlab('Intervention')+ylab('Total Infections Averted per 1000 Additional Vaccine Doses (Ref=SQ)')+facet_wrap(~X2)+labs(fill='Intervention')+scale_fill_manual(values=c(col[2], col[3], col[4]))+theme_bw())
dev.off()

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/infecavertedperadddose.freescale.pdf', height=8.5, width =11 )
print(ggplot(add.dose.l, aes(x=X3, y=(-as.numeric(as.character(X1)))*1000, fill=X3)) +
        geom_violin(draw_quantiles = c(0.5))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+ xlab('Intervention')+ylab('Total Infections Averted per 1000 Additional Vaccine Doses (Ref=SQ)')+facet_wrap(~X2, scales='free')+labs(fill='Intervention')+scale_fill_manual(values=c(col[2], col[3], col[4]))+theme_bw())
dev.off()


#add.ds<-add.dose.l
#colnames(add.ds)<-c('value','country','Strategy')
#levels(add.ds$Strategy)<-c('S1', 'S2', 'S3')
#box.add.doses<-ggplot(add.ds, aes(as.factor(country), 1000*(-as.numeric(as.character(value))), fill=as.factor(country))) +
 # geom_boxplot() +facet_grid(~Strategy)+theme_bw()+coord_flip()+labs(fill='Contact Matrix', y='Infections Averted/1000 additional vaccine doses', x='Contact Matrix')+ggtitle('Infections Averted per 1000 Additional Vaccine Doses')+geom_hline(yintercept = 0, color='dark gray')+scale_y_continuous(breaks=scales::pretty_breaks(n=8))
#ggsave('boxplot.add.doses.pdf',plot=box.add.doses,device = 'pdf')


####################Averted per dose
avt.dose.l<-do.call(rbind,avt.dose)
colnames(avt.dose.l)<-c('value','country','Strategy')
#levels(Hosp.costs.l$Strategy)<-c('S1', 'S2', 'S3', 'SQ')
levels(avt.dose.l$Strategy)<-c('S1: 2-4 year olds', 'S2: 5-11 year olds','S3: 12-16 year olds')
avt.dose.l$country<-factor(avt.dose.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#avt.dose.l$country<-relevel(avt.dose.l$country, "GB")

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/infecaverted.samescale.pdf', height=8.5, width =11 )
print(ggplot(avt.dose.l, aes(x=Strategy, y=(as.numeric(as.character(value)))*1000, fill=Strategy)) +
        geom_violin(draw_quantiles = c(0.5))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+ xlab('Intervention')+ylab('Infections Averted per 1000 Vaccine Doses (Ref=SQ)')+facet_wrap(~country, ncol=5, nrow=2)+labs(fill='Intervention')+theme(legend.position="bottom")+labs(fill='Intervention')+scale_fill_manual(values=c(col[2], col[3], col[4]))+theme_bw())
dev.off()

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/infecaverted.freescale.pdf', height=11, width =8.5 )
print(ggplot(avt.dose.l, aes(x=country, y=(as.numeric(as.character(value)))*1000, fill=country)) +
        geom_violin(scale = 'width', draw_quantiles = c(0.5))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+ xlab('Intervention')+ylab('Infections Averted per 1000 Vaccine Doses (Ref=SQ)')+facet_wrap(~Strategy, scales='free', ncol=1, nrow=3)+labs(fill='Intervention')+theme(legend.position="bottom")+labs(fill='Intervention')+theme_bw()+geom_hline(aes(yintercept=0)))
      
      #+scale_fill_manual(values=c(col[2], col[3], col[4]))+theme_bw())
dev.off()



#################################################################
#Age-stratified barplots
##################################################################
age.strat.loader<-function(strain, i.country, version, vout)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #strategy.summary<-NULL
  age.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'agestrat.',vout,strain.name[strain],'all.programs',version)))
  load(age.tables, envir = .GlobalEnv)
} #ou

#########
#function for epi outcomes
#######

epi.pulls.allcountries<-function(i.vout, i.version, strainpull)
{
  agg.hosp<-c()
for(i.country in 1:(length(fname)-1))
{
age.strat.loader(strain=strainpull,i.country,version=i.version,vout=i.vout)
  
 age.small<-list()
 risk.age<-list()
 prop.small<-list()
 
for(iter in 1:4)
{
  redis<-get(paste0('aged.', i.vout), envir = .GlobalEnv)[[iter]]
  
  #redis<-aged.cases[[iter]]
 
age.small[[iter]]<-cbind(redis[,1]+redis[,2]+redis[,12]+redis[,13], #0-2, 
      redis[,3]+redis[,14], #2-4,
      redis[,4]+redis[,15], #5-11
      redis[,5]+redis[,6]+redis[,16]+redis[,17],
      redis[,7]+redis[,18],
      redis[,8]+redis[,19]+redis[9]+redis[,20]+redis[10]+redis[,21]+redis[11]+redis[,22]#65+
)



#prop.small[[iter]]<-age.small[[iter]]/rowSums(age.small[[iter]])

risk.age[[iter]]<-cbind(redis[,1]+redis[,12],
redis[,2]+redis[,13], 
redis[,3]+redis[,14],
redis[,4]+redis[,15],
redis[,5]+redis[,16],
redis[,6]+redis[,17],
redis[,7]+redis[,18],
redis[,8]+redis[,19],
redis[,9]+redis[,20],
redis[,10]+redis[,21],
redis[,11]+redis[,22]
)

colnames(age.small[[iter]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25+ year olds')

#colnames(prop.small[[iter]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25-44 yr olds', '45-64 yr olds', '65-74 yr olds', '75+ year olds')
}

test<-melt(age.small)
test$L1<-factor(test$L1)
levels(test$L1)<-c('SQ','S1', 'S2','S3')
test2<-cbind(test, rep(fname[i.country],dim(test)[1]))
ifelse(i.country==1, agg.hosp<-test2, agg.hosp<-rbind(agg.hosp, test2))
}
  
  colnames(agg.hosp)<-c('X1','Age','value','Strategy','CM')
  
  hosp.m2<-ddply(agg.hosp[,2:5], .(Age, CM, Strategy), summarise, m1=mean(value), LB=quantile(value, 0.025, na.rm=F), UB=quantile(value, 0.975, na.rm=F))
  hosp.m2$CM<-factor(hosp.m2$CM, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
  hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c('SQ','S1','S2','S3'))
  hosp.m2$Age<-factor(hosp.m2$Age, levels=c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25+ year olds'))
  
  
#graph  
  if(i.vout=='death'| i.vout=='hosp')
  {
    if(i.vout=='hosp')
    {
      epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=CM, group=CM)) + geom_pointrange(aes(ymax=UB, ymin=LB), position=position_dodge(0.85))+theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n Hospitalizations'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')

    }else{
epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=CM, group=CM))+ geom_pointrange(aes(ymax=UB, ymin=LB), position=position_dodge(0.85)) +theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n', paste0(i.vout,'s')), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')
  }}else{
    if(i.vout=='cases') {i.vout='Symptomatic Cases'}
    epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=CM, group=CM))+ 
      geom_pointrange(aes(ymax=UB, ymin=LB), position=position_dodge(0.85)) +theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n', i.vout, '(Thousands)'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')
  }
  
  ##savers
  setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
  ggsave(paste0('epidots',strain.name[strainpull], i.vout,'3.5',i.version,".png"),plot =epio , width=11, height=8.5, units='in',device='png')
  
  
  
  if(i.vout=='death'| i.vout=='hosp')
  {
    if(i.vout=='hosp')
    {
      epio2<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=CM, group=CM)) + geom_point()+theme_bw()+
        geom_line()+facet_wrap(~Age, ncol=3, scales='free_y')+scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n Hospitalizations'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')
      
    }else{
      epio2<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=CM, group=CM))+ geom_point()+
        geom_line() +theme_bw()+facet_wrap(~Age, ncol=3, scales='free_y')+scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n', paste0(i.vout,'s')), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')
    }}else{
      if(i.vout=='cases') {i.vout='Symptomatic Cases'}
      epio2<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=CM, group=CM))+ 
        geom_point()+geom_line()+theme_bw()+facet_wrap(~Age, ncol=3, scales='free_y')+
        scale_color_manual(values=palette, name='Contact Matrix')+labs(y=paste('Average Number of \n', i.vout, '(Thousands)'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'bottom')
    }
  
  
  ##savers
  setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")
  ggsave(paste0('epitrends',strain.name[strainpull], i.vout,'3.5',i.version,".png"),plot =epio2 , width=11, height=8.5, units='in',device='png')
  print(epio2)
  
}


epi.pulls.allcountries('death', 'v1', strainpull = 2)
epi.pulls.allcountries('cases', 'v1', strainpull = 2)
epi.pulls.allcountries('hosp', 'v1', strainpull = 2)
epi.pulls.allcountries('GP', 'v1', strainpull = 2)




#agg.hosp$CM<-relevel(agg.hosp$CM, 'GB')

##Violin plots
pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/agestrat.deaths.freescale.pdf', height=8.5, width =11 )
for(jj in 1:10)
{
  print(ggplot(agg.hosp,aes(x=CM,y=as.numeric(as.character(value)), fill=CM))+geom_violin(draw_quantiles = c(0.5))+facet_wrap_paginate(Age~Strategy, ncol=2, nrow = 2, page=jj)+scale_fill_manual(values=palette)+theme_bw())
}
dev.off()


####Plot of Means
#colnames(agg.hosp)<-c('X1','Age','value','Strategy','CM')
#hosp.m2<-ddply(agg.hosp[,2:5], .(Age, CM, Strategy), summarise, m1=mean(value), LB=quantile(value, 0.025, na.rm=F), UB=quantile(value, 0.975, na.rm=F))
 #hosp.m2$CM<-factor(hosp.m2$CM, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
#hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c('SQ','S1','S2','S3'))

#ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=CM, group=CM)) + geom_point()+theme_bw()+facet_wrap(~Age, ncol=3, scales='free_y')+geom_line()+scale_color_manual(values=palette, name='Contact Matrix')+labs(y='Average Number of \nSymptomatic Cases (Thousands)', x='Intervention Strategy')+theme_bw()


#####################

#geom_errorbar(aes(ymin =LB , ymax = UB), width = 0.2)+

q<-ggplot(inv.costs.f, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var],',', 'strain=', strain.name[strain]), x = 'Intervention Strategy', y = 'Average Annual Costs(£)', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
  axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")


#pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/cum.incidence.freescale.pdf', height=11, width =8.5)
print(ggplot(cumi.l, aes(as.factor(country), (as.numeric(as.character(value)))/1e6, fill=as.factor(country)))+geom_violin(draw_quantiles = c(0.5))+scale_y_continuous(breaks=scales::pretty_breaks(n=8))+facet_wrap(~Strategy, scale='free', ncol=2, nrow=2)+theme_bw()+labs(fill='Intervention', y='Average Annual Cumulative Incidence (Millions)', x='Intervention')+theme(legend.position="bottom"))





#################Cumulative Incidence
cumi.l<-do.call(rbind.data.frame, cumi.l)
colnames(cumi.l)<-c('value','country','Strategy')
levels(cumi.l$Strategy)<-c('S1: 2-4 year olds', 'S2: 5-11 year olds','S3: 12-16 year olds', 'SQ: High Risk Groups Only')
#cumi.l$Strategy<-relevel(cumi.l$Strategy, 'SQ: High Risk Groups Only')
cumi.l$country<-factor(cumi.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(cumi.l$country)

cumi.sum<-ddply(cumi.l, .(country, Strategy), summarise, m1=mean(round(as.numeric(as.character(value))))/1e6,
                LB=quantile(round(as.numeric(as.character(value))), 0.025)/1e6, UB=quantile(round(as.numeric(as.character(value)))/1e6, 0.975))
setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens")

#pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/cum.incidence.freescale.pdf', height=11, width =8.5)
cumi<-ggplot(cumi.l, aes(as.factor(country), (as.numeric(as.character(value)))/1e6, fill=as.factor(country)))+geom_boxplot()+scale_y_continuous(breaks=scales::pretty_breaks(n=6))+facet_wrap(~Strategy, scale='free_y', ncol=2, nrow=2)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Cumulative Incidence (Millions)', x='Contact Matrix')+theme(legend.position="bottom")+scale_fill_manual(values = palette)
#+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))
dev.off()
ggsave(paste0('tot.infection',version,".png"),plot =cumi , width=11, height=8.5, units='in',device='png')






############Hospital Costs
sum.costs.l<-do.call(rbind.data.frame, sum.costs.l)
colnames(sum.costs.l)<-c('value','country','Strategy')
levels(sum.costs.l$Strategy)<-c('S1: 2-4 year olds', 'S2: 5-11 year olds','S3: 12-16 year olds', 'SQ: High Risk Groups Only')
sum.costs.l$country<-factor(sum.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(sum.costs.l$country)
sum.costs<-ddply(sum.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))), UB=quantile(as.numeric(as.character(value)), 0.975), LB=quantile(as.numeric(as.character(value)), 0.025))

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/total.cost.freescale.pdf', height=11, width =8.5 )
print(ggplot(sum.costs, aes(x=as.factor(country), y=(m1/1e6), fill=as.factor(country)))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=8))+facet_wrap(~Strategy, scale='free', ncol=2, nrow=2)+theme_bw()+labs(fill='Intervention', y='Average Annual Costs from Influenza (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(Strategy), ymin=m1/1e6, ymax=UB/1e6), width=.2,  position='dodge', inherit.aes = F))
#+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))
dev.off()


GP.costs.l<-do.call(rbind.data.frame, GP.costs.l)
colnames(GP.costs.l)<-c('value','country','Strategy')
levels(GP.costs.l$Strategy)<-c('S1', 'S2', 'S3', 'SQ')
GP.costs.l$country<-factor(GP.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(GP.costs.l$country)
GP.costs<-ddply(GP.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))), UB=quantile(as.numeric(as.character(value)), 0.975), LB=quantile(as.numeric(as.character(value)), 0.025))

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/GP.cost.freescale.pdf', height=11, width =8.5 )
print(ggplot(GP.costs, aes(x=as.factor(Strategy), y=(m1/1e6), fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=8))+facet_wrap(~country, scale='free', ncol=4, nrow=3)+theme_bw()+labs(fill='Intervention', y='Average Annual GP Costs from Influenza (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(Strategy), ymin=m1/1e6, ymax=UB/1e6), width=.2,  position='dodge', inherit.aes = F))
dev.off()



cum.l<-do.call(rbind.data.frame, cum.l)
colnames(Hosp.costs.l)<-c('value','country','Strategy')
levels(Hosp.costs.l$Strategy)<-c('S1', 'S2', 'S3', 'SQ')
Hosp.costs.l$country<-factor(Hosp.costs.l$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels(Hosp.costs.l$country)
hosp.sum<-ddply(Hosp.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))), UB=quantile(as.numeric(as.character(value)), 0.975), LB=quantile(as.numeric(as.character(value)), 0.025))

pdf('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/hosp.outcomes.freescale.pdf', height=8.5, width =11 )
print(ggplot(hosp.sum, aes(x=as.factor(Strategy), y=(m1/1e6), fill=as.factor(Strategy)))+scale_fill_manual(values=c(col[2], col[3], col[4], col[1]))+geom_bar(stat='identity', color="black", position=position_dodge())+scale_y_continuous(breaks=scales::pretty_breaks(n=8))+facet_wrap(~country, scale='free', ncol=4, nrow=3)+theme_bw()+labs(fill='Contact Matrix', y='Average Annual Hospitalization Costs (£ Millions)', x='Intervention')+geom_errorbar(aes(x=as.factor(Strategy), ymin=m1/1e6, ymax=UB/1e6), width=.2,  position='dodge', inherit.aes = F))
dev.off()


#################Table DATA#################
lapply(1:14, colnames(Hosp.costs.l)<-c('value','country', 'Strategy'))

hosp.sum<-ddply(Hosp.costs.l, .(country, Strategy), summarise, m1=round(median(as.numeric(as.character(value)))),LB=round(quantile(as.numeric(as.character(value)), 0.025)), UB=round(quantile(as.numeric(as.character(value)), 0.975)))
avt.sum<-ddply(avt.dose.l, .(X2, X3), summarise, m1=mean(as.numeric(as.character(X1))),LB=quantile(as.numeric(as.character(X1)), 0.025), UB=quantile(as.numeric(as.character(X1)), 0.975))
cum.sum<-ddply(cumi.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value))),LB=quantile(as.numeric(as.character(value)), 0.025), UB=quantile(as.numeric(as.character(value)), 0.975))
costs.sum<-ddply(sum.costs.l, .(country, Strategy), summarise, m1=mean(as.numeric(as.character(value)))/1e6,LB=quantile(as.numeric(as.character(value)), 0.025)/1e6, UB=quantile(as.numeric(as.character(value)), 0.975)/1e6)
add.sum<-ddply(add.dose.l, .(X2, X3), summarise, m1=mean(as.numeric(as.character(X1))),LB=quantile(as.numeric(as.character(X1)), 0.025), UB=quantile(as.numeric(as.character(X1)), 0.975))

QALY.tot.cum<-ddply(QALY.tot.l, .(country, Strategy), summarise, m1=round(mean(as.numeric(as.character(value)))),LB=round(quantile(as.numeric(as.character(value)), 0.025)), UB=round(quantile(as.numeric(as.character(value)), 0.975)))
QALY.death.cum<-ddply(QALY.death.l, .(country, Strategy), summarise, m1=round(mean(as.numeric(as.character(value)))),LB=round(quantile(as.numeric(as.character(value)), 0.025)), UB=round(quantile(as.numeric(as.character(value)), 0.975)))
GP.cum<-ddply(GP.costs.l, .(country, Strategy), summarise, m1=round(mean(as.numeric(as.character(value)))),LB=round(quantile(as.numeric(as.character(value)), 0.025)), UB=round(quantile(as.numeric(as.character(value)), 0.975)))



mean(Reduce('+',lapply(vaccine.cost.statquo[1:14], rowSums))/tab.diff)/1e6
quantile(Reduce('+',lapply(vaccine.cost.statquo[1:14], rowSums))/tab.diff, c(0.025, 0.975))/1e6

mean(Reduce('+',lapply(vaccine.cost2[1:14], rowSums))/tab.diff)/1e6
quantile(Reduce('+',lapply(vaccine.cost2[1:14], rowSums))/tab.diff, c(0.025, 0.975))/1e6

mean(Reduce('+',lapply(vaccine.cost3[1:14], rowSums))/tab.diff)/1e6
quantile(Reduce('+',lapply(vaccine.cost3[1:14], rowSums))/tab.diff, c(0.025, 0.975))/1e6

mean(Reduce('+',lapply(vaccine.cost4[1:14], rowSums))/tab.diff)/1e6
quantile(Reduce('+',lapply(vaccine.cost4[1:14], rowSums))/tab.diff, c(0.025, 0.975))/1e6




##########
barplot<-function(strain, cov.var, dcount)
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', cov.f[cov.var]))
  
  #make status quo category
  for(tt in 1:2)
  {
    discount.loader(strain,cov.f[cov.var], dcount)
    inv.costs<-get(paste0('inv.costs',strain.name[strain], dcount.l[dcount]))
    
    inv.costs$variable <- factor(inv.costs$variable, levels=c(paste0('vaccine cost',tt), paste0('Healthcare costs',tt), paste0('GP costs',tt), paste0('hospital costs',tt),paste0('QALYs lost',tt),paste0('QALYs lost case',tt), paste0('QALYs lost death',tt), paste0('QALYs lost hosp', tt)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost (All Sources)', 'Years of Life Lost (Symptomatic Cases)', 'Years of Life Lost (Deaths)','Years of Life Lost (Hospitalization)'))
    if(tt==1) inv.costs$strategy<-'Status Quo'
    ifelse(tt==1, inv.costs.f<-na.omit(inv.costs)[1:6,], inv.costs.f<-rbind(inv.costs.f, na.omit(inv.costs)))
  }
  #cols<-my.col
  
  q<-ggplot(inv.costs.f, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var],',', 'strain=', strain.name[strain]), x = 'Intervention Strategy', y = 'Average Annual Costs(£)', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
    axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")
  
  
  ggsave(paste0('costsplot',strain.name[strain],cov.f[cov.var],dcount.l[dcount],".png"),plot = q, width=11, height=7.5, units='in',device='png')
  return(q)
}

barplot(1, cov.var=2, dcount=3)
  

##############################################






intl.ICER.loader<-function(strain, i.country, dcount, intervent)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  strategy.summary<-NULL
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.intl', strain.name[strain],'p', '*', 'd',dcount)))
  
  cer<<-load(cea.tables[[intervent]])
} #outputs as ICER.set2 list, QALY diff is difference from status Quo

incidence<-function(i.country, program)
{
cov.sens.loader(strain=1,i.country=i.country, end.cov=0.55, version='v1')

dataset<-unlist(lapply(intl.tableH1N1[[program]], rowSums))
yet<-c(fname[i.country],inv.names2[program], mean(dataset), quantile(dataset, c(0.025, 0.975)))
return(yet)
}

incidence(2,1)


ddply(Hosp.costs.l, .(country, Strategy), summarise, m1=round(median(as.numeric(as.character(value)))),LB=round(quantile(as.numeric(as.character(value)), 0.025)), UB=round(quantile(as.numeric(as.character(value)), 0.975)))

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')


  intl.ICER.loader(1,1,3.5,2)

  load(cea.tables[[2]])
inv.names<-c('Preschool', 'Primary', 'Secondary')
cos.f<-qds.f<-qds.m<-cos.m<-NULL
strain<-1


library(scales)
llg<-na.omit(llg)
pdf(file='cepintvstratbyage1.pdf',onefile = T, height=8.5, width=11)
for(j in 1:17)
{
   print(ggplot(llg, aes(as.numeric(as.character(qaly)), cost, color=program))+ theme_grey()+
      geom_point(alpha=.5)+stat_ellipse(type='norm', alpha = .9, aes(color = program), level=0.90)+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+scale_x_continuous(breaks=scales::pretty_breaks(n=8))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, color='blue')+geom_abline(slope = 20000, intercept=0, color='light blue')+facet_wrap_paginate(Age~program, ncol = 1, nrow = 3, page=j, scales='free')+xlab('Contact Matrix')+ylab('Cost per 1 QALY Gained')+ labs(fill='Intervention'))
  
}
dev.off()



####total CEP

pdf(file='cep.total.countrybyprogram.pdf',onefile = T, height=8.5, width=11)
for(j in 1:17)
{
  print(ggplot(llg[llg$Age=='total',], aes(as.numeric(as.character(qaly)), cost, color=program))+ theme_grey()+
          geom_point(alpha=.5)+stat_ellipse(type='norm', alpha = .9, aes(color = program), level=0.90)+scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+scale_x_continuous(breaks=scales::pretty_breaks(n=8))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, color='blue')+geom_abline(slope = 20000, intercept=0, color='light blue')+facet_wrap_paginate(program~country, ncol = 2, nrow = 3, page=j, scales='free')+xlab('Contact Matrix')+ylab('Cost per 1 QALY')+ labs(fill='Intervention')+ scale_fill_viridis('viridis'))
  
}
dev.off()




#######################################
#Probability curves/ Acceptability curves
########################################

##########################################
####Acceptability curve, probability each intervention has optimal cost-effectiveness
#################################################
  
  inv.names2<-c('SQ','I1', 'I2', 'I3')
  
  setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
  
  accept.curve<-function(threshold.val, qaly, cost)
  {
      qaly2<- Map('*', -1, Map('*',qaly,threshold.val))
      nmb<-Map('-',qaly2, cost)
      #nmb3<-Reduce('+', nmb)/tab.diff
    
    #curve.data.strat<-c(as.character(fname[i.country]),as.character(inv.names2[intv]), threshold.val, as.numeric(colMeans(nmb3)), as.numeric(mean(rowMeans(nmb3))))
    
    return(nmb)
  }

#accept.curve(20000,tot.QALY.loss[[1]], tot.costs.loss[[1]])
  
  
  
c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
  
outerf<-list()
accept.out1<-accept.out2<-accept.out3<-accept.out4<-list()
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  
  
  load(qal.tables[[1]])
  load(cos.tables[[1]])
  n.samples<-dim(tot.QALY.loss[[1]][[1]])[1]
  tab.diff<-length(tot.QALY.loss[[1]])
  ss<-seq(0,30000,by = 300)
  for(hh in 1:length(ss)) {
    
    
    accept.out1<-accept.curve(ss[hh], tot.QALY.loss[[1]], tot.costs.loss[[1]]) 
    accept.out2<-accept.curve(ss[hh],tot.QALY.loss[[2]], tot.costs.loss[[2]])
    accept.out3<-accept.curve(ss[hh],tot.QALY.loss[[3]], tot.costs.loss[[3]])
    accept.out4<-accept.curve(ss[hh],tot.QALY.loss[[4]], tot.costs.loss[[4]])
    
    pat1<-data.frame(as.factor(rep(inv.names2[1],n.samples)), Reduce('+',accept.out1)/tab.diff, rowSums(Reduce('+',accept.out1)/tab.diff))
    pat2<-data.frame(as.factor(rep(inv.names2[2],n.samples)), Reduce('+',accept.out2)/tab.diff, rowSums(Reduce('+',accept.out2)/tab.diff))
    pat3<-data.frame(as.factor(rep(inv.names2[3],n.samples)), Reduce('+',accept.out3)/tab.diff, rowSums(Reduce('+',accept.out3)/tab.diff))
    pat4<-data.frame(as.factor(rep(inv.names2[4],n.samples)), Reduce('+',accept.out4)/tab.diff, rowSums(Reduce('+',accept.out4)/tab.diff))
    
    colnames(pat1)<-c.names
    
    
    colnames(pat2)<-colnames(pat1)
    colnames(pat3)<-colnames(pat1)
    colnames(pat4)<-colnames(pat1)
    
    #am.1[[hh]]<-pat1
    #am.2[[hh]]<-pat2
    #am.3[[hh]]<-pat3
    #am.4[[hh]]<-pat4
    
    
    for(aa in 2:dim(pat1)[2]){
      e.frame<-data.frame(pat1[,aa], pat2[,aa],pat3[,aa],pat4[,aa])
      xx<-apply(e.frame, 1, function(z) {which.max(z)})  
      ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=4)
      for(cc in 1:length(xx))
      {ex.frame[cc,xx[cc]]<-1}
      
      
      outp<-data.frame(fname[i.country], c.names[aa], ss[hh], t(colMeans(ex.frame)))
      
      
      ifelse(aa==2&hh==1&i.country==1, outerf<-outp, outerf<-rbind(outerf, outp))
  }
    
  }
  }


setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
colnames(outerf)<-c('country', 'Age','threshold','SQ','I1','I2','I3')


dev.off()

outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
#m<-ggplot(data=outer.g, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=variable, linetype=variable))  
 # pdf(file='OC.3program.pdf',onefile = T, height=11, width=8.75)
  #for(gg in 1:32){
   # print(m+ facet_wrap_paginate(country~Age, ncol = 2, nrow = 4, scales='free', page=gg)+
 #           geom_line(aes(color=variable, group=variable, linetype=variable))+labs(title='Probability a strategy is Optimal', x='Willing-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
            #scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_v#line(xintercept = c(15000,20000), linetype=3))
 # }
  #dev.off()
  
  
  ##################totals only
  
  outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
  mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=variable))  
  #pdf(file='OC.totals.3program.pdf',onefile = T, height=11, width=8.75)
optim.mk<-mk+facet_wrap(~country, ncol = 2, nrow = 5, scales='free_y')+geom_line(aes(color=variable, group=variable))+labs(x='Willingness-to-pay Threshold (GBP per QALY)', y='Probability that strategy is Optimal',color='Strategy') +
            scale_y_continuous(breaks=scales::pretty_breaks(n=7))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"))+theme_bw()+theme(legend.position="bottom")

 
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
ggsave(paste0('OCall',dcount,".png"),plot = optim.mk, width=8.5, height = 11, units='in',device='png')

  
###GB only
  
  
  gb.grob<-outer.g[outer.g$Age=='total'&outer.g$country=='GB',]
           
  col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
  oc.gb<-ggplot(data= gb.grob, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), group=variable,color=variable))+scale_color_manual(values=c(col))+
    geom_line()+labs(title='Optimal Acceptability Curve \n per Strategy', x='Willingness-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
    scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+ theme_bw()+
    annotate("text", x = 29000, y = 0.05, label = paste0('GB'))
  
  
  dev.off()
  


##########################################################################################################
# ACCEPTIBILITY CURVE for single intervention probability it is cost effective on WTP scale
#################################################################################################################################

  prob.curve<-function(threshold.val, int, i.country)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.intl', strain.name[strain],'p', '*', 'd',dcount)))
  
  
  load(cea.tables[[int-1]]);
  ICER.set.n<-get(paste0('ICER.set', int))
  
  qal2<-unlist(ICER.set.n$`net QALY`)*threshold.val
  nmb<-Map('-',qal2, unlist(ICER.set.n$`net cost`))
  #nmb2<-as.matrix(Reduce('+',nmb)/tab.diff)
  
  #nmb2<-Reduce('+', ICER.set2$`icer`)/tab.diff

  ex.zeros<-matrix(rep(0,length(nmb)), ncol=1)
  
  ex.zeros[nmb > 0]<-1  
  
  curve.data.strat<-c(as.character(fname[i.country]), threshold.val, as.numeric(colMeans(ex.zeros)))
  
  return(curve.data.strat)
}

outer<-list()
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
  ss<-seq(1000,30000,by = 1000)
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out2<-prob.curve(ss[hh], 2, i.country ), 
                                  accept.out2<-rbind(accept.out2, prob.curve(ss[hh], 2, i.country )))}
  
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out3<-prob.curve(ss[hh], 3, i.country ), 
                                  accept.out3<-rbind(accept.out3, prob.curve(ss[hh], 3, i.country )))}
  
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out4<-prob.curve(ss[hh], 4, i.country ), 
                                  accept.out4<-rbind(accept.out4, prob.curve(ss[hh], 4, i.country )))}

  
  #colnames(accept.out2)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
  
  #colnames(accept.out2)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')

  colnames(accept.out2)<-c('country','WTP','Probability')
colnames(accept.out2)<-colnames(accept.out3)
colnames(accept.out2)<-colnames(accept.out4)

outer[[i.country]]<-list(accept.out2, accept.out3, accept.out4)
}

names(outer)<-fname[1:10]

#install.packages("ggthemes") # Install 


for(gg in 1:length(outer))
{
  for(hh in 1:length(outer[[4]]))
  {
dtm<-outer[[gg]][[hh]]

dtm<-data.frame(dtm)
#colnames(dtm)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
colnames(dtm)<-c('country','WTP','var')

#accept.m<-melt(dtm,id.vars = c('country', 'threshold'))

accept.m<-dtm
program<-rep(inv.names2[hh+1], dim(accept.m)[1])

accept.r<-cbind(accept.m, program)

#ifelse(hh==1, accept.out<-accept.r, accept.out<-rbind(accept.out, accept.r))
ifelse(gg==1&hh==1, accept.out<-accept.r, accept.out<-rbind(accept.out, accept.r))
  }
  }

####
#pdf(file='AC.total.3program.pdf',onefile = T, height=11, width=8.75)
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')

prob.m<-ggplot(data=accept.out, aes(x=as.numeric(as.character(WTP)), y=as.numeric(as.character(var)), color=program))+facet_wrap(~country, ncol = 2, nrow = 5, scales='free')+scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
          geom_line(aes(color=program, group=program))+labs(title='Acceptability Curve per Strategy', x='Willingness-to-pay (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=6))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)


ggsave(paste0('ACall',dcount, version,".png"),plot = prob.m, width=8.5, height = 11, units='in',device='png')


##strat by country
accept.out$country<-factor(accept.out$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))

prob.nl<-ggplot(data=accept.out, aes(x=as.numeric(as.character(WTP)), y=as.numeric(as.character(var)), color=country))+facet_wrap(~program, nrow = 3, scales='free_y')+scale_color_manual(values=palette)+
  geom_line(aes(color=country, group=country))+labs(title='Acceptability Curve per Country', x='Willingness-to-pay (GB£ per QALY)', y='Probability Cost-Effective',color='Contact Matrix') +
  scale_y_continuous(breaks=scales::pretty_breaks(n=6))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)

ggsave(paste0('ACpercountry',dcount, version,".png"),plot = prob.nl, width=8.5, height = 11, units='in',device='png')


#######

ac.gb<-ggplot(data=accept.out[accept.out$country=='GB',], aes(x=as.numeric(as.character(WTP)), y=as.numeric(as.character(var)), color=program, group=program))+scale_color_manual(values=col[2:4])+
  geom_line()+labs(title='Acceptability Curve \n per Strategy', x='Willingness-to-pay (£)', y='Probability Cost-Effective',color='Strategy') +
  scale_y_continuous(breaks=scales::pretty_breaks(n=7))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme_bw()+theme(legend.position="none")+
  annotate("text", x = 29000, y = 0.05, label = paste0('GB'))


########Optimal strategy##################################

optimal.strat.intl<-function(i.country, strain, collect)
{
  qal.all<-list()
  cost.all<-list()
  outerf<-list()
  
  qal.tables<-cos.tables<-c()

  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  
  load(qal.tables[[1]])
  load(cos.tables[[1]])
  n.samples<-dim(tot.QALY.loss[[1]][[1]])[1]
  tab.diff<-length(tot.QALY.loss[[1]])
  
  #for(ii in 1:4) 
  #{
    #optimal strategy compares all strategy outcomes to each other and sums to 1
    #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[ii]))
   # 
    #if(ii==1){
     # qal.all[[ii]]<-tot.QALY.loss[[1]];
      #cost.all[[ii]]<-tot.costs.loss[[1]];
  #  }
   # qal.all[[ii+1]]<-tot.QALY.loss[[2]]
    #cost.all[[ii+1]]<-tot.costs.loss[[2]]
  #}
   # }
  #Reduce('+',lapply(tot.QALY.loss[[oo]], rowSums))/tab.diff
  intv.number<-length(tot.QALY.loss)
  qal.m<-lapply(lapply(1:length(tot.QALY.loss), function(oo) Reduce('+', tot.QALY.loss[[oo]])/tab.diff), rowSums)
  cost.m<-lapply(lapply(1:length(tot.QALY.loss), function(oo) Reduce('+', tot.costs.loss[[oo]])/tab.diff), rowSums)
  
  ss<-seq(1,40000,by = 500)
  for(hh in 1:length(ss)) {
    
    #cov.name3<-c('0%','30%','55%','70%')
    #AC.all<-lapply(1:length(qal.all), function(oo) accept.curve(ss[hh], qal.all[[oo]], cost.all[[oo]]) )
    #AC.all2<-lapply(1:length(qal.all), function(pp) data.frame(rep(cov.name2[pp],n.samples), Reduce('+',AC.all[[pp]])/tab.diff, rowSums(Reduce('+',AC.all[[pp]])/tab.diff)))
    
    ttest0<-(-qal.m[[1]]+0)-(-cost.m[[1]]-0)/ss[hh]
    ttesti<-list()
    #QALYS gained-additional cost/WTP
    for(oo in 2:intv.number) {ttesti[[oo]]<-(-qal.m[[oo]]+qal.m[[1]])-(-cost.m[[oo]]+cost.m[[1]])/ss[hh]}
    #for(oo in 2:intv.number) {ttesti[[oo]]<-(qal.m)-(cost.m)/ss[hh]}
    
    OC.all<-list(ttest0, ttesti[[2]], ttesti[[3]], ttesti[[4]])
    
    #OC.all2<- lapply(OC.all, rowSums)
    #n.samples<-length(OC.all[[1]])
    OC.all3<-lapply(lapply(1:length(OC.all), function(pp) data.frame(rep(as.character(inv.names2[pp],n.samples)), OC.all[[pp]])), unname)
    
    
    new.ac<-array(unlist(OC.all3), dim=c(n.samples, dim(OC.all3[[1]])[2], intv.number))
    
    #for(aa in 2:24){
    e.frame<-data.frame(new.ac[,2,1:intv.number])
    xx<-apply(e.frame, 1, function(z) {which.max(z)})  
    ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=intv.number)
    for(cc in 1:length(xx))
    {ex.frame[cc,xx[cc]]<-1}
    outp<-data.frame(fname[i.country], ss[hh], t(colMeans(ex.frame)))
    
    ifelse(hh==1, outerf<-outp, outerf<-rbind(outerf, outp))
    #}
  }
  
  colnames(outerf)<-c('country', 'threshold',inv.names2)
  #colnames(outerf)<-c('country','threshold','0%','30%','55%','70%')
  
  #outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
  outer.g<-melt(outerf, id.vars=c('country', 'threshold'))
  
  #colnames(outer.g)<-c('country', 'Age', 'Threshold','Coverage','value')
  colnames(outer.g)<-c('country', 'Threshold','Intervention','value')
  #mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=Coverage, linetype=Coverage))
  if(collect==T){ outer.all[[i.country]]<<-outer.g}
     
  
 if(collect==F)
   {mk<-ggplot(data= outer.g, aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value))))

    #final graph
  OC3<-mk+geom_line(aes(color=Intervention))+scale_color_manual(values = col)+scale_y_continuous(breaks=scales::pretty_breaks(n=6), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Coverage"))+theme_bw()+labs(title=fname[i.country])
  
  #OC3<-OC3+scale_color_manual(values=cols[1:8])
  return(OC3)
 }
}

#arrange output by coverage and discount
optimal3<-NULL
optimal3<-ggarrange(
  optimal.strat.intl(4,1, F)+rremove("xy.title"), optimal.strat.intl(1,1, F)+rremove("xy.title"), optimal.strat.intl(2,1,F)+rremove("xy.title"), 
  optimal.strat.intl(3,1, F)+rremove("xy.title"), optimal.strat.intl(10,1, F)+rremove("xy.title"), optimal.strat.intl(5,1,F)+rremove("xy.title"),
  optimal.strat.intl(6,1,F)+rremove("xy.title"), optimal.strat.intl(7,1,F)+rremove("xy.title"), optimal.strat.intl(8,1, F)+rremove("xy.title"),optimal.strat.intl(9,1,F)+rremove("xy.title"),
  ncol=2, nrow=5, common.legend = T, legend='bottom')

graphs.ocpost2<-NULL
graphs.ocpost2<-annotate_figure(optimal3,
                                bottom = text_grob('Willingness-to-pay Threshold (GBP per QALY)'),
                                left = text_grob('Probability that Strategy is Optimal', rot=90))

ggsave(paste0("OC3.all.intl.",version,".png"),plot=graphs.ocpost2, width=8.5, height=11, units='in',device='png')


outer.all<-list()
for(i.country in 1:length(fname))
{
  optimal.strat.intl(i.country,1,collect=T)
  
}


allopts<-do.call(rbind.data.frame, outer.all)

a.fronter<-ddply(allopts,.(country, Threshold, Intervention), summarize, mm=max(value)) 

mk4<-ggplot(data=allopts, aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), group=country))+facet_wrap(~Intervention, nrow=5, ncol=2)

#final graph
OC4<-mk4+geom_line(aes(color=country))+scale_color_manual(values =palette)+scale_y_continuous(breaks=scales::pretty_breaks(n=6), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Coverage"))+theme_bw()



######
doSummary = function(df) {
  mm=max(df$value)
  Intervention = df$Intervention
  data.frame(mm=mm, Intervention=Intervention)
}

a.fronter2<-ddply(allopts,.(country, Threshold), doSummary) 
a.fronter2<-ddply(allopts,.(country, Threshold),summarise, mm=max(value)) 

a.fronter3<-merge(a.fronter2, allopts, by.x = c('country', 'Threshold', 'mm'), by.y = c('country','Threshold', 'value'))


a.fronter3$country<-factor(a.fronter3$country, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))

mk5<-ggplot(data=a.fronter3, aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(mm))))+geom_line(aes(color=country))+scale_color_manual(values =palette)+geom_line(aes(color=country))+scale_y_continuous(breaks=scales::pretty_breaks(n=10), limits = c(0.3,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000))+theme_bw()+labs(y='Pr(Intervention is Optimal)', x='Willingness-to-pay \nThreshold (GBP per QALY)')

#+geom_path(data=allopts, aes(as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=country, group=country), linetype='dashed')

mk5<-mk5+geom_point(aes(shape=Intervention, color=country))+scale_shape_manual(values=c(3, 16, 17))+theme(legend.position="bottom") 

ggsave(paste0("CEACcc.",version,".png"),plot=mk5, width=11, height=8.5, units='in',device='png')

                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                            ###############################################################                                                                                              
#final graph
OC4<-mk4+geom_line(aes(color=country))+scale_y_continuous(breaks=scales::pretty_breaks(n=6), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Coverage"))+theme_bw()



optimal.strat.intl<-function(i.country, strain, collect)
{
  qal.all<-list()
  cost.all<-list()
  outerf<-list()
  
  qal.tables<-cos.tables<-c()
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  
  load(qal.tables[[1]])
  load(cos.tables[[1]])
  n.samples<-dim(tot.QALY.loss[[1]][[1]])[1]
  tab.diff<-length(tot.QALY.loss[[1]])
  
  #for(ii in 1:4) 
  #{
  #optimal strategy compares all strategy outcomes to each other and sums to 1
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[ii]))
  # 
  #if(ii==1){
  # qal.all[[ii]]<-tot.QALY.loss[[1]];
  #cost.all[[ii]]<-tot.costs.loss[[1]];
  #  }
  # qal.all[[ii+1]]<-tot.QALY.loss[[2]]
  #cost.all[[ii+1]]<-tot.costs.loss[[2]]
  #}
  # }
  #Reduce('+',lapply(tot.QALY.loss[[oo]], rowSums))/tab.diff
  intv.number<-length(tot.QALY.loss)
  qal.m<-lapply(lapply(1:length(tot.QALY.loss), function(oo) Reduce('+', tot.QALY.loss[[oo]])/tab.diff), rowSums)
  cost.m<-lapply(lapply(1:length(tot.QALY.loss), function(oo) Reduce('+', tot.costs.loss[[oo]])/tab.diff), rowSums)
  
  ss<-seq(1,40000,by = 1000)
  for(hh in 1:length(ss)) {
    
    #cov.name3<-c('0%','30%','55%','70%')
    #AC.all<-lapply(1:length(qal.all), function(oo) accept.curve(ss[hh], qal.all[[oo]], cost.all[[oo]]) )
    #AC.all2<-lapply(1:length(qal.all), function(pp) data.frame(rep(cov.name2[pp],n.samples), Reduce('+',AC.all[[pp]])/tab.diff, rowSums(Reduce('+',AC.all[[pp]])/tab.diff)))
    
    ttest0<-(-qal.m[[1]]+0)-(-cost.m[[1]]-0)/ss[hh]
    ttesti<-list()
    #QALYS gained-additional cost/WTP
    for(oo in 2:intv.number) {ttesti[[oo]]<-(-qal.m[[oo]]+qal.m[[1]])-(-cost.m[[oo]]+cost.m[[1]])/ss[hh]}
    #for(oo in 2:intv.number) {ttesti[[oo]]<-(qal.m)-(cost.m)/ss[hh]}
    
    OC.all<-list(ttest0, ttesti[[2]], ttesti[[3]], ttesti[[4]])
    
    #OC.all2<- lapply(OC.all, rowSums)
    #n.samples<-length(OC.all[[1]])
    OC.all3<-lapply(lapply(1:length(OC.all), function(pp) data.frame(rep(as.character(inv.names2[pp],n.samples)), OC.all[[pp]])), unname)
    
    
    new.ac<-array(unlist(OC.all3), dim=c(n.samples, dim(OC.all3[[1]])[2], intv.number))
    
    #for(aa in 2:24){
    e.frame<-data.frame(new.ac[,2,1:intv.number])
    xx<-apply(e.frame, 1, function(z) {which.max(z)})  
    ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=intv.number)
    for(cc in 1:length(xx))
    {ex.frame[cc,xx[cc]]<-1}
    outp<-data.frame(fname[i.country], ss[hh], t(colMeans(ex.frame)))
    
    ifelse(hh==1, outerf<-outp, outerf<-rbind(outerf, outp))
    #}
  }
  
  colnames(outerf)<-c('country', 'threshold',inv.names2)
  #colnames(outerf)<-c('country','threshold','0%','30%','55%','70%')
  
  #outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
  outer.g<-melt(outerf, id.vars=c('country', 'threshold'))
  
  #colnames(outer.g)<-c('country', 'Age', 'Threshold','Coverage','value')
  colnames(outer.g)<-c('country', 'Threshold','Intervention','value')
  #mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=Coverage, linetype=Coverage))
  if(collect==T){ outer.all[[i.country]]<<-outer.g}
  
  
  mk<-ggplot(data= outer.g, aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value))))
  
  #final graph
  OC3<-mk+geom_line(aes(color=Intervention))+scale_color_manual(values = col)+scale_y_continuous(breaks=scales::pretty_breaks(n=6), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Coverage"))+theme_bw()+labs(title=fname[i.country])
  
  #OC3<-OC3+scale_color_manual(values=cols[1:8])
  return(OC3)
}













outerf<-list()
accept.out1<-accept.out2<-accept.out3<-accept.out4<-list()
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  
  
  load(qal.tables[[1]])
  load(cos.tables[[1]])
  n.samples<-dim(tot.QALY.loss[[1]][[1]])[1]
  tab.diff<-length(tot.QALY.loss[[1]])
  ss<-seq(0,30000,by = 300)
  for(hh in 1:length(ss)) {
    
    
    accept.out1<-accept.curve(ss[hh], tot.QALY.loss[[1]], tot.costs.loss[[1]]) 
    accept.out2<-accept.curve(ss[hh],tot.QALY.loss[[2]], tot.costs.loss[[2]])
    accept.out3<-accept.curve(ss[hh],tot.QALY.loss[[3]], tot.costs.loss[[3]])
    accept.out4<-accept.curve(ss[hh],tot.QALY.loss[[4]], tot.costs.loss[[4]])
    
    pat1<-data.frame(as.factor(rep(inv.names2[1],n.samples)), Reduce('+',accept.out1)/tab.diff, rowSums(Reduce('+',accept.out1)/tab.diff))
    pat2<-data.frame(as.factor(rep(inv.names2[2],n.samples)), Reduce('+',accept.out2)/tab.diff, rowSums(Reduce('+',accept.out2)/tab.diff))
    pat3<-data.frame(as.factor(rep(inv.names2[3],n.samples)), Reduce('+',accept.out3)/tab.diff, rowSums(Reduce('+',accept.out3)/tab.diff))
    pat4<-data.frame(as.factor(rep(inv.names2[4],n.samples)), Reduce('+',accept.out4)/tab.diff, rowSums(Reduce('+',accept.out4)/tab.diff))
    
    colnames(pat1)<-c.names
    
    
    colnames(pat2)<-colnames(pat1)
    colnames(pat3)<-colnames(pat1)
    colnames(pat4)<-colnames(pat1)
    
    #am.1[[hh]]<-pat1
    #am.2[[hh]]<-pat2
    #am.3[[hh]]<-pat3
    #am.4[[hh]]<-pat4
    
    
    for(aa in 2:dim(pat1)[2]){
      e.frame<-data.frame(pat1[,aa], pat2[,aa],pat3[,aa],pat4[,aa])
      xx<-apply(e.frame, 1, function(z) {which.max(z)})  
      ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=4)
      for(cc in 1:length(xx))
      {ex.frame[cc,xx[cc]]<-1}
      
      
      outp<-data.frame(fname[i.country], c.names[aa], ss[hh], t(colMeans(ex.frame)))
      
      
      ifelse(aa==2&hh==1&i.country==1, outerf<-outp, outerf<-rbind(outerf, outp))
    }
    
  }
}

