#R0 and Re calculations
local({r <- getOption("repos")
      r["CRAN"] <- "http://cran.r-project.org"
     options(repos=r)
})
.libPaths('/home/nwenzel/intl')

install.packages("ggplot2")
install.packages("scales")
install.packages("gridExtra")
install.packages("ggthemes")
install.packages("ggforce")
install.packages("abind")
install.packages("cowplot")
install.packages("RColorBrewer")
install.packages("R.utils")
install.packages("randomcolor")
install.packages("plotly")



library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(ggforce)       # kable : prettier data.frame output
library(abind)
library(R.utils)
library(cowplot)
library(vioplot)
require(RColorBrewer)
require(randomcoloR)
require(MASS)

library(plyr)
library(reshape)
library(data.table) 
library(fluEvidenceSynthesis)
library(plotly)
# MOST IMPORTANT PACKAGE SO THIS CODE ISN'T SLOW
#"/Users/Natasha/Library/R/3.3/library" 
#library(xtable)
library(dplyr)
library(tidyr)
#.libPaths("/Users/Natasha/Library/R/3.3/library" )
#devtools::install_github('richfitz/dde', upgrade = FALSE)
#devtools::install_github('mrc-ide/cinterpolate', upgrade = FALSE)
#devtools::install_github('mrc-ide/odin', upgrade = T)
rm(list = ls())
library(odin)
###############################################################################
#Load in input data
################################################################################
##################################################
####### FILE LOADING FUNCTIONS FOR GRAPHS ########
#################################################
#loaders for all files
setwd('~/intl')
source('SUPPORT_R0_loaders.R')
source('Output_vacc_program_simulation.R')
source('SUPPORT_ICER_loaders.R')
vstrategy<-dget('FUNC_cov_strategy_uk.R')

#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/OUTPUT_ICER_warp_sampler.R')
load('Intl_general_parameters.RData')

##load general parameters---------------------
p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving

strain.name<-c('H1N1','H3N2','B')

load('ili.counts.rda')
load('virological.rda')

#load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/ili2009.rda')

fiftyfive<-c(coverageB[[18]]$V12[4:23],0.5500)
fiftyfive2<-c(fiftyfive, rep(0.55,3))
fiftyfive2[fiftyfive2>0.55]<-0.55
fiftyfive.short<-fiftyfive2[1:24][seq(1,dim(coverageB[[18]])[1],4)]
fiftyfive.cov<-list(fiftyfive.short, fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short, fiftyfive2, fiftyfive2, fiftyfive2, fiftyfive2, fiftyfive2)


####INTL FUNCTION Start=================================================
num.samp<-1000
new.cov<-fiftyfive.cov
r.complete<-vector('list',length(fname))

strain<-1
for(i.country in 1:length(fname))
{
  #Sampler for R0 and Re----------------------
  
  vacc.program.simulation(program = 1, i.country=i.country, strain.choice = strain, Reff=T, version='u1')
  vacc.program.simulation(program = 2, i.country=i.country, strain.choice = strain,Reff=T,version='u1')
  vacc.program.simulation(program = 3, i.country=i.country, strain.choice = strain, Reff=T, version='u1')
  vacc.program.simulation(program = 4, i.country=i.country, strain.choice = strain, Reff=T, version='u1')
  
  vacc.program.simulation(program = 1, i.country=i.country, strain.choice = strain, Reff=F,version='u1')
  vacc.program.simulation(program = 2, i.country=i.country, strain.choice = strain, Reff=F,version='u1')
  vacc.program.simulation(program = 3, i.country=i.country, strain.choice = strain, Reff=F,version='u1')
  vacc.program.simulation(program = 4, i.country=i.country, strain.choice = strain, Reff=F,version='u1')
}#end country loop







###################Attack Rate and R0 sensitivity to input parameters

cov.sens.loader<-function(strain, i.country, end.cov, version)
{
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
  #Rearrange tables to correct order of interventions------------------
  sq.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStratStatusQuo', end.cov, strain.name[strain], version)))
  sens.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStrat*','*chool', end.cov, strain.name[strain], version)))
  
  if(length(sens.tables)==0){
    sens.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStrat*','*chool', end.cov, strain.name[strain])));
  }
  if(length(sq.tables)==0)
  {sq.tables<-list.files(pattern=glob2rx(paste0(fname[i.country],'SensStratStatusQuo', end.cov, strain.name[strain])));
  print(paste(fname[i.country], 'base version used'))}  
  
  
  load(sens.tables[1]); avg.incidence1<-total.keep; #preschool
  load(sens.tables[2]); avg.incidence2<-total.keep; #primary
  load(sens.tables[3]); avg.incidence3<-total.keep; #secondary
  load(sq.tables[1]); avg.incidence4<-total.keep; #status quo
  
  
  #Reduce('+',total.keep)/length(total.keep); 
  
  #"UKAgeStratPreschool+Primary SchoolB"    "UKAgeStratPreschool+Primary+SecondaryB"
  #[3] "UKAgeStratPreschool+SecondaryB"         "UKAgeStratPreschoolB"                  
  #[5] "UKAgeStratPrimary+SecondaryB"           "UKAgeStratPrimarySchoolB"              
  #[7] "UKAgeStratSecondarySchoolB"             "UKAgeStratStatusQuoB" 
  #vaccines.all<-list(vaccines.statquo,vaccines.preschool,vaccines.prime,vaccines.sec,vaccines.preprimary,vaccines.preprimesec,vaccines.preprimesec,vaccines.presec,vaccines.primesec)
  var.name<-paste0('intl.table',strain.name[strain])
  
  assign(var.name, list(avg.incidence4,avg.incidence1,avg.incidence2, avg.incidence3), envir = .GlobalEnv)
  #save(tableH1,file=paste0('AnnualIncH1'))
}



########################Graphs--------------

strain.name<-c('H1N1','H3N2','B');
p.name2<-c('Preschool','PrimarySchool', 'SecondarySchool', "StatusQuo"  )   

y.name<<-c(1995:2014)
epidemic.outcomes<-list()
epidemic.means<-list()


strain<-i.strain<-1
for(i.country in 1:length(fname))
{
  #if(i.country==4) next
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  #allocate memory
  R0.tablesB<-list.files(pattern=glob2rx(paste0(fname[i.country],'R0*B')));
  R0.tablesH1<-list.files(pattern=glob2rx(paste0(fname[i.country],'R0*H1N1')));
  R0.tablesH3<-list.files(pattern=glob2rx(paste0(fname[i.country],'R0*H3N2')));
  
  RE.tablesB<-list.files(pattern=glob2rx(paste0(fname[i.country],'RE*B')));
  RE.tablesH1<-list.files(pattern=glob2rx(paste0(fname[i.country],'RE*H1N1')));
  RE.tablesH3<-list.files(pattern=glob2rx(paste0(fname[i.country],'RE*H3N2')));
  
  A.tablesB<-list.files(pattern=glob2rx(paste0(fname[i.country],'Assort*B')));
  A.tablesH1<-list.files(pattern=glob2rx(paste0(fname[i.country],'Assort*H1N1')));
  A.tablesH3<-list.files(pattern=glob2rx(paste0(fname[i.country],'Assort*H3N2')));
  
  Tr.tablesB<-list.files(pattern=glob2rx(paste0(fname[i.country],'transmissibility*B')));
  Tr.tablesH1<-list.files(pattern=glob2rx(paste0(fname[i.country],'transmissibility*H1N1')));
  Tr.tablesH3<-list.files(pattern=glob2rx(paste0(fname[i.country],'transmissibility*H3N2')));
  
  p.name2<-c('Preschool','PrimarySchool', 'SecondarySchool', "StatusQuo"  )   
  graphics.off()
  
  #program<-1
  for(program in 1:4)
  {
    #if(program==1|program==3|program==4|program==5|program==6) next
    
    #load in data for specific programs
    load(R0.tablesH1[program]); r0.programh1<-total.keep;
    #load(R0.tablesH3[program]); r0.programh3<-total.keep;
    #load(R0.tablesB[program]); r0.programb<-total.keep;
    
    
    load(RE.tablesH1[program]); re.programh1<-total.keep;
    #load(RE.tablesH3[program]); re.programh3<-total.keep;
    #load(RE.tablesB[program]); re.programb<-total.keep;
    
    
    load(A.tablesH1[program]); A.programh1<-total.akeep;
    #load(A.tablesH3[program]); A.programh3<-total.akeep;
    #load(A.tablesB[program]); A.programb<-total.akeep;
    
    
    #arrange
    #r0.allprogs<-data.frame(cbind(r0.programh1[[1]],r0.programh3[[1]], r0.programb[[1]]))
    #colnames(r0.allprogs)<-strain.name
    
    #re.allprogs<-data.frame(cbind(re.programh1[[1]],re.programh3[[1]], re.programb[[1]]))
    #colnames(re.allprogs)<-strain.name
    
    #A.allprogs<-data.frame(cbind(A.programh1[[1]],A.programh3[[1]], A.programb[[1]]))
    #colnames(A.allprogs)<-strain.name
    
    r0.allprogs<-data.frame(r0.programh1[[1]])
    colnames(r0.allprogs)<-strain.name[1]
    
    re.allprogs<-data.frame(re.programh1[[1]])
    colnames(re.allprogs)<-strain.name[1]
    
    A.allprogs<-data.frame(A.programh1[[1]])
    colnames(A.allprogs)<-strain.name[1]
    
    
    r0.frame<-rE.frame<-rA.frame<-NULL
    r0.list<-rE.list<-A.list<-list()
    
    #i.season<-4
    
    #for(i.strain in 1:length(strain.name))
    #{
    if(i.strain==1) {r0.program<-r0.programh1; rE.program<-re.programh1; A.program<-A.programh1}
    if(i.strain==2) {r0.program<-r0.programh3; rE.program<-re.programh3; A.program<-A.programh3}
    if(i.strain==3) {r0.program<-r0.programb;  rE.program<-re.programb ; A.program<-A.programb}
    
    #final size
    cov.sens.loader(strain,i.country,end.cov=0.55, version='u1')
    fs.pull<-intl.tableH1N1[[program]]
    fs.sums<-lapply(fs.pull, rowSums)
    
    
    for (i.season in 1:length(total.keep))
    {
      r0.short<-cbind(r0.program[[i.season]], as.numeric(r0.program[[i.season]][,2])>=1)
      
      #rE.short.e<-matrix(as.numeric(rE.program[[i.season]][as.numeric(rE.program[[i.season]][,2])>=1,]), ncol=2)
      rE.short<-cbind(rE.program[[i.season]], as.numeric(rE.program[[i.season]][,2])>=1)
      A.short<- cbind(A.program[[i.season]], as.numeric(r0.program[[i.season]][,2])>=1)
      
      year.filler<-function(df) {matrix(as.numeric(rep(y.name[i.season], dim(df)[1]), ncol=1))}
      
      #convert to data frame for ggplot2
      
      nf.r0.e<-cbind(r0.short, A.short[,1], rE.short[,2], fs.sums[[i.season]], year.filler(r0.short))
      nf.rE.e<-cbind(rE.short, year.filler(rE.short))
      nf.A.e<-cbind(A.short, year.filler(A.short))
      
      #if(length(nf.r0.e)==0) nf.r0.e<-c(0,0,0,y.name[i.season])
      #if(length(nf.rE.e)==0) nf.rE.e<-c(0,0,y.name[i.season])
      
      #colnames(nf)<-strain.name; #add variable names
      ifelse(i.season==1, r0.frame<-nf.r0.e, r0.frame<-rbind(r0.frame, nf.r0.e));
      ifelse(i.season==1, rE.frame<-nf.rE.e, rE.frame<-rbind(rE.frame, nf.rE.e));
      ifelse(i.season==1, A.frame<-nf.A.e,   A.frame<-rbind(A.frame, nf.A.e));
      
    }
    
    colnames(r0.frame)<-c('timestep','R0','epidemic','Assortativity','RE', 'Final Size','Season')
    colnames(rE.frame)<-c('timestep','RE','epidemic','Season')
    colnames(A.frame)<-c('Assortativity','epidemic','Season')
    r0.list[[i.strain]]<-r0.frame
    rE.list[[i.strain]]<-rE.frame
    A.list[[i.strain]]<-A.frame
    #}
    
    #names(r0.list)<-c('H1N1','H3N2','B')
    #names(rE.list)<-c('H1N1','H3N2','B')
    #names(A.list)<-c('H1N1','H3N2','B')
    
    names(r0.list)<-c('H1N1')
    names(rE.list)<-c('H1N1')
    names(A.list)<-c('H1N1')
    
    
    
    rearrange<-function(list, piece) {yy<-data.frame(list[[piece]], rep(strain.name[piece],length(list[[piece]]))); 
    setnames(yy, names(yy[dim(yy)[2]]), "strain")
    return(yy)}
    
    #r0.df<-rbind(rearrange(r0.list,1), rearrange(r0.list,2), rearrange(r0.list,3))
    #rE.df<-rbind(rearrange(rE.list,1), rearrange(rE.list,2), rearrange(rE.list,3))
    #A.df<-rbind(rearrange(A.list,1), rearrange(A.list,2), rearrange(A.list,3))
    
    r0.df<-rbind(rearrange(r0.list,1))
    rE.df<-rbind(rearrange(rE.list,1))
    A.df<-rbind(rearrange(A.list,1))
    
    #add a bind for each country for a massive data frame for the tornado type plots
    c.repeat<-rep(as.factor(fname[i.country]), dim(r0.df)[1])
    p.repeat<-rep(as.factor(p.name2[program]), dim(r0.df)[1])
    ifelse(program==1, r0.output<-cbind(c.repeat,p.repeat, r0.df), r0.output<-rbind(r0.output, cbind(c.repeat,p.repeat, r0.df)))
    c.repeat<-rep(as.factor(fname[i.country]), dim(rE.df)[1])
    p.repeat<-rep(as.factor( p.name2[program]), dim(rE.df)[1])
    ifelse(program==1, rE.output<-cbind(c.repeat,p.repeat, rE.df), rE.output<-rbind(rE.output, cbind(c.repeat,p.repeat, rE.df)))
    c.repeat<-rep(as.factor(fname[i.country]), dim(A.df)[1])
    p.repeat<-rep(as.factor( p.name2[program]), dim(A.df)[1])
    ifelse(program==1, A.output<-cbind(c.repeat,p.repeat, A.df), A.output<-rbind(A.output, cbind(c.repeat,p.repeat, A.df)))
    
  } 
  
  #epidemic years only
  ifelse(i.country==1|program==1, r0.outcomes<-r0.output, r0.outcomes<-rbind(r0.output, r0.outcomes))
  ifelse(i.country==1|program==1, rE.outcomes<-rE.output, rE.outcomes<-rbind(rE.output, rE.outcomes))
  ifelse(i.country==1|program==1, A.outcomes<-A.output, A.outcomes<-rbind(A.output, A.outcomes))
}

# Run the functions length, mean, and sd on the value of "change" for each group, 
r0.cdata <- ddply(r0.outcomes, c('c.repeat','p.repeat',"Season",'strain'), summarise,
                  N    = length(R0),
                  mean1 = mean(R0),
                  mean2= mean(Assortativity),
                  sd2   = sd(R0),
                  sd22 = sd(Assortativity),
                  se1   = sd2 / sqrt(N),
                  se2 = sd22/sqrt(N),
                  me1=abs(se1*1.96),
                  me2=abs(se2*1.96))


r0.cdata <- ddply(r0.outcomes, c('c.repeat','p.repeat',"Season",'strain'), summarise,
                  N    = length(R0),
                  mean1 = mean(R0),
                  mean2= mean(Assortativity),
                  sd2   = sd(R0),
                  sd22 = sd(Assortativity),
                  se1   = sd2 / sqrt(N),
                  se2 = sd22/sqrt(N),
                  me1=abs(se1*1.96),
                  me2=abs(se2*1.96))

A.cdata <- ddply(r0.outcomes, c('c.repeat',"Season",'strain'), summarise,
                 N    = length(R0),
                 mean2= mean(Assortativity),
                 sd2 = sd(Assortativity),
                 se1 = sd2/sqrt(N),
                 me1=abs(se1*1.96))





rE.cdata <- ddply(rE.outcomes, c('c.repeat','p.repeat',"Season", "strain"), summarise,
                  N    = length(RE),
                  mean = mean(RE),
                  sd2   = sd(RE),
                  se   = sd2 / sqrt(N),
                  me=abs(se*1.96))

rA.cdata <- ddply(A.df, c('c.repeat','p.repeat', "Season", "strain"), summarise,
                  N    = length(Assortativity),
                  mean = mean(Assortativity),
                  sd2   = sd(Assortativity),
                  se   = sd2 / sqrt(N),
                  me=abs(se*1.96))

epidemic.means[[i.country]]<-list(r0.cdata, rE.cdata, rA.cdata)
names(epidemic.means[[i.country]])<-c('R0 & A','RE','Assortativity')


FS.data <- ddply(r0.outcomes, c('c.repeat','p.repeat','Final.Size','strain'), summarise,
                 N    = length(Final.Size),
                 mean1 = mean(Final.Size),
                 sd1 = sd(Final.Size),
                 me1=abs(sd1*1.96))




#####Tornado of countries and R0, RE, Assort

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'IntlSens'))
pdf(file='R0.Assort.bycountry.pdf',onefile = T, height=8.5, width=11)
for(j in 1:4)
{
  print(ggplot(r0.outcomes[r0.outcomes$R0>0,], aes(y=R0, x=Assortativity, color=p.repeat)) + geom_point()+theme_grey()+stat_ellipse(type='norm', alpha = .9, aes(color = p.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(breaks=scales::pretty_breaks(n=7) )+geom_vline(aes(xintercept=1))+facet_wrap_paginate(p.repeat~c.repeat, ncol = 11, page=j, scales='free')+ylab('Reproductive Number')+ labs(fill='Intervention'))
}
dev.off()



pdf(file='R0.Assort.allcountries.pdf',onefile = T, height=8.5, width=11)
print(ggplot(r0.outcomes, aes(y=R0, x=Assortativity, color=c.repeat)) + geom_point(alpha=.85)+stat_ellipse(type='norm', alpha = .9, aes(color = c.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(limits=c(0.7,1),breaks=scales::pretty_breaks(n=5))+geom_hline(aes(yintercept=1))+ylab('Reproductive Number')+ labs(fill='Intervention'))
dev.off()


#Assortativity vs R0 graphs all countries together
mj<-ggplot(r0.outcomes, aes(y=R0, x=Assortativity, color=c.repeat))
pdf(file='R0.Assort.allcountries.pdf',onefile = T, height=8.5, width=11)
print( mj+ geom_point(alpha=.85)+stat_ellipse(type='norm', alpha = .9, aes(color = c.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(limits=c(0.7,1),breaks=scales::pretty_breaks(n=5))+geom_hline(aes(yintercept=1))+ylab('Reproductive Number')+ labs(fill='Intervention'))
dev.off()





#Assortativity vs R0 graphs faceted by country 
mj<-ggplot(r0.outcomes, aes(y=R0, x=Assortativity, color=c.repeat))
pdf(file='R0.Assort.bycountry.pdf',onefile = T, height=8.5, width=11)
for(gg in 1:3){
  print( mj+facet_wrap_paginate(~c.repeat, ncol = 2, nrow = 2, scales='free', page=gg)+geom_point(alpha=.85)+stat_ellipse(type='norm', alpha = .9, aes(color = c.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(limits=c(0.7,1),breaks=scales::pretty_breaks(n=5))+geom_hline(aes(yintercept=1))+ylab('Reproductive Number')+ labs(title='Assortativity vs R)', fill='Intervention'))}
dev.off()


#Assortativity vs RE graphs faceted by country 
mj<-ggplot(r0.outcomes, aes(y=RE, x=Assortativity, color=c.repeat))
pdf(file='RE.Assort.bycountry.pdf',onefile = T, height=8.5, width=11)
for(gg in 1:3){
  print( mj+facet_wrap_paginate(~c.repeat, ncol = 2, nrow = 2, scales='free', page=gg)+geom_point(alpha=.85)+stat_ellipse(type='norm', alpha = .9, aes(color = c.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(limits=c(0.7,1),breaks=scales::pretty_breaks(n=5))+geom_hline(aes(yintercept=1))+ylab('Reproductive Number')+ labs(title=expression('Assortativity vs R'[e]), fill='Intervention'))}
dev.off()


#Assortativity vs RE graphs faceted by country 
mj<-ggplot(r0.outcomes, aes(y=RE, x=Assortativity, color=c.repeat))
pdf(file='RE.Assort.bycountryintervention.pdf',onefile = T, height=8.5, width=11)
for(gg in 1:13){
  print( mj+facet_wrap_paginate(~c.repeat, ncol = 2, nrow = 2, scales='free', page=gg)+geom_point(alpha=.85)+stat_ellipse(type='norm', alpha = .9, aes(color = c.repeat))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(limits=c(0.7,1),breaks=scales::pretty_breaks(n=5))+geom_hline(aes(yintercept=1))+ylab('Reproductive Number')+ labs(title=expression('Assortativity vs R'[e]), fill='Intervention'))}
dev.off()

#Box plots of R0, Re, and Associativity

mh<-ggplot(r0.outcomes, aes(x=c.repeat, y=R0, fill=c.repeat)) +
  geom_boxplot()+labs(fill='Contact Matrix')
mh1<-mh+xlab('Country Contact Matrix')+ylab(expression('R'[0]))+geom_hline(aes(yintercept=1))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+theme(legend.position="none")

mh2<-ggplot(r0.outcomes, aes(x=c.repeat, y=RE, fill=c.repeat)) +
  geom_boxplot()+labs(fill='Contact Matrix')+xlab('Country Contact Matrix')+ylab(expression('R'[e]))+geom_hline(aes(yintercept=1))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+theme(legend.position="none")

mh3<-ggplot(r0.outcomes, aes(x=c.repeat, y=Assortativity, fill=c.repeat)) +
  geom_boxplot()+scale_y_continuous(limits=c(0.7,1))+xlab('Country Contact Matrix')+ylab('Assortativity')+theme(legend.position="none")

mh4<-ggplot(r0.outcomes, aes(x=c.repeat, y=Final.Size, fill=c.repeat)) +
  geom_boxplot()+xlab('Country Contact Matrix')+ylab('Cumulative Incidence')+theme(legend.position="none")

pdf(file='Repro.output.fs.pdf',onefile = T, height=12, width=11) #for all interventions
print(plot_grid(mh1, mh2, mh3, mh4, ncol=1, nrow=4))
dev.off()

plot_grid(mh1, mh2, mh3, mh4, ncol=1, nrow=4, labels = 'AUTO', align = 'h', label_size = 12)
ggsave(paste0('Repro.output.fs.pdf','.tiff'),device='tiff',plot=last_plot(), width=30,height=23,units='cm')



###Same Graph as above but stratified by interventions. Which interventions limit Re, are some more effective depending on Assortativity.


pdf(file='Repro.output.byintv.pdf',onefile = T, height=8.5, width=11) #for all interventions
mh<-ggplot(r0.outcomes, aes(x=c.repeat, y=R0, fill=c.repeat)) +
  geom_boxplot()+labs(fill='Contact Matrix')
pdf(file='cep.total.overlap.pdf',onefile = T, height=8.5, width=11)
mh<-ggplot(r0.outcomes, aes(x=c.repeat, y=R0, fill=c.repeat)) +
  geom_boxplot()
mh1<-mh+xlab('Country Contact Matrix')+ylab(expression('R'[0]))+geom_hline(aes(yintercept=1))+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+theme(legend.position="none")+ ggtitle(expression('R'[0]), 'for all interventions')

mh2<-ggplot(r0.outcomes, aes(x=c.repeat, y=RE, fill=c.repeat)) 
pdf(file='Re.boxplot.byintv.pdf',onefile = T, height=8.5, width=11)
for(j in 1:3)
{
  print(mh2+geom_boxplot()+labs(fill='Contact Matrix')+xlab('Country Contact Matrix')+ylab(expression('R'[e]))+geom_hline(aes(yintercept=1))+scale_y_continuous(limits=c(0,2.5), breaks=scales::pretty_breaks(n = 8))+theme(legend.position="none")+facet_wrap_paginate(~p.repeat, ncol = 1, nrow = 2, page=j, scales='free'))
}
dev.off()

mh3<-ggplot(r0.outcomes, aes(x=c.repeat, y=Assortativity, fill=c.repeat)) +
  geom_boxplot()+scale_y_continuous(limits=c(0.7,1))+xlab('Country Contact Matrix')+ylab('Assortativity')+theme(legend.position="none")

pdf(file='Re.boxplot.byintv.grouped.pdf',onefile = T, height=8.5, width=11)
print(plot_grid(mh1, mh2, mh3, ncol=1, nrow=3))

dev.off()
