library(fluEvidenceSynthesis)
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(data.table)  # MOST IMPORTANT PACKAGE SO THIS CODE ISN'T SLOW
library(abind)
library(xtable)
library(plyr)
library(reshape2)
library(ggplot2)
library(cowplot)

rm(list = ls())

alp<-list.files(path='/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata', pattern=glob2rx('*table11.csv'))
allpolymod<-c(alp[1], alp[2],alp[3],alp[5],alp[6], alp[7], alp[8],alp[10], alp[9],alp[4],alp[11])
fname<-c('BE', 'DE', 'FI', 'GB', 'IT', 'LU', 'NL', 'PL', 'PE', 'FR', 'ZI')
y.name<<-c(1995:2014)
age.group.limits<<-c(1,2,5,12,15,17,25,45,65,75) #upper limits

risk.ratios.ce<<-matrix(c(0.021,0.055,0.055,0.098,0.098,0.098,0.087,0.092,0.183,0.45,0.45,0,0,0,0,0,0,0,0,0,0,0),ncol=(length(age.group.limits)+1), byrow=TRUE)  

age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)
strain.name<-c('H1N1','H3N2','B');
####multiplot function load-----------------------------------------------------
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
}


####Start Function------------------------------------------------------------------
post.means<-array()
prior.means<-array()
strain<-1
country<-4
i.country<-country
version<-'v1'

for(strain in 1:length(strain.name))
{
#for (i.country in 1:(length(fname)))
#{
polymod<-fread(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Inputdata', paste(allpolymod[i.country])))
polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds

#Reconstruct contract matricies from posterior contact data \\
#if(i.country==4) 
#{setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
 # flu.tablesB<-list.files(pattern=glob2rx('UKflu.p8*B'));
  #flu.tablesH1<-list.files(pattern=glob2rx('UKflu.p8*H1'));
  #flu.tablesH3<-list.files(pattern=glob2rx('UKflu.p8*H3'));
  
  #ct.tablesB<-list.files(pattern=glob2rx('UKct.p8*B'))
  #ct.tablesH1<-list.files(pattern=glob2rx('UKct.p8*H1'))
  #ct.tablesH3<-list.files(pattern=glob2rx('UKct.p8*H3'))
#} else  {
#if(country<9) {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version,'*H1', sep='')))
  ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version,'*H3', sep='')))
  ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version,'*B', sep='')))
  
  flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.', version,'*B',sep='')));
  flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version, '*H1',sep='')));
  flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version,'*H3',sep='')));
#}

#if(country>8) {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
 # ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.k2*H1', sep='')))
  #ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.k2*H3', sep='')))
  #ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.k2*B', sep='')))
  
 # flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.k2*B',sep='')));
  #flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.k2*H1',sep='')));
  #flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.k2*H3',sep='')));
}

  
  
  SENS.post.pars<-function(strain, i.country, num.samp, version)
  {
    mcmc.samp<-function(i.season, total.samp)
    {
      #attach(jj)
      #load GB mcmcbatch, and other country ct.table
      setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
      if(strain==1) {load(ct.tablesH1[i.season])};
      if(strain==2) {load(ct.tablesH3[i.season])};
      if(strain==3) {load(ct.tablesB[i.season])};
      
      if(strain==1) {load(flu.tablesH1[i.season])}; 
      if(strain==2) {load(flu.tablesH3[i.season])}; 
      if(strain==3) {load(flu.tablesB[i.season])}; 
      
      mcmc.post<-mcmc.result$batch #rename mcmc output to simple data.table
      
      iteration<-c(round(runif(total.samp,1,dim(mcmc.post)[1]*0.05))) 
      #take only last 20% of the runs to use for posterior distribution for calculations
      
      post.sample[[i.season]]<<-mcmc.post[iteration,] 
      rand.contact.samp[[i.season]]<<-contact.ids[iteration]
      #remove(jj)
      
    }
    post.sample<-list()
    rand.contact.samp<-list()
    total.samp<-num.samp
    #remove(jj)
    #if(i.country==4)
    #{setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', '0.55'))
     # flu.tablesB<-list.files(pattern=glob2rx('GBflu.u5*B'));
      #flu.tablesH1<-list.files(pattern=glob2rx('GBflu.u5*H1'));
      #flu.tablesH3<-list.files(pattern=glob2rx('GBflu.u5*H3'));
      
      #ct.tablesB<-list.files(pattern=glob2rx('GBct.u5*B'))
      #ct.tablesH1<-list.files(pattern=glob2rx('GBct.u5*H1'))
      #ct.tablesH3<-list.files(pattern=glob2rx('GBct.u5*H3'))
    #}else
    #{
      setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
      ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version, '*H1', sep='')))
      ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version, '*H3', sep='')))
      ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version, '*B', sep='')))
      
      flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.',version, '*B',sep='')));
      flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.',version, '*H1',sep='')));
      flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version, '*H3',sep='')));
    #}
    
    if(strain==1){year.length<-length(ct.tablesH1)}
    if(strain==2){year.length<-length(ct.tablesH3)}
    if(strain==3){year.length<-length(ct.tablesB)}
    
    #jj<-where('year.length')
    sapply(1:year.length, FUN= function(x) mcmc.samp(x, num.samp))
    #rm(i.country)
    p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving
    
    ###Save output to file
    kname<-paste0(fname[i.country],'PriorSample',version, strain.name[strain])  
    save(post.sample,file=kname) 
    jname<-paste0(fname[i.country],'PriorContactSample',version,strain.name[strain])  
    save(rand.contact.samp, file=jname)
  }
  
 
  for(i.country in 1:length(fname)){ SENS.prior.pars(strain=1, i.country, num.samp=300, version='v1')}
  
  
  
  
  
  post.ct.sampler<-function(program, i.country, strain.choice, version)
  {
    
    if(i.country==4)
    {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
      ct.samp<-list.files(pattern=glob2rx(paste0('UK','PostContactSample',version, strain.name[strain.choice])))
      
      }else{
      setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
      ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PostContactSample',version, strain.name[strain.choice])))
      
    }
    
    
    if(length(ct.samp)==0)
    {
      stop(paste(fname[i.country], 'base version used'))
    }
    
    load(ct.samp);
    
    polymod<<-NULL
    polymod<-as.data.frame(fread(file.path('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata',paste(allpolymod[i.country]))))
    polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
    
    year.length<-length(rand.contact.samp);
    post.means<-list(); post.sd<-list()
    
    for(i.season in 1:year.length)
    {
      if(strain==1) {contact.set<-rand.contact.samp[[i.season]]; }
      if(strain==2) {contact.set<-rand.contact.samp[[i.season]]; }
      
      p.contacts<-NULL
      
      for (jj in 1:length(contact.set))
      {
        p.contacts[[jj]]<-contact_matrix(as.matrix(polymod[contact.set[[jj]],]),
                                         age_sizes[,1], age.group.limits)
      }
      post.means[[i.season]]<-Reduce('+', p.contacts)/length(p.contacts)
      post.sd[[i.season]]<- sqrt(Reduce('+', lapply(Map('-',  p.contacts, post.means[[i.season]]), function(x) x^2))/length(p.contacts))
    }
    
    #version<-'v1'
    save(post.means,file=paste0('post.means',fname[i.country], version, strain.name[strain]))
    save(post.sd,file=paste0('post.sd',fname[i.country], version, strain.name[strain]))
    
  }



  
  
for(i.country in 1:length(fname)){post.ct.sampler(1, i.country, strain.choice = 1, version='v1')}
  
  

#for(i.country in 4:length(fname))
#{if(i.country==4) {post.ct.sampler(1, i.country, strain.choice = 1, version='u5')
 # }else{
#  post.ct.sampler(1, i.country, strain.choice = 1, version='u3')}
#}



######################################
#######PRIOR samples###################
######################################


SENS.prior.pars<-function(strain, i.country, num.samp, version)
{
  mcmc.samp<-function(i.season, total.samp)
  {
    #attach(jj)
    #load GB mcmcbatch, and other country ct.table
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
    if(strain==1) {load(ct.tablesH1[i.season])};
    if(strain==2) {load(ct.tablesH3[i.season])};
    if(strain==3) {load(ct.tablesB[i.season])};
    
    if(strain==1) {load(flu.tablesH1[i.season])}; 
    if(strain==2) {load(flu.tablesH3[i.season])}; 
    if(strain==3) {load(flu.tablesB[i.season])}; 
    
    mcmc.post<-mcmc.result$batch #rename mcmc output to simple data.table
    
    iteration<-c(round(runif(total.samp,1,dim(mcmc.post)[1]*0.05))) 
    #take only last 20% of the runs to use for posterior distribution for calculations
    
    post.sample[[i.season]]<<-mcmc.post[iteration,] 
    rand.contact.samp[[i.season]]<<-contact.ids[iteration]
    #remove(jj)
    
  }
  post.sample<-list()
  rand.contact.samp<-list()
  total.samp<-num.samp
  #remove(jj)
  #if(i.country==4)
  #{setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', '0.55'))
   # flu.tablesB<-list.files(pattern=glob2rx('GBflu.u5*B'));
  #  flu.tablesH1<-list.files(pattern=glob2rx('GBflu.u5*H1'));
   # flu.tablesH3<-list.files(pattern=glob2rx('GBflu.u5*H3'));
    
    #ct.tablesB<-list.files(pattern=glob2rx('GBct.u5*B'))
    #ct.tablesH1<-list.files(pattern=glob2rx('GBct.u5*H1'))
    #ct.tablesH3<-list.files(pattern=glob2rx('GBct.u5*H3'))
  #}else
    {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
      ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version,'*H1', sep='')))
      ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.',version,'*H3', sep='')))
      ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.', version,'*B', sep='')))
      
      flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.',version,'*B',sep='')));
      flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.',version, '*H1',sep='')));
      flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.', version, '*H3',sep='')));
    }
  
  if(strain==1){year.length<-length(ct.tablesH1)}
  if(strain==2){year.length<-length(ct.tablesH3)}
  if(strain==3){year.length<-length(ct.tablesB)}
  
  #jj<-where('year.length')
  sapply(1:year.length, FUN= function(x) mcmc.samp(x, num.samp))
  #rm(i.country)
  p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving
  
  ###Save output to file
  kname<-paste0(fname[i.country],'PriorSample',version, strain.name[strain])  
  save(post.sample,file=kname) 
  jname<-paste0(fname[i.country],'PriorContactSample',version,strain.name[strain])  
  save(rand.contact.samp, file=jname)
}

n.samples<-500
for(i.country in 1:length(fname))
{
SENS.prior.pars(1,i.country, n.samples, version = 'v1')
}


prior.ct.sampler<-function(program, i.country, strain.choice, version)
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PriorContactSample',version, strain.name[strain.choice])))
  
  if(length(ct.samp)==0)
  {
    ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PriorContactSample', strain.name[strain.choice])))
    #flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample', strain.name[strain.choice])))
    print(paste(fname[i.country], 'base version used'))
  }
  
  load(ct.samp);
  
  polymod<<-NULL
  polymod<-as.data.frame(fread(file.path('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata',paste(allpolymod[i.country]))))
  polymod<-cbind(polymod,polymod$V12) #add column for 75+ year olds
  
  year.length<-length(rand.contact.samp);
  post.means<-list(); post.sd<-list()
  
  for(i.season in 1:year.length)
  {
    if(strain==1) {contact.set<-rand.contact.samp[[i.season]]; }
    if(strain==2) {contact.set<-rand.contact.samp[[i.season]]; }
    
    p.contacts<-NULL
    
    for (jj in 1:length(contact.set))
    {
      p.contacts[[jj]]<-contact_matrix(as.matrix(polymod[contact.set[[jj]],]),
                                       age_sizes[,1], age.group.limits)
    }
    post.means[[i.season]]<-Reduce('+', p.contacts)/length(p.contacts)
    post.sd[[i.season]]<- sqrt(Reduce('+', lapply(Map('-',  p.contacts, post.means[[i.season]]), function(x) x^2))/length(p.contacts))
  }
  
  save(post.means,file=paste0('prior.means',fname[i.country], strain.name[strain]))
  save(post.sd,file=paste0('prior.sd',fname[i.country], strain.name[strain]))
  
}


for(i.country in 1:length(fname))
{prior.ct.sampler(1, i.country, strain.choice = 1, version='u1')}



##################################################################
# Posterior and Prior Contact Matrix in Discrete and Continuous 
##############################################################
country.ct.list.m<-list()
country.ct.list.sd<-list()

version<-'v1'
for(i.country in 1:length(fname))
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
  load(paste0('post.means',fname[i.country], version, 'H1N1')); 
  load(paste0('post.sd',fname[i.country], version, 'H1N1')); 
  
  #load(paste0('prior.means',fname[i.country], 'H1N1')); prior.meansH1N1<-prior.means
 
  
   post.m.list<-list()
   post.sd.list<-list()
  
  for(i.season in 1:14)
  {
    #pull means for each season 
    post.mat.s<-post.means[[i.season]]*age.group.sizes
    post.mat.sd<-post.sd[[i.season]]*age.group.sizes
    
    colnames(post.mat.s)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    rownames(post.mat.s)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    
    colnames(post.mat.sd)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    rownames(post.mat.sd)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    
    post.m.list[[i.season]]<-post.mat.s 
    post.sd.list[[i.season]]<-post.mat.sd 
    
  }
  
  
  country.ct.list.m[[i.country]]<-Reduce('+', post.m.list)/length(post.m.list)
  country.ct.list.sd[[i.country]]<-Reduce('+', post.sd.list)/length(post.sd.list)
  
}

names(country.ct.list.m)<-fname[1:10]
names(country.ct.list.sd)<-fname[1:10]





###############Canberra distance



for(i.country in 1:length(fname))
{
  if(i.country==4) next
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[4])))
  load(paste0('post.means','GB', version, 'H1N1')); gbpost.means<-post.means
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  load(paste0('post.means',fname[i.country], version, 'H1N1')); 

  top<-lapply(Map('-', gbpost.means[1:14], post.means[1:14]), function(x) abs(x))
  bottom<-Map('+', gbpost.means[1:14], post.means[1:14])
  
  season.can.dist<-lapply(lapply(Map('/', top, bottom), na.omit), sum)
  allseason.can.dist<-sum(unlist(lapply(season.can.dist, sum)))
  
  season.bind<-cbind(rep(fname[i.country], 14), c(1:14), unlist(season.can.dist))
  all.bind<-cbind(rep(fname[i.country], 1), unlist(allseason.can.dist))
  
  ifelse(i.country==1, sbind.out<-season.bind, sbind.out<-rbind(season.bind, sbind.out))
  ifelse(i.country==1, allbind.out<-all.bind, allbind.out<-rbind(all.bind, allbind.out))
}
         

  
  

############make paged average plots of contact matrix and comparison of rates via boxplots
library(tidyverse)
library(ggpubr)
library(gridExtra)


for(i.country in 1:(length(fname)-1))
{
dat<-melt(country.ct.list.m[[i.country]])
## reshape data (tidy/tall form)
#levels(dat$X2)<-factor(dat$X2, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
levels(dat$X1)<-factor(dat$X1, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
dat<-dat%>% arrange(X1)
levels(dat$X2)<-factor(dat$X2, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
## plot data
jhk<-ggplot(dat, aes(X1, X2)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red")+labs(title=paste('Average Number \n of Contacts per day', fname[i.country] ), x='Participants', y='Contacts')+theme(text = element_text(size=16), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 

o.name<-paste0('grid.mean', i.country)
assign(o.name, jhk)
}

for(i.country in 1:(length(fname)-1))
{
  dat<-melt(country.ct.list.sd[[i.country]])
  ## reshape data (tidy/tall form)
  levels(dat$X1)<-factor(dat$X1, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  dat<-dat%>% arrange(X1)
  levels(dat$X2)<-factor(dat$X2, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  ## plot data
  jhk<-ggplot(dat, aes(X1, X2)) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(value, 1))) +
    scale_fill_gradient(low = "white", high = "red")+labs(title=paste('Standard Deviation \n Among Number of Contacts per day', fname[i.country] ), x='Participants', y='Contacts')+theme(text = element_text(size=16), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  
  o.name<-paste0('grid.sd', i.country)
  assign(o.name, jhk)
}

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')

#graph CM as matrix with values like heatmap
mg<-lapply(1:10, function(yy) {ggarrange(get(paste0('grid.mean',yy)),get(paste0('grid.sd',yy)),  ncol = 2)})

graphs.allpost<-marrangeGrob(mg,
                             #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                             #labels=c('A'),
                             nrow = 2, ncol=1)
ggsave(paste0("mean.sd.mtx",strain.name[strain], version,".pdf"),device='pdf',plot=graphs.allpost, width=50,height=28,units='cm')


#####################Relative Contacts
country.rel.ct.m<-list()
country.rel.ct.sd<-list()

for(i.country in 1:length(fname))
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  
  load(paste0('post.means',fname[i.country], version,'H1N1')); 
  load(paste0('post.sd',fname[i.country], version, 'H1N1')); 
  
  #load(paste0('prior.means',fname[i.country], 'H1N1')); prior.meansH1N1<-prior.means
  
  rel.post.m.list<-list()
  rel.post.sd.list<-list()
  
  for(i.season in 1:14)
  {
    #pull means for each season 
    rel.post.mat.s<-(post.means[[i.season]])/colSums(post.means[[i.season]])
    rel.post.mat.sd<-(post.sd[[i.season]])/colSums(post.sd[[i.season]])
    
    colnames(rel.post.mat.s)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    rownames(rel.post.mat.s)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    
    colnames(rel.post.mat.sd)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    rownames(rel.post.mat.sd)<-c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
    
    rel.post.m.list[[i.season]]<-rel.post.mat.s 
    rel.post.sd.list[[i.season]]<-rel.post.mat.sd 
    
  }
  
  country.rel.ct.m[[i.country]]<-Reduce('+', rel.post.m.list)/length(rel.post.m.list)
  country.rel.ct.sd[[i.country]]<-Reduce('+', rel.post.sd.list)/length(rel.post.sd.list)
}

names(country.rel.ct.m)<-fname[1:10]
names(country.rel.ct.sd)<-fname[1:10]


for(i.country in 1:(length(fname)-1))
{
  dat<-melt(country.rel.ct.m[[i.country]])
  ## reshape data (tidy/tall form)
  levels(dat$X1)<-factor(dat$X1, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  dat<-dat%>% arrange(X1)
  levels(dat$X2)<-factor(dat$X2, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  ## plot data
  jhk.rel<-ggplot(dat, aes(X1, X2)) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(value, 1))) +
    scale_fill_gradient(low = "white", high = "blue")+labs(title=paste('Proportion of Contacts \n per Age group ', fname[i.country] ), x='Participants', y='Contacts')+theme(text = element_text(size=16), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  
  op.name<-paste0('rel.mean', i.country)
  assign(op.name, jhk.rel)
}

for(i.country in 1:(length(fname)-1))
{
  dat<-melt(country.rel.ct.sd[[i.country]])
  ## reshape data (tidy/tall form)
  levels(dat$X1)<-factor(dat$X1, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  dat<-dat%>% arrange(X1)
  levels(dat$X2)<-factor(dat$X2, levels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))
  ## plot data
  jhk.ff<-ggplot(dat, aes(X1, X2)) +
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(value, 1))) + 
    scale_fill_gradient(low = "white", high = "blue")+labs(title=paste('Standard Deviation Among Proportion \n of Contact per Age Group', fname[i.country] ), x='Participants', y='Contacts')+theme(text = element_text(size=16), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  
  ou.name<-paste0('rel.sd', i.country)
  assign(ou.name, jhk.ff)
}

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')

#graph CM as matrix with values like heatmap
mg<-lapply(1:10, function(yy) {ggarrange(get(paste0('rel.mean',yy)),get(paste0('rel.sd',yy)),  ncol = 2)})

graphs.allpost<-marrangeGrob(mg,
                             #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                             #labels=c('A'),
                             nrow = 2, ncol=1)
ggsave(paste0("realtive.mean.sd.mtx",strain.name[strain],version,".pdf"),device='pdf',plot=graphs.allpost, width=50,height=28,units='cm')




####

#graph CM as matrix with values like heatmap
mgg<-lapply(1:10, function(yy) {ggarrange(get(paste0('grid.mean',yy)),get(paste0('rel.mean',yy)),  ncol = 2)})

graphs.allpost2<-marrangeGrob(mgg,
                             #ggarrange(get(paste0('post',yy)), nrow=1, labels=c('C')),
                             #labels=c('A'),
                             nrow = 2, ncol=1)
ggsave(paste0("realtive.ct.mean.mtx",strain.name[strain],version,".pdf"),device='pdf',plot=graphs.allpost2, width=50,height=28,units='cm')











########boxplots

for(i.country in 1:(length(fname)-1))
{
  dat<-melt(country.ct.list.sd[[i.country]])

  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  load(paste0('prior.means',fname[i.country], 'H1N1')); prior.means<-post.means
  
  load(paste0('post.means',fname[i.country], version, 'H1N1')); 

  #post.mat.s<-post.means[[i.season]]*age.group.sizes
  #post.mat.sd<-post.sd[[i.season]]*age.group.sizes

#load('BEpost.meansu3H1N1')
#load('BEprior.meansu3H1N1')
pt.means<-Map('*',post.means,age.group.sizes)
pr.means<-Map('*',prior.means,age.group.sizes)
post.mat.mean<-apply(array(unlist(pt.means), c(11,11,2000)), c(1,2), mean)
prior.mat.mean<-apply(array(unlist(pr.means), c(11,11,300)), c(1,2), mean)
post.m.ggplot <- as.data.frame(melt(pt.means)) #for 11 groups
pre.m.ggplot<- as.data.frame(melt(pr.means)) #for 11 groups


if(dim(post.mat.mean)[1]==7)
{
#Code for 7 compartment plot
col.post<-cbind(post.mat.mean[,1],(post.mat.mean[,2]+post.mat.mean[,3]),(post.mat.mean[,4]+post.mat.mean[,5]),(post.mat.mean[,6]+post.mat.mean[,7]), post.mat.mean[,8],post.mat.mean[,9], (post.mat.mean[,10]+post.mat.mean[,11]))
post.mat.s<-cbind(col.post[1,],(col.post[2,]+col.post[3,]),(col.post[4,]+col.post[5,]),(col.post[6,]+col.post[7,]), col.post[8,],col.post[9,], (col.post[10,]+col.post[11,]))

post.m.ggplot <- melt(post.mat.s)
#prior
#Code for 7 compartment plot
col.prior<-cbind(prior.mat.mean[,1],(prior.mat.mean[,2]+prior.mat.mean[,3]),(prior.mat.mean[,4]+prior.mat.mean[,5]),(prior.mat.mean[,6]+prior.mat.mean[,7]), prior.mat.mean[,8],prior.mat.mean[,9], (prior.mat.mean[,10]+prior.mat.mean[,11]))
prior.mat.s<-cbind(col.prior[1,],(col.prior[2,]+col.prior[3,]),(col.prior[4,]+col.prior[5,]),(col.prior[6,]+col.prior[7,]), col.prior[8,],col.prior[9,], (col.prior[10,]+col.prior[11,]))
pre.m.ggplot<- melt(prior.mat.s)


#7 compartment plots
comp1<-ggplot(data=post.m.ggplot, aes(x=X1,y=X2,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Posterior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day",strain.name[strain]))+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

#comp2<-ggplot(data=pre.m.ggplot, aes(x=X1,y=X2,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Prior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day",strain.name[strain]))+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

#section data into discrete intervals
discrete1<-cut(post.m.ggplot$value, breaks=c(0,0.1e-5,1,5,7,10, 12, 15,18,20), right=FALSE, ordered_result = FALSE)
post.discrete<-cbind(post.m.ggplot, discrete1)
#discrete2<-cut(pre.m.ggplot$value, breaks=c(0,0.1e-5,0.2,0.4,0.6,1,1.25,3,5,7), right=FALSE, ordered_result = FALSE)
#prior.discrete<-cbind(pre.m.ggplot, discrete2)
#dev.off()
#dev.off()

p3<-ggplot(data=post.discrete, aes(x=X1,y=X2,fill=factor(discrete1)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Posterior Contact Matrix\n (Discrete Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-5,1,5,7,10, 12, 15,18,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-1)','[1-5)','[5-7)','[7-10)', '[10-12)','[12-15)','[15-18)','[18,20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

#p4<-ggplot(data=prior.discrete, aes(x=X1,y=X2,fill=factor(discrete2)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Prior Contact Matrix\n (Discrete Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-5,0.2,0.4,0.6,1,1.25,3,5,7), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-0.2)','[0.2,0.4)','[0.4,0.6)','[0.6-1)', '[1-1.25)','[1.25,3)','[3,5)','[5,7)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))


#plot_grid(comp2, comp1, p4,p3, nrow=2, labels = "AUTO", align = 'h', label_size = 14)
#ggsave(paste0("Image3",fname[i.country],strain.name[strain],".tiff"),device='tiff',plot=last_plot(), width=45,height=45,units='cm')
#dev.off()

xx<-plot_grid(comp1, p3, nrow=1, labels = "AUTO", align = 'h', label_size = 14)
ggsave(filename=paste0("Image.exmat",fname[i.country],strain.name[strain],'7','.png'), plot=xx, width=45,height=35,units='cm')
}

################
comp2<-comp1<-p4<-p3<-NULL

#comp1<-p3<-NULL

#11 compartment plots------------------
comp1<-ggplot(data=post.m.ggplot, aes(x=X2,y=X1,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Posterior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day"))+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

comp2<-ggplot(data=pre.m.ggplot, aes(x=X2,y=X1,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Prior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day"))+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

#section data into discrete intervals
discrete1<-cut(post.m.ggplot$value, breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.5,3,5,7,10,20), right=FALSE, ordered_result = FALSE)
post.discrete<-as.data.frame(cbind(post.m.ggplot, discrete1))
discrete2<-cut(pre.m.ggplot$value, breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.5,3,5,7,10,20), right=FALSE, ordered_result = FALSE)
prior.discrete<-as.data.frame(cbind(pre.m.ggplot, discrete2))
#dev.off()
graphics.off()

p3<-ggplot(data=post.discrete, aes(x=X2,y=X1,fill=factor(discrete1)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Posterior Contact Matrix\n (Discrete Scale)',strain.name[strain]))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.25,3,5,7,10,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-0.2)','[0.2,0.4)','[0.4,0.6)','[0.6-1)', '[1-1.25)','[1.25-3)','[3-5)','[5-7)','[7-10)','[10-20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

p4<-ggplot(data=prior.discrete, aes(x=X2,y=X1,fill=factor(discrete2)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Prior Contact Matrix\n (Discrete Scale)',strain.name[strain]))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.25,3,5,7,10,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-0.2)','[0.2,0.4)','[0.4,0.6)','[0.6-1)', '[1-1.25)','[1.25-3)','[3-5)','[5-7)','[7-10)','[10-20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))

all.tog<-plot_grid(comp2, comp1, p4,p3, nrow=2, labels = "AUTO", align = 'h', label_size = 13)

ggsave(filename=paste0("Image.exmat",fname[i.country],strain.name[strain],'11','.png'), plot=all.tog, width=45,height=35,units='cm')
#dev.off()


#all.tog<-plot_grid(comp1, p3, nrow=1, labels = "AUTO", align = 'h', label_size = 14)
#ggsave(filename=paste0("Image.exmat",fname[i.country],strain.name[strain],'11','.tiff'), plot=all.tog, width=35,height=17,units='cm')
}
#dev.off()
#end country loop
 #end strain loop


###################################

for(i.country in 1:(length(fname)-1))
{
  dat<-melt(country.ct.list.sd[[i.country]])
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste0(fname[i.country])))
  load(paste0('prior.means',fname[i.country], 'H1N1')); prior.means<-post.means
  
  load(paste0('post.means',fname[i.country], version, 'H1N1')); 
  
  
  pt.means<-Map('*',post.means,age.group.sizes)
  pr.means<-Map('*',prior.means,age.group.sizes)
  post.mat.mean<-apply(array(unlist(pt.means), c(11,11,2000)), c(1,2), mean)
  prior.mat.mean<-apply(array(unlist(pr.means), c(11,11,300)), c(1,2), mean)
  post.m.ggplot <- as.data.frame(melt(pt.means)) #for 11 groups
  pre.m.ggplot<- as.data.frame(melt(pr.means)) #for 11 groups
  
  
  if(dim(post.mat.mean)[1]==7)
  {
    #Code for 7 compartment plot
    col.post<-cbind(post.mat.mean[,1],(post.mat.mean[,2]+post.mat.mean[,3]),(post.mat.mean[,4]+post.mat.mean[,5]),(post.mat.mean[,6]+post.mat.mean[,7]), post.mat.mean[,8],post.mat.mean[,9], (post.mat.mean[,10]+post.mat.mean[,11]))
    post.mat.s<-cbind(col.post[1,],(col.post[2,]+col.post[3,]),(col.post[4,]+col.post[5,]),(col.post[6,]+col.post[7,]), col.post[8,],col.post[9,], (col.post[10,]+col.post[11,]))
    
    post.m.ggplot <- melt(post.mat.s)
    #prior
    #Code for 7 compartment plot
    col.prior<-cbind(prior.mat.mean[,1],(prior.mat.mean[,2]+prior.mat.mean[,3]),(prior.mat.mean[,4]+prior.mat.mean[,5]),(prior.mat.mean[,6]+prior.mat.mean[,7]), prior.mat.mean[,8],prior.mat.mean[,9], (prior.mat.mean[,10]+prior.mat.mean[,11]))
    prior.mat.s<-cbind(col.prior[1,],(col.prior[2,]+col.prior[3,]),(col.prior[4,]+col.prior[5,]),(col.prior[6,]+col.prior[7,]), col.prior[8,],col.prior[9,], (col.prior[10,]+col.prior[11,]))
    pre.m.ggplot<- melt(prior.mat.s)
    
    
    #7 compartment plots
    comp1<-ggplot(data=post.m.ggplot, aes(x=X1,y=X2,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Posterior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day",strain.name[strain]))+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))
    
    #section data into discrete intervals
    discrete1<-cut(post.m.ggplot$value, breaks=c(0,0.1e-5,1,5,7,10, 12, 15,18,20), right=FALSE, ordered_result = FALSE)
    post.discrete<-cbind(post.m.ggplot, discrete1)
    #discrete2<-cut(pre.m.ggplot$value, breaks=c(0,0.1e-5,0.2,0.4,0.6,1,1.25,3,5,7), right=FALSE, ordered_result = FALSE)
    #prior.discrete<-cbind(pre.m.ggplot, discrete2)
    #dev.off()
    #dev.off()
    
    p3<-ggplot(data=post.discrete, aes(x=X1,y=X2,fill=factor(discrete1)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste0('Posterior Contact Matrix\n (Discrete Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-5,1,5,7,10, 12, 15,18,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-1)','[1-5)','[5-7)','[7-10)', '[10-12)','[12-15)','[15-18)','[18,20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+scale_y_discrete(limits=c(1:7),labels=c('6m-1','1-4','5-14', '15-24','25-44','45-64','65+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))
    
    xx<-plot_grid(comp1, p3, nrow=1, labels = "AUTO", align = 'h', label_size = 14)
    ggsave(filename=paste0("Image.exmat",fname[i.country],strain.name[strain],'7','.png'), plot=xx, width=45,height=35,units='cm')
  }
  
  ################
  comp2<-comp1<-p4<-p3<-NULL
  
  #comp1<-p3<-NULL
  
  #11 compartment plots------------------
  #comp1<-ggplot(data=post.m.ggplot, aes(x=X2,y=X1,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Posterior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day"))+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))
  
  #comp2<-ggplot(data=pre.m.ggplot, aes(x=X2,y=X1,fill=value))+labs(x='Age Group of Participant',y='Age Group of Contact',title=paste('Prior Contact Matrix\n (Continuous Scale)'))+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name=paste("Contacts per\n Person per Day"))+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=13),axis.title=element_text(size=13))
  
  #section data into discrete intervals
  discrete1<-cut(post.m.ggplot$value, breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.5,3,5,7,10,20), right=FALSE, ordered_result = FALSE)
  post.discrete<-as.data.frame(cbind(post.m.ggplot, discrete1))
  discrete2<-cut(pre.m.ggplot$value, breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.5,3,5,7,10,20), right=FALSE, ordered_result = FALSE)
  prior.discrete<-as.data.frame(cbind(pre.m.ggplot, discrete2))
  #dev.off()
  graphics.off()
  
  p3<-ggplot(data=post.discrete, aes(x=X2,y=X1,fill=factor(discrete1)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste('Posterior Contact Matrix\n (Discrete Scale)',strain.name[strain]))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.25,3,5,7,10,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-0.2)','[0.2,0.4)','[0.4,0.6)','[0.6-1)', '[1-1.25)','[1.25-3)','[3-5)','[5-7)','[7-10)','[10-20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=14),axis.title=element_text(size=14))
  
  p4<-ggplot(data=prior.discrete, aes(x=X2,y=X1,fill=factor(discrete2)))+labs(x='Age Groups of Participant',y='Age Group of Contact',title=paste('Prior Contact Matrix\n (Discrete Scale)',strain.name[strain]))+ geom_tile(color="white", size=0.1)+coord_equal()+ scale_color_viridis(breaks=c(0,0.1e-6,0.2,0.4,0.6,1,1.25,3,5,7,10,20), discrete=TRUE)+scale_fill_viridis(option='inferno',labels=c('0','(0-0.2)','[0.2,0.4)','[0.4,0.6)','[0.6-1)', '[1-1.25)','[1.25-3)','[3-5)','[5-7)','[7-10)','[10-20)'), name=paste("Contacts per\n Person per Day"), discrete=TRUE)+scale_x_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+scale_y_discrete(limits=c(1:11),labels=c('6m-1','1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+'))+theme(plot.title=element_text(hjust=0),plot.margin = unit(c(0, 0, 0, 1), "cm"),axis.text=element_text(size=14),axis.title=element_text(size=14))
  
  all.tog<-plot_grid(p4,p3, nrow=1, labels = "AUTO", align = 'h', label_size = 13)
  
  ggsave(filename=paste0("Image.exmat",fname[i.country],strain.name[strain],'11','.png'), plot=all.tog, width=45,height=25,units='cm')
  #dev.off()
}




##################################################################
# 3d plot?
##############################################################

png(file=heatpathall, width=3000,height=5000, bg='white')
#par(mfrow=c(1,2))
#multiplot(p1, p2, cols=2)
library(grid)
library(gridExtra)
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

underhill<-grid_arrange_shared_legend(comp2, p3)

###Plotly 3d plot
devtools::install_github("ropensci/plotly")
library(plotly)
Sys.setenv('plotly_username'='taj.wenzel')
Sys.setenv("plotly_api_key"='MVsC3gXUhwCakszsk14D')
set_credentials_file(username = 'your_username', key = 'your_key')
#get SD
post.mat.sd<-apply(array(unlist(pt.means), c(11,11,out)), c(1,2), sd)
prior.mat.sd<-apply(array(unlist(pr.means), c(11,11,out)), c(1,2), sd)
z2<-prior.mat.mean+1.96*prior.mat.sd
z3<-prior.mat.mean-1.96*prior.mat.sd
labels <- c('75+','64-75','45-64','25-44','17-24','15-16','12-14','5-11','2-4', '6m-1')
p <- plot_ly(showscale = F) %>%
  add_surface(z = ~prior.mat.mean) %>%
  add_surface(z = ~z2, opacity = 0.90) %>%
  add_surface(z = ~z3, opacity = 0.90) %>%
  layout(
    scene = list(
    xaxis = list(title = 'Age Group Participant', tickvals=seq(1,11,1), ticktext=rev(labels), tickangle=30),
    yaxis = list(title = "Age Group Contact", tickvals=seq(1,11,1), ticktext=rev(labels)),
    zaxis = list(title = "Number of Contacts")))

all.tog<-plot_grid(comp2, p, nrow=1, labels = "AUTO", align = 'h', label_size = 13)

chart_link = api_create(p, filename="surface-3")
chart_link

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}
#library(plotly)
post.means<-array()
prior.means<-array()


##########3d Barplot of contact matrix---------------
# draw each bar: from left to right ...

library(graphics)
data<-(abs(GBPL[[1]]))
data[upper.tri(data,diag = FALSE)]<-NA
data<-replace(data,data==Inf,NA)

mar.default <- c(1,4,4,2) + 0.2
par(mar = mar.default + c(2, 0, 0, 0)) 
# generate 'empty' persp plot
pmat = persp(x=c(0,10), y=c(0,10), z=matrix(c(0,0,0,0), nrow=2), 
             xlim=c(0,10), ylim=c(0,10), zlim=c(0,3), 
             xlab='Age Groups', ylab='Age Groups', zlab='Proportion increase', cex.lab=2,
             theta=60, phi=40, d=1, box=T, main=paste0(fname[i.country],'/UK Ratio')) 

# define color ramp
my_cols = heat.colors(20)

# generate color matrix (values between 1 and 10, corresponding to 10 values my_cols
#for GB
colmat = matrix(data = 1, ncol = 10, nrow = 10)
data<-replace(data,data<0,NA)
colmat<-replace(colmat,data==0,1)
colmat<-replace(colmat, data <0, 2)
colmat<-replace(colmat, 0 < data & data <=0.1, 3)
colmat<-replace(colmat, 0.2 < data & data <=0.25, 4)
colmat<-replace(colmat, 0.25 < data & data <=0.3, 5)
colmat<-replace(colmat, 0.3 < data & data <=0.4, 6)
colmat<-replace(colmat, 0.4 < data & data <=0.5, 7)
colmat<-replace(colmat, 0.5 < data & data <=0.6, 8)
colmat<-replace(colmat, 0.07 < data & data <=0.08, 9)
colmat<-replace(colmat, 0.08 < data & data <=0.09, 10)
colmat<-replace(colmat, 0.09 < data & data <=0.1, 11)
colmat<-replace(colmat, 0.1 < data & data <=0.11, 12)
colmat<-replace(colmat, 0.11 < data & data <=0.12, 13)

#for difference in prior/posterior

i<-NULL
j<-NULL
for (i in 1:ncol(data))
{
  # ... and back to front 
  #j<-9
  for (j in nrow(data):1)
  {
    xy <- which(data == data[i,j], arr.ind=TRUE)
    
    # side facing y
    x = rep(xy[1],4)
    y = c(xy[2]-1,xy[2],xy[2],xy[2]-1)
    z = c(0,0,data[i,j],data[i,j])
    polygon(trans3d(x, y, z, pmat), col=my_cols[colmat[i,j]], border=1)
    
    #  side facing x
    x = c(xy[1]-1,xy[1],xy[1],xy[1]-1)
    y = rep(xy[2]-1,4)
    z = c(0,0,data[i,j],data[i,j])
    polygon(trans3d(x, y, z, pmat), col=my_cols[colmat[i,j]], border=1)
    
    # top side
    x = c(xy[1]-1,xy[1],xy[1],xy[1]-1)
    y = c(xy[2]-1,xy[2]-1,xy[2],xy[2])
    z = rep(data[i,j],4)
    polygon(trans3d(x, y, z, pmat), col=my_cols[colmat[i,j]], border=1)
  }
}

# define axis ranges etc
x.axis <- 1:ncol(data) - 0.5
min.x <- 0
max.x <- 11
y.axis <- 1:nrow(data) - 0.5 
min.y <- 0
max.y <- 11
z.axis <- seq(0, 1, by=0.1)
min.z <- 0
max.z <- 1

# add some distance between tick labels and the axis
xoffset = 2
yoffset = 1
zoffset = 1
ticklength = 0.2

# x axis ticks
tick.start <- trans3d(x.axis, min.y, min.z, pmat)
tick.end <- trans3d(x.axis, (min.y - ticklength), min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

# y axis ticks
tick.start <- trans3d(max.x, y.axis, min.z, pmat)
tick.end <- trans3d(max.x + ticklength, y.axis, min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

# z axis ticks
tick.start <- trans3d(min.x, min.y, z.axis, pmat)
tick.end <- trans3d(min.x, (min.y - ticklength), z.axis, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

# x labels
labels <- c('75+','64-75','45-64','25-44','17-24','15-16','12-14','5-11','2-4', '6m-1')
label.pos <- trans3d(x.axis-0.5, (min.y - xoffset), min.z, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=1, cex=2)

# y labels
labels <- c('6m-1','2-4','5-11', '12-14', '15-16', '17-24','25-44','45-64','65-74', '75+')
label.pos <- trans3d((max.x + yoffset), y.axis+0.5, min.z, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=1, cex=2)

# z labels
labels <- as.character(z.axis)
label.pos <- trans3d(min.x, (min.y - zoffset), z.axis, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), srt=0, cex=2) 
