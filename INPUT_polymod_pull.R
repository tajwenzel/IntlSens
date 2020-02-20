#PolyMod DataManip
rm(list=ls())
setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata")
library(dplyr)
library(fluEvidenceSynthesis)
library(compare)
library(tidyr)
library(plyr)
#data("age_sizes")
#touchtype=1
#data("polymod_uk")
rm(list = ls())


polymodpull<-function(agesize, touchtype)
  {
  rawpolymod<-read.table("contacts_final.txt",header=T)
  rawsubject<-read.table("participants_final.txt",header=T)
  
  merged<-merge(rawpolymod,rawsubject, by=c("global_id"))
  #merged<-merge(rawpolymod,rawsubject, by=c("part_id"))
  clean<-!is.na(as.numeric(as.character(merged$participant_age)))&!is.na(as.numeric(as.character(merged$dayofweek)))
  merged<-merged[clean,]
  
  d.merged<-select(merged,global_id,local_id.x,country.x,participant_age,dayofweek, cnt_count.x,cnt_age_mean,cnt_home,cnt_leisure,cnt_work,cnt_school,cnt_transport,cnt_otherplace,cnt_touch, holiday) 
  
  #country levels BE DE FI GB IT LU NL PL
  BE<-d.merged[d.merged$country.x=="BE",]
  DE<-d.merged[d.merged$country.x=='DE',]
  FI<-d.merged[d.merged$country.x=='FI',]
  GB<-d.merged[d.merged$country.x=='GB',]
  IT<-d.merged[d.merged$country.x=='IT',]
  LU<-d.merged[d.merged$country.x=='LU',]
  NL<-d.merged[d.merged$country.x=='NL',]
  PL<-d.merged[d.merged$country.x=='PL',]
  
  combo<-list(BE, DE, FI, GB, IT, LU, NL, PL)
  #######################GENERATE WORK DATASET######################################################
  #lapply(combo,attributes)
  #kk=4
 
  
  for (kk in 1:length(combo))
  {
    merge<-combo[[kk]];
    workdata<-!is.na(merge$dayofweek)&!is.na(merge$cnt_age_mean)&!is.na(merge$participant_age)&!is.na(merge$holiday)&merge$holiday==0&merge$cnt_touch==as.numeric(touchtype);
    workdata<-as.data.frame(merge[workdata,])
    
    #lapply(set, class)
    #check
    #which(is.na(as.numeric(as.character(workdata$participant_age))))
    
    
    ###For AGE upper/lower bound
    #for (i in 1:length(workdata$cnt_age_l))
    #    { 
    #ifelse (is.na(workdata$cnt_age_l[i])==TRUE, workdata$age_group_l[i]<-NA, workdata$age_group_l[i]<-as.age.group(workdata$cnt_age_l[i], limits = as.numeric(c(1, 5, 15, 25, 45, 65))));
    #ifelse (is.na(workdata$cnt_age_r[i])==TRUE, workdata$age_group_r[i]<-NA,workdata$age_group_r[i]<-as.age.group(workdata$cnt_age_r[i], limits = as.numeric(c(1, 5, 15, 25, 45, 65))));
    #}
    
    ###For AGE mean
    for (i in 1:length(workdata$cnt_age_mean))
    { options(na.action='na.pass')
      workdata$age_group_mean[i]<-as.age.group(workdata$cnt_age_mean[i], limits = as.numeric(c(1,2,5,12,15,16,25,45,65,86))) #upper limits);
     }
    #check for age group mismatch
    #for (i in 1:length(workdata$cnt_age_l))
    #{ 
    #  ifelse (is.na(workdata$cnt_age_l[i])==TRUE, workdata$age_group_l[i]<-NA, workdata$age_group_l[i]<-as.age.group(workdata$cnt_age_l[i], limits = as.numeric(c(1, 5, 15, 25, 45, 65))));
    #  ifelse (is.na(workdata$cnt_age_r[i])==TRUE, workdata$age_group_r[i]<-NA,workdata$age_group_r[i]<-as.age.group(workdata$cnt_age_r[i], limits = as.numeric(c(1, 5, 15, 25, 45, 65))));
    #}
    #for (i in 1:length(workdata$age_group_l)) {
    #  if (is.na(workdata$age_group_r[i]) || workdata$age_group_r[i] == workdata$age_group_l[i]) {
    #   workdata$age_group[i] <- workdata$age_group_l[i]
    #} else {
    # workdata$age_group[i] <- NA
    #}
    #}
    #which(is.na(workdata$holiday))
    #summary(workdata$holiday)
    #workdata<-subset(workdata,workdata$holiday==0)
    #head(workdata)
    
    #which(is.na(workdata$dayofweek))
    #summary(workdata$dayofweek)
    #workdata$day<-workdata$dayofweek
    #workdata$day[workdata$day==6] <- 0
    
    #for (j in 1:5) {
     # workdata$day[workdata$day==j ] <-1
    #}
    #table(workpart$dayofweek==0|workpart$dayofweek==6)
    
    options(na.action='na.pass')
    agegrp<-(as.factor(as.matrix(workdata$age_group_mean)))
    B <- as.data.frame(model.matrix(workdata$local_id.x ~ 0+agegrp))
    ###special conditions for italy
    #if(kk==5) 
     # {agegrp1<-as.data.frame(rep(0,dim(B)[1]));
      #  B<-cbind(agegrp1,B);
       #   colnames(B)=c('agegrp1','agegrp2','agegrp3','agegrp4','agegrp5','agegrp6','agegrp7','agegrp8','agegrp9')}
    #####
    day<-as.numeric(workdata$dayofweek==0|workdata$dayofweek==6)
    workdata<-na.omit(cbind(workdata,B,day))
    
    GBfinal<-select(workdata,local_id.x,participant_age,agegrp1,agegrp2,agegrp3,agegrp4,agegrp5,agegrp6,agegrp7,agegrp8,agegrp9,agegrp10,day)
    
    #options(na.action='na.omit')
    set<-ddply(GBfinal,c("local_id.x","participant_age","day"),numcolwise(sum))
    set<-set[complete.cases(set),];
    set$participant_age<-na.pass(as.numeric(as.character(set$participant_age)))
    set$day<-na.omit(as.numeric(as.character(set$day)))
    set<-na.omit(set)
    colnames(set)=c('V0','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
    
    newpolymod<-set[,-1]
    colnames(newpolymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
    final.polymod<-arrange(newpolymod,V1)
    
    sname<-c('BEtable11', 'DEtable11', 'FItable11', 'GBtable11', 'ITtable11', 'LUtable11', 'NLtable11', 'PLtable11')
    
    write.csv(final.polymod,paste0(sname[kk],".csv"),row.names = FALSE)
  }
}

polymodpull(age_sizes,touchtype=1)
#set mismatch check
hist(polymod$V1, col=2, breaks=seq(0,100,5))
hist(polymod_uk$V1, col=3, breaks=seq(0,100,5), xlab='Age Distribution of Contacts', main='Differences in polymod data', add=T)
#trial<-read.csv('PLtable.csv',header=TRUE,sep=',')

#sum(is.na(workdata$age_group))
#unmatched<-mutate(workdata,new=workdata$age_group_r+workdata$age_group_l);
#is.odd <- function(x) x %% 2 != 0
#yes<-subset((test$new),is.odd(test$new)==TRUE)
#hist(yes,freq=TRUE, xlim=range(0:14))
#hist(test$new, freq=TRUE,xlim=range(0:14))

#check<-setdiff(polymod$V1,polymod_uk$V1)
#compare(polymod$V1,polymod_uk$V1)


####--------------------------------------------------------------------------------
setwd("/Users/Natasha/Downloads")

if(FR=T)
{
rawpolymod<-read.csv('2015_Beraud_France_contact_common.csv')
rawsubject<-read.csv('2015_Beraud_France_participant_common.csv')
rawdays<-read.csv('2015_Beraud_France_sday.csv')
}

merged1<-merge(rawpolymod,rawsubject, by=c("part_id"))
merged2<-merge(merged1,rawdays, by=c("part_id") )
clean<-!is.na(as.numeric(as.character(merged2$part_age)))&!is.na(as.numeric(as.character(merged2$dayofweek)))
merged<-merged2[clean,]

d.merged<-select(merged,part_id,part_age,dayofweek, cnt_age_exact, cnt_age_est_min, cnt_age_est_max, cnt_home,cnt_leisure,cnt_work,cnt_school,cnt_transport,cnt_otherplace,phys_contact, holiday) 


cnt_age_mean<-round((d.merged$cnt_age_est_min+d.merged$cnt_age_est_max)/2)
#cnt_age_mean[is.na(cnt_age_mean)]<-0
d.merged$cnt_age_exact[is.na(d.merged$cnt_age_exact)]<-cnt_age_mean[is.na(d.merged$cnt_age_exact)]
cnt_age_all<-d.merged$cnt_age_exact


e.merged<-cbind(d.merged$part_id,d.merged$part_age,d.merged$dayofweek, cnt_age_all, d.merged$cnt_home,d.merged$cnt_leisure,d.merged$cnt_work,d.merged$cnt_school,d.merged$cnt_transport,d.merged$cnt_otherplace,d.merged$phys_contact, d.merged$holiday)

colnames(e.merged)<-c(colnames(d.merged)[1:3],'cnt_age_all',colnames(d.merged)[7:dim(d.merged)[2]])

f.merged<-e.merged[,colSums(is.na(e.merged))!=nrow(e.merged) ]

touchtype<-1
merge<-as.data.table(f.merged)

merge<-merge[!is.na(merge$dayofweek),]
merge<-merge[!is.na(merge$cnt_age_all),]
merge<-merge[!is.na(merge$part_age),]
merge<-merge[!is.na(merge$holiday),]
merge<-merge[merge$holiday==1,]

workdata<-!is.na(merge$dayofweek)&!is.na(merge$cnt_age_all)&!is.na(merge$part_age)&!is.na(merge$holiday)&merge$holiday==1&merge$phys_contact==as.numeric(touchtype);
workdata<-as.data.frame(merge[workdata,])

for (i in 1:length(workdata$cnt_age_all))
{ options(na.action='na.pass')
  workdata$age_group_mean[i]<-as_age_group(workdata$cnt_age_all[i], limits = as.numeric(c(1,2,5,12,15,16,25,45,65,86))) #upper limits);
}

options(na.action='na.pass')
agegrp<-(as.factor(as.matrix(workdata$age_group_mean)))
B <- as.data.frame(model.matrix(workdata$part_id ~ 0+agegrp))

day<-as.numeric(workdata$dayofweek==0|workdata$dayofweek==6)
workdata<-na.omit(cbind(workdata,B,day))

GBfinal<-select(workdata,part_id,part_age,agegrp1,agegrp2,agegrp3,agegrp4,agegrp5,agegrp6,agegrp7,agegrp8,agegrp9,agegrp10,day)

#options(na.action='na.omit')
set<-ddply(GBfinal,c("part_id","part_age","day"),numcolwise(sum))
set<-set[complete.cases(set),];
set$part_age<-na.pass(as.numeric(as.character(set$part_age)))
set$day<-na.omit(as.numeric(as.character(set$day)))
set<-na.omit(set)
colnames(set)=c('V0','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')

newpolymod<-set[,-1]
colnames(newpolymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
final.polymod<-arrange(newpolymod,V1)

sname<-c('BEtable11', 'DEtable11', 'FItable11', 'GBtable11', 'ITtable11', 'LUtable11', 'NLtable11', 'PLtable11','FRtable11')

write.csv(final.polymod,paste0(sname[9],".csv"),row.names = FALSE)

polymod<-cbind(final.polymod,final.polymod$V12)

contact_matrix(as.matrix(polymod), age_sizes[,1], age_group_limits = age.group.limits)*age.group.sizes

###------------
rm(list = ls())
setwd("/Users/Natasha/Downloads")

if(ZI=T)
{
  rawpolymod<-read.csv('2017_Melegaro_Zimbabwe_contact_common.csv')
  rawsubject<-read.csv('2017_Melegaro_Zimbabwe_participant_common.csv')
  rawdays<-read.csv('2017_Melegaro_Zimbabwe_sday.csv')
}

merged1<-merge(rawpolymod,rawsubject, by=c("part_id"))
merged2<-merge(merged1,rawdays, by=c("part_id") )
clean<-!is.na(as.numeric(as.character(merged2$part_age)))&!is.na(as.numeric(as.character(merged2$dayofweek)))
merged<-merged2[clean,]

d.merged<-select(merged,part_id,part_age,dayofweek, cnt_age_exact, cnt_age_est_min, cnt_age_est_max, cnt_home,cnt_leisure,cnt_work,cnt_school,cnt_transport,cnt_otherplace,phys_contact) 


cnt_age_mean<-round((d.merged$cnt_age_est_min+d.merged$cnt_age_est_max)/2)
#cnt_age_mean[is.na(cnt_age_mean)]<-0
d.merged$cnt_age_exact[is.na(d.merged$cnt_age_exact)]<-cnt_age_mean[is.na(d.merged$cnt_age_exact)]
cnt_age_all<-d.merged$cnt_age_exact


e.merged<-cbind(d.merged$part_id,d.merged$part_age,d.merged$dayofweek, cnt_age_all, d.merged$cnt_home,d.merged$cnt_leisure,d.merged$cnt_work,d.merged$cnt_school,d.merged$cnt_transport,d.merged$cnt_otherplace,d.merged$phys_contact, d.merged$holiday)

colnames(e.merged)<-c(colnames(d.merged)[1:3],'cnt_age_all',colnames(d.merged)[7:dim(d.merged)[2]])

f.merged<-e.merged[,colSums(is.na(e.merged))!=nrow(e.merged) ]

touchtype<-1
merge<-as.data.table(f.merged)

merge<-merge[!is.na(merge$dayofweek),]
merge<-merge[!is.na(merge$cnt_age_all),]
merge<-merge[!is.na(merge$part_age),]
merge<-merge[!is.na(merge$holiday),]
merge<-merge[merge$holiday==1,]

workdata<-!is.na(merge$dayofweek)&!is.na(merge$cnt_age_all)&!is.na(merge$part_age)&merge$phys_contact==as.numeric(touchtype);
workdata<-as.data.frame(merge[workdata,])

for (i in 1:length(workdata$cnt_age_all))
{ options(na.action='na.pass')
  workdata$age_group_mean[i]<-as_age_group(workdata$cnt_age_all[i], limits = as.numeric(c(1,2,5,12,15,16,25,45,65,86))) #upper limits);
}

options(na.action='na.pass')
agegrp<-(as.factor(as.matrix(workdata$age_group_mean)))
B <- as.data.frame(model.matrix(workdata$part_id ~ 0+agegrp))

day<-as.numeric(workdata$dayofweek==0|workdata$dayofweek==6)
workdata<-na.omit(cbind(workdata,B,day))
workdata$agegrp1<-rep(0,nrow(workdata))

GBfinal<-select(workdata,part_id,part_age,agegrp1,agegrp2,agegrp3,agegrp4,agegrp5,agegrp6,agegrp7,agegrp8,agegrp9,agegrp10,day)

#options(na.action='na.omit')
set<-ddply(GBfinal,c("part_id","part_age","day"),numcolwise(sum))
set<-set[complete.cases(set),];
set$part_age<-na.pass(as.numeric(as.character(set$part_age)))
set$day<-na.omit(as.numeric(as.character(set$day)))
set<-na.omit(set)
colnames(set)=c('V0','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')

newpolymod<-set[,-1]
colnames(newpolymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
final.polymod<-arrange(newpolymod,V1)

sname<-c('BEtable11', 'DEtable11', 'FItable11', 'GBtable11', 'ITtable11', 'LUtable11', 'NLtable11', 'PLtable11','FRtable11', 'ZItable11')

write.csv(final.polymod,paste0(sname[10],".csv"),row.names = FALSE)

polymod<-cbind(final.polymod,final.polymod$V12)

contact_matrix(as.matrix(polymod), age_sizes[,1], age_group_limits = age.group.limits)
contact_matrix(as.matrix(polymod), age_sizes[,1], age_group_limits = age.group.limits)*age.group.sizes
age.group.limits<<-c(1,2,5,12,15,16,25,45,65,75) #upper limits
age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)


####-------------------------Peru

rm(list = ls())
setwd("/Users/Natasha/Downloads")

if(PE=T)
{
  rawpolymod<-read.csv('2015_Grijalva_Peru_contact_common.csv')
  rawsubject<-read.csv('2015_Grijalva_Peru_participant_common.csv')
  rawdays<-read.csv('2015_Grijalva_Peru_sday.csv')
}

merged1<-merge(rawpolymod,rawsubject, by=c("part_id"))
merged2<-merge(merged1,rawdays, by=c("part_id") )
clean<-!is.na(as.numeric(as.character(merged2$part_age)))&!is.na(as.numeric(as.character(merged2$dayofweek)))
merged<-merged2[clean,]

d.merged<-select(merged,part_id,part_age,dayofweek, cnt_age_exact, cnt_age_est_min, cnt_age_est_max, cnt_home,cnt_leisure,cnt_work,cnt_school,cnt_transport,cnt_otherplace,phys_contact) 


cnt_age_mean<-round((d.merged$cnt_age_est_min+d.merged$cnt_age_est_max)/2)
#cnt_age_mean[is.na(cnt_age_mean)]<-0
d.merged$cnt_age_exact[is.na(d.merged$cnt_age_exact)]<-cnt_age_mean[is.na(d.merged$cnt_age_exact)]
cnt_age_all<-d.merged$cnt_age_exact


e.merged<-cbind(d.merged$part_id,d.merged$part_age,d.merged$dayofweek, cnt_age_all, d.merged$cnt_home,d.merged$cnt_leisure,d.merged$cnt_work,d.merged$cnt_school,d.merged$cnt_transport,d.merged$cnt_otherplace,d.merged$phys_contact, d.merged$holiday)

colnames(e.merged)<-c(colnames(d.merged)[1:3],'cnt_age_all',colnames(d.merged)[7:dim(d.merged)[2]])

f.merged<-e.merged[,colSums(is.na(e.merged))!=nrow(e.merged) ]

touchtype<-1
merge<-as.data.frame(f.merged)

merge<-merge[!is.na(merge$dayofweek),]
merge<-merge[!is.na(merge$cnt_age_all),]
merge<-merge[!is.na(merge$part_age),]
merge<-merge[!is.na(merge$holiday),]
merge<-merge[merge$holiday==1,]

workdata<-!is.na(merge$dayofweek)&!is.na(merge$cnt_age_all)&!is.na(merge$part_age)&merge$phys_contact==as.numeric(touchtype);
workdata<-as.data.frame(merge[workdata,])

for (i in 1:length(workdata$cnt_age_all))
{ options(na.action='na.pass')
  workdata$age_group_mean[i]<-as_age_group(workdata$cnt_age_all[i], limits = as.numeric(c(1,2,5,12,15,16,25,45,65,86))) #upper limits);
}

options(na.action='na.pass')
agegrp<-(as.factor(as.matrix(workdata$age_group_mean)))
B <- as.data.frame(model.matrix(workdata$part_id ~ 0+agegrp))

day<-as.numeric(workdata$dayofweek==0|workdata$dayofweek==6) #as.numeric for TRUE/FALSE, weekends are 1
workdata<-na.omit(cbind(workdata,B,day))
#workdata$agegrp1<-rep(0,nrow(workdata))

GBfinal<-select(workdata,part_id,part_age,agegrp1,agegrp2,agegrp3,agegrp4,agegrp5,agegrp6,agegrp7,agegrp8,agegrp9,agegrp10,day)

#options(na.action='na.omit')
set<-ddply(GBfinal,c("part_id","part_age","day"),numcolwise(sum))
set<-set[complete.cases(set),];
set$part_age<-na.pass(as.numeric(as.character(set$part_age)))
set$day<-na.omit(as.numeric(as.character(set$day)))
set<-na.omit(set)
colnames(set)=c('V0','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')

newpolymod<-set[,-1]
colnames(newpolymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
final.polymod<-arrange(newpolymod,V1)

sname<-c('BEtable11', 'DEtable11', 'FItable11', 'GBtable11', 'ITtable11', 'LUtable11', 'NLtable11', 'PLtable11','FRtable11', 'ZItable11', 'PEtable11')

write.csv(final.polymod,paste0(sname[11],".csv"),row.names = FALSE)

polymod<-cbind(final.polymod,final.polymod$V12)

contact_matrix(as.matrix(polymod), age_sizes[,1], age_group_limits = age.group.limits)
contact_matrix(as.matrix(polymod), age_sizes[,1], age_group_limits = age.group.limits)*age.group.sizes
age.group.limits<<-c(1,2,5,12,15,16,25,45,65,75) #upper limits
age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)


####-------
