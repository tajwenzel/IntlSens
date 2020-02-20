#parameters draws per country for all interventions
#not dependent on program



SENS.samp.pars<-function(strain, i.country, num.samp, version)
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
    
    iteration<-c(round(runif(total.samp,dim(mcmc.post)[1]*0.85,dim(mcmc.post)[1]))) 
    #take only last 20% of the runs to use for posterior distribution for calculations
    
    post.sample[[i.season]]<<-mcmc.post[iteration,] 
    rand.contact.samp[[i.season]]<<-contact.ids[iteration]
    #remove(jj)
    
  }
  post.sample<-list()
  rand.contact.samp<-list()
  total.samp<-num.samp
  #remove(jj)
  if(i.country==4)
  {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
    flu.tablesB<-list.files(pattern=glob2rx('GBflu.u2*B'));
    flu.tablesH1<-list.files(pattern=glob2rx('GBflu.u2*H1'));
    flu.tablesH3<-list.files(pattern=glob2rx('GBflu.u2*H3'));
    
    ct.tablesB<-list.files(pattern=glob2rx('GBct.u2*B'))
    ct.tablesH1<-list.files(pattern=glob2rx('GBct.u2*H1'))
    ct.tablesH3<-list.files(pattern=glob2rx('GBct.u2*H3'))
  }else{
  if(i.country==3)
  {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
    ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.r2*H1', sep='')))
    ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.r2*H3', sep='')))
    ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.r2*B', sep='')))
    
    flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.r2*B',sep='')));
    flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.r2*H1',sep='')));
    flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.r2*H3',sep='')));
  } else {
  if(i.country==9)
  {setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
    ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.r1*H1', sep='')))
    ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.r1*H3', sep='')))
    
    flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.r1*H1',sep='')));
    flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.r1*H3',sep='')));
  }
  else {
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
    
    ct.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.u1*H1', sep='')))
    ct.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'ct.u1*H3', sep='')))
    ct.tablesB<- list.files(pattern=glob2rx(paste(fname[i.country],'ct.u1*B', sep='')))
    
    flu.tablesB<-list.files(pattern=glob2rx(paste(fname[i.country], 'flu.u1*B',sep='')));
    flu.tablesH1<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.u1*H1',sep='')));
    flu.tablesH3<-list.files(pattern=glob2rx(paste(fname[i.country],'flu.u1*H3',sep='')));
  }}}
  
if(strain==1){year.length<-length(ct.tablesH1)}
if(strain==2){year.length<-length(ct.tablesH3)}
if(strain==3){year.length<-length(ct.tablesB)}
  
#jj<-where('year.length')
sapply(1:year.length, FUN= function(x) mcmc.samp(x, num.samp))
#rm(i.country)
p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving
  
  ###Save output to file
  kname<-paste0(fname[i.country],'PostSample',version, strain.name[strain])  
  save(post.sample,file=kname) 
  jname<-paste0(fname[i.country],'PostContactSample',version,strain.name[strain])  
  save(rand.contact.samp, file=jname)
}

