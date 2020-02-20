###Acceptability curve, probability each intervention has optimal cost-effectiveness
##############
accept.curve<-function(threshold.val, ICERlist)
{
  all.icer<-ICERlist
  for(qq in 1:length(all.icer))
  {
    icer<-all.icer[[qq]]
    nmb<-icer[[1]]*threshold.val-icer[[2]]
    
    ifelse(qq==1, all.nmb<-nmb,all.nmb<-cbind(all.nmb,nmb))
  }
  
  ex.zeros<-cbind(rep(0,dim(all.nmb)[1]),rep(0,dim(all.nmb)[1]),rep(0,dim(all.nmb)[1]))
  
  xx<-apply(all.nmb, 1, function(z) {which.max(z)})  
  
  for(cc in 1:length(xx))
  {ex.zeros[cc,xx[cc]]<-1}
  
  curve.data<-colMeans(ex.zeros)
  return(c(threshold.val, curve.data))
}

icer.l2<-list(unname(Reduce('+',lapply(ICER.set2[[1]], rowSums))/tab.diff),
              unname(Reduce('+',lapply(ICER.set2[[2]], rowSums))/tab.diff))


icer.l3<-list(unname(Reduce('+',lapply(ICER.set3[[1]], rowSums))/tab.diff),
              unname(Reduce('+',lapply(ICER.set3[[2]], rowSums))/tab.diff))


icer.l4<-list(unname(Reduce('+',lapply(ICER.set4[[1]], rowSums))/tab.diff),
              unname(Reduce('+',lapply(ICER.set4[[2]], rowSums))/tab.diff))

nn<-list(icer.l2, icer.l3, icer.l4)





ss<-seq(0,30000,by = 500)
for(hh in 1:length(ss)) {ifelse(hh==1, accept.out<-accept.curve(ss[hh], nn), accept.out<-rbind(accept.out,accept.curve(ss[hh], nn)))}

colnames(accept.out)<-c('WTP','Preschool','Primary','Secondary')
accept.m<-melt(data.frame(accept.out),id.vars = 'WTP')
AC.compare<-ggplot(data=accept.m, aes(x=WTP, y=value, group=variable)) +
  geom_line(aes(color=variable))+labs(title='Optimal Acceptability Curve', x='Willing-to-pay Threshold (£)', y='Probability Optimally Cost-Effective',color='Strategy')+theme_linedraw()

AC.compare+theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, .1))+scale_x_continuous(breaks = seq(0,30000, 2000))+geom_vline(xintercept = c(15000,20000), linetype=3)


ggsave(paste0('ACcompare',dcount,".png"),plot = last_plot(), width=10, units='in',device='png')

#####
#for single intervention probability it is cost effective on WTP scale
#####

greg<-lapply(tot.QALY.loss, rowSums)

greg2<-cbind(greg[[1]], greg[[2]], greg[[3]], greg[[4]])

ex.zeros<-cbind(rep(0,dim(greg2)[1]),rep(0,dim(greg2)[1]),rep(0,dim(greg2)[1]), rep(0,dim(greg2)[1]))

xx<-unname(apply(greg2, 1, function(z) {which.max(z)}))
for(cc in 1:length(xx))
{ex.zeros[cc,(xx[cc])]<-1}



qaly.l<-do.call(cbind.data.frame, lapply(lapply(tot.QALY.loss, unname), rowSums))
for(cc in 1:length(xx))
{ex.zeros[cc,xx[cc]]<-1}


pname<-c(rep('Status Quo', dim(tot.QALY.loss[[1]])[1]),rep('Preschool', dim(tot.QALY.loss[[1]])[1]),rep('Primary', dim(tot.QALY.loss[[1]])[1]),rep('Secondary', dim(tot.QALY.loss[[1]])[1]))

qaly.ls<-cbind(rowSums(qaly.l), pname)
ex.zeros<-cbind(rep(0,dim(all.nmb)[1]),rep(0,dim(all.nmb)[1]),rep(0,dim(all.nmb)[1]))
curve.data<-colMeans(ex.zeros)

prob.curve<-function(threshold.val, icer)
{
  nmb<-icer[[1]]*threshold.val-icer[[2]]
  
  ex.zeros<-rep(0,length(nmb))
  
  ex.zeros[nmb >= threshold.val]<-1  
  
  curve.data<-mean(ex.zeros)
  return(c(threshold.val, curve.data))
}

for(iii in 1:length(nn))
{
  ice.choice<-nn[[iii]]
  
  ss<-seq(0,20000,by = 250)
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out<-prob.curve(ss[hh],ice.choice), accept.out<-rbind(accept.out,prob.curve(ss[hh], ice.choice)))}
  
  ifelse(iii==1, prob.out<-accept.out, prob.out<-cbind(prob.out,accept.out[,2]))
}

colnames(prob.out)<-c('WTP','Preschool','Primary','Secondary')

#install.packages("ggthemes") # Install 
library(ggthemes) # Load

accept.m<-melt(data.frame(prob.out),id.vars = 'WTP')
AC.single<-ggplot(data=accept.m, aes(x=WTP, y=value, group=variable)) +
  geom_line(aes(color=variable))+labs(title='Acceptability Curve per Strategy', x='Willing-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy')+theme_linedraw()

AC.single+theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, .1))+scale_x_continuous(breaks = seq(1000,20000, 1000))+geom_vline(xintercept = c(15000,20000), linetype=3)


ggsave(paste0('ACsingle',dcount,".png"),plot = AC.single, width=10, units='in',device='png')

