#.libPaths('/home/nwenzel/UKflu')
library(fluEvidenceSynthesis)
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(ggforce)       # kable : prettier data.frame output
library(data.table)  # MOST IMPORTANT PACKAGE SO THIS CODE ISN'T SLOW
library(abind)
library(R.utils)
library(reshape)
library(cowplot)
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
#setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/SUPPORT_R0_loaders.R')
source('Output_vacc_program_simulation.R')
#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/SUPPORT_ICER_loaders.R')
#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/IntlSens/OUTPUT_ICER_sampler.R')
vstrategy2<-dget('FUNC_cov_strategy_uk.R')

#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/OUTPUT_ICER_warp_sampler.R')
load('Intl_general_parameters.RData')

##load general parameters---------------------
p.name<-c('StatusQuo', 'Preschool', 'PrimarySchool','SecondarySchool','Preschool+Primary School'  ,'Preschool+Primary+Secondary','Preschool+Secondary','Primary+Secondary') #intervention names for saving

strain.name<-c('H1N1','H3N2','B')

load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/ili.counts.rda')
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/virological.rda')

#load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/ili2009.rda')
load('cov.function55.RData')
#load('cov.function70.RData')
#load('cov.function30.RData')

####INTL FUNCTION Start=================================================
num.samp<-200
new.cov<-mL55
r.complete<-vector('list',length(fname))

strain<-1
for(i.country in 1:(length(fname)-1))
{
  #Sampler for R0 and Re----------------------
  
  vacc.program.simulation(program = 1, i.country=i.country, strain.choice = strain, Reff=T, version='v1')
  vacc.program.simulation(program = 2, i.country=i.country, strain.choice = strain,Reff=T,version='v1')
  vacc.program.simulation(program = 3, i.country=i.country, strain.choice = strain, Reff=T, version='v1')
  vacc.program.simulation(program = 4, i.country=i.country, strain.choice = strain, Reff=T, version='v1')
  
  vacc.program.simulation(program = 1, i.country=i.country, strain.choice = strain, Reff=F,version='v1')
  vacc.program.simulation(program = 2, i.country=i.country, strain.choice = strain, Reff=F,version='v1')
  vacc.program.simulation(program = 3, i.country=i.country, strain.choice = strain, Reff=F,version='v1')
  vacc.program.simulation(program = 4, i.country=i.country, strain.choice = strain, Reff=F,version='v1')
}#end country loop

