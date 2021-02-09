rm(list=ls())

# Load libraries and ggplot formatting
library(reshape2)
source("code/setTheme_BB.R")
# Load input data for model
ydata=readRDS("data/inputdataList-current.RDS")
# Load custom function for estimating likelihood
source("models/threeM_lakeSpecific_likelihood_2021-02-09.R")
# Load data for comparing to simulation
gc=read.csv("data/3M_gcData_forsimulations.csv", stringsAsFactors = F)

#### Estimate lake specific parameters ####

# Get all lakeIDs and set up output vectors
lakes=names(ydata$datam)
r.ba = rep(NA, length=length(lakes)) # conversion of bulk C to acetate
r.bc = rep(NA, length=length(lakes)) # conversion of bulk C to CO2
r.aa = rep(NA, length=length(lakes)) # conversion of bulk C to acetate
r.ac = rep(NA, length=length(lakes)) # conversion of bulk C to CO2
u.a  = rep(NA, length=length(lakes)) # conversion acetate to methane
u.c  = rep(NA, length=length(lakes)) # conversion CO2 to methane
std.m = rep(NA, length=length(lakes))
std.c = rep(NA, length=length(lakes))
std.a = rep(NA, length=length(lakes))
nll = rep(NA, length=length(lakes))
converge = rep(NA, length=length(lakes))

# Loop through every lake
for(j in 1:length(lakes)){
  # put lakeID into input data
  ydata$lakeID=lakes[j]
  # parameter guesses
  parameters=log(c(0.005, 0.001, 0.001, 0.005, 0.001, 0.005,10, 10, 10)) # Parameter guesses
  # fit model
  fit1=optim(parameters, fn=lakespecific.nll, ydata, control = list(maxit=1e4)) # optimizer
  # put outputs into place
  r.ba[j]=fit1$par[1]
  r.bc[j]=fit1$par[2]
  r.aa[j]=fit1$par[3]
  r.ac[j]=fit1$par[4]
  u.a[j]=fit1$par[5]
  u.c[j]=fit1$par[6]
  std.m[j]=fit1$par[7]
  std.c[j]=fit1$par[7]
  std.a[j]=fit1$par[7]
  nll[j]=fit1$value
  converge[j]=fit1$convergence
}
# combine parameter estimates into a single dataframe 
lakefits=data.frame(lakes, r.ba, r.bc, r.aa, r.ac, u.a, u.c, std.m, std.c, std.a, nll, converge)
# get rid of individual parameter vectors
rm(r.ba, r.bc, r.aa, r.ac, u.a, u.c, std.m, std.c, std.a, nll, converge)

# simulate
simulateModel1=function(data=ydata){
  library(reshape2)
  CH4simData=data.frame()
  CO2simData=data.frame()

  for(i in 1:length(lakes)){
    # CH4 observations
    m = data$datam[names(data$datam)==lakes[i]]
    # CO2 observations 
    c = data$datac[names(data$datac)==lakes[i]]
    # acetate observations
    a = data$dataa[names(data$dataa)==lakes[i]]
    
    # predictor variables
    bulkC = data$OMdata$OM_umol[data$OMdata$LakeID==lakes[i]] # umol of carbon
    BES = c(1,1,1,0,0,0,0,0,0,1,1,1) # BES binary indicator
    algae = c(1,1,1,1,1,1,0,0,0,0,0,0) # algae binary indicator
    Calgae=(25/12)/1000*1000000 # algal treatment amount C
    
    # set up yhat
    yhat=matrix(NA, ncol=12, nrow=nrow(m[[1]]))
    colnames(yhat)=colnames(m[[1]])[-1:-2]
    # initial observation of y
    yhat[1,]=t(m[[1]][1,-1:-2]) 
    
    # chat
    chat=matrix(NA, ncol=12, nrow=nrow(c[[1]]))
    colnames(chat)=colnames(c[[1]])[-1:-2]
    # initial observation of c
    chat[1,]=t(c[[1]][1,-1:-2]) 
    
    # set up ahat
    ahat=matrix(NA, ncol=12, nrow=nrow(m[[1]]))
    colnames(ahat)=colnames(a[[1]])[-1]
    # initial observation of a
    ahat[1,]=ydata$init.a$acetate_umol[ydata$init.a$lakeID==lakes[i]]
    
    # parameters 
    r.ba = lakefits$r.ba[lakefits$lakes==lakes[i]] # conversion of bulk C to acetate
    r.bc = lakefits$r.bc[lakefits$lakes==lakes[i]] # conversion of bulk C to CO2
    r.aa = lakefits$r.aa[lakefits$lakes==lakes[i]] # conversion of algae to acetate
    r.ac = lakefits$r.ac[lakefits$lakes==lakes[i]] # conversion of algae to CO2
    u.a = lakefits$u.a[lakefits$lakes==lakes[i]] # conversion acetate to methane
    u.c = lakefits$u.c[lakefits$lakes==lakes[i]] # conversion CO2 to methane
    
    std.m = lakefits$std.m[lakefits$lakes==lakes[i]]
    std.c = lakefits$std.c[lakefits$lakes==lakes[i]]
    std.a = lakefits$std.a[lakefits$lakes==lakes[i]]
    
    # simulate process model 
    tp=nrow(yhat)
    # simulate in a matrix so that there are 12 columns (one for every bottle)
    for(i in 2:tp){
      ahat[i,] = ahat[i-1,] + r.ba*bulkC + algae*r.aa*Calgae - BES*u.a*ahat[i-1,]
      chat[i,] = chat[i-1,] + r.bc*bulkC + algae*r.ac*Calgae - BES*u.c*chat[i-1,]
      yhat[i,] = yhat[i-1,] + BES*u.a*ahat[i-1,] + BES*u.c*chat[i-1,]
    }
    
    # Get methane data
    temp=data.frame(yhat)
    temp$incub_days=m[[1]]$incub_days
    temp$lakeID=m[[1]]$lakeID
    temp=melt(temp, id.vars = c("lakeID", "incub_days"), variable.name = "treatment", value.name = "sim1CH4_umol")
    CH4simData=rbind(CH4simData, temp)
    
    # Get CO2 data
    temp=data.frame(chat)
    temp$incub_days=m[[1]]$incub_days
    temp$lakeID=m[[1]]$lakeID
    temp=melt(temp, id.vars = c("lakeID", "incub_days"), variable.name = "treatment", value.name = "sim1CO2_umol")
    CO2simData=rbind(CO2simData, temp)
  }
  
  CH4simData$rep=gsub(".*_", "", CH4simData$treatment,ignore.case = T)
  CH4simData$treatment=gsub("_[0-9]", "", CH4simData$treatment,ignore.case = T)
  CO2simData$rep=gsub(".*_", "", CO2simData$treatment,ignore.case = T)
  CO2simData$treatment=gsub("_[0-9]", "", CO2simData$treatment,ignore.case = T)
  
  simData=merge(CH4simData, CO2simData, by=c("lakeID", "treatment", "incub_days", "rep"))
  return(simData)
}
simData=simulateModel1()

# merge with data
sim=merge(simData, gc, by=c("lakeID", "treatment", "incub_days", "rep"))

ggplot(sim, aes(x=incub_days, y=sim1CH4_umol))+
  geom_line(aes(group=rep),color="red")+
  facet_grid(treatment~lakeID)+
  geom_point(aes(x=incub_days, y=totalCH4_umol), color="black")+
  ylab("Methane (umol)")+xlab("Day")

ggplot(sim, aes(x=incub_days, y=sim1CO2_umol))+
  geom_line(aes(group=rep),color="red")+
  facet_grid(treatment~lakeID)+
  geom_point(aes(x=incub_days, y=totalCO2_umol), color="black")+
  ylab("Methane (umol)")+xlab("Day")






