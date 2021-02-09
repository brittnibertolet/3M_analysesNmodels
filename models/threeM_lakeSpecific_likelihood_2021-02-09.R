# Lake specific model 
## 1st-order reactions from bulk OM to acetate and carbon dioxide to CH4 
## 7 parameters
lakespecific.nll = function(param=parameters, data=ydata){
  
  # parameters are estimate across all lakes, so not inside the for loop
  
  # parameters
  r.ba = exp(param[1]) # conversion of bulk C to acetate
  r.bc = exp(param[2]) # conversion of bulk C to CO2
  r.aa = exp(param[3]) # conversion of bulk C to acetate
  r.ac = exp(param[4]) # conversion of bulk C to CO2
  u.a  = exp(param[5]) # conversion acetate to methane
  u.c  = exp(param[6]) # conversion CO2 to methane
  
  std.m = exp(param[7])
  std.c = exp(param[8])
  std.a = exp(param[9])
  
  # Get lake ID
  lakei=ydata$lakeID
  
  # CH4 observations
  m = data$datam[names(data$datam)==lakei]
  # CO2 observations 
  c = data$datac[names(data$datac)==lakei]
  # acetate observations
  a = data$dataa[names(data$dataa)==lakei]
    
  # predictor variables
  bulkC = data$OMdata$OM_umol[data$OMdata$LakeID==lakei] # umol of carbon
  BES = c(1,1,1,0,0,0,0,0,0,1,1,1) # BES binary indicator
  algae = c(1,1,1,1,1,1,0,0,0,0,0,0) # algae binary indicator
  Calgae=(25/12)/1000*1000000/2 # algal treatment amount C
    
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
  ahat[1,]=ydata$init.a$acetate_umol[ydata$init.a$lakeID==lakei]
  
  
  # simulate process model 
  tp=nrow(yhat)
  # simulate in a matrix so that there are 12 columns (one for every bottle)
  for(i in 2:tp){
    ahat[i,] = ahat[i-1,] + r.ba*bulkC + algae*r.aa*Calgae - BES*u.a*ahat[i-1,]
    chat[i,] = chat[i-1,] + r.bc*bulkC + algae*r.ac*Calgae - BES*u.c*chat[i-1,]
    yhat[i,] = yhat[i-1,] + BES*u.a*ahat[i-1,] + BES*u.c*chat[i-1,]
  }  
  
  # Calculate log likelihoods
  # make sure data and estimates are lining up properly 
  # how many NAs in the dnorm call and is that forcing a low nll 
  ## if there are lots of NA, artificially make nll very high
  nllm = -sum(dnorm(x = as.matrix(m[[1]][-1,-1:-2]), mean = yhat[-1,], sd = std.m, log = T), na.rm=T)
  nllc = -sum(dnorm(x = as.matrix(c[[1]][-1,-1:-2]), mean = chat[-1,], sd = std.c, log = T), na.rm=T)
  nlla = -sum(dnorm(x = as.matrix(a[[1]][-1]), mean = ahat[nrow(ahat),], sd = std.a, log = T), na.rm=T)
  
  # sum the loglikelihoods and weight by the number of observations 
  nll = sum(c((nllm/length(nllm)),(nllc/length(nllm)),(nlla/length(nllm)))) 
  return(nll)
}
