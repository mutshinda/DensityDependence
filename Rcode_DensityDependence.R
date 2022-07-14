

###########################################################################################

##  R code for simulating data from Gompertz model                                        #															

###########################################################################################

dataGompertz<-function(n, r, k, sdr, init=k){ 
  # n is the sample size, sdr the standard deviation of the noise
  # init is the initial population size, set to k
  y<-rep(0,n)
  epsillon<-rep(0,n)
  gr<-rep(0,n)
  y[1]<-init
  for(t in 2: n){
    epsillon[t]<-rnorm(1, 0, sdr)
    gr[t]<-exp(r*(1-log(y[t-1])/k) + epsillon[t])
    y[t]<-y[t-1]*gr[t]
  }
  y
}


###########################################################################################

# OpenBUGS code for fitting the stochastic Gompertz and stochastic Ricker models to data  #

############################################################################################

gompertzModel=function() {
  for(t in 2: n){
    y[t]~dnorm(m[t],tau.y)
    m[t]<- r + beta*y[t-1]
    ypred[t]~dnorm(m[t],tau.y)
    #ypred is drawn from the PPD at time t 
    squerr[t]<-pow((ypred[t]-y[t]),2)
  }
  k<-r/(1-beta)
  r~dnorm(0,1)%_%I(0,)
  beta~dnorm(0,1)
  tau.y~dgamma(1,1)
  sigma2.y<-1/tau.y
  rmse<-mean(squerr[2:n])
  # rmse is the root mean squared error
}


rickerModel=function() {
  for(t in 2: n){
    m[t]<- y[t-1] + r*(1-exp(y[t-1])/K)
    y[t]~dnorm(m[t],tau.y)
    ypred[t]~dnorm(m[t],tau.y)
    #ypred is drawn from the PPD at time t 
    squerr[t]<-pow(ypred[t]-y[t],2)
  }
  K~dgamma(1,1)
  r~dnorm(0,1)%_%I(0,)
  tau.y~dgamma(1,1)
  sigma2.y<-1/tau.y
  rmse<-mean(squerr[2:n])
  # rmse is the root mean squared error
} 



###########################################################################################

# R code for the simulation study                                                         #

###########################################################################################

# All data are simulated from the stochastic Gompertz model

# We'll simulate 300 observations starting from k and drop the first 200 samples to ensure that the last n=100 observations come from the
# stationary distribution

## Root Mean Squared Errors under fitted stochastic Gompertz (RMSE1) and stochastic Ricker (RMSE2) models

RMSE1<-NULL

RMSE2<-NULL


## Deviance Information Criteria under fitted stochastic Gompertz (DIC1) and stochastic Ricker (DIC2) models

DIC1<-NULL 

DIC2<-NULL 


# vs=stationary variance of simulated popualtion trajectories 

# sigma2y1=estimated environmental variance from stochastic Gompertz model

# sigma2y1=estimated environmental variance from stochastic Ricker model


vs = NULL

sigma2y1 = NULL

sigma2y2 = NULL




# Fitting the Gompertz model (model1) and the Ricker model (model2) to the simulated data

# Simulation setup

m=300

k_values<-3:6; r=1.2; sdr=sqrt(0.55)

# beta=(1-r/k) # corresponding value of the AR(1) parameter


# Requiring the the R library BRugs for fitting OpenBUGS from within R


library(BRugs)



nsim=100 # number of replications for each combination of parameters


for (k in k_values){
  
  for(i in 1:nsim){ 
    
    #simulate a dataset
    
    simData<- dataGompertz(m, 1.2, k, sdr)
    
    ## stationary variance of the ith data replicate
    
    vs[i]<-var(log(simData)) 
    
    # Formatting the data for OpenBUGS
    
    DataToBUGS=list(y=log(simData[201:300]), n=100)
    
    writeModel(gompertzModel, "model1.bug")
    
    writeModel(rickerModel, "model2.bug")
    
    bugsData(DataToBUGS, "Data1.bug")
    
    
    ## Fitting the stochastic Gompertz (model1) and the stochastic Ricker (model 2) to simulated data replicates
    
    thing1=BRugsFit("model1.bug", "Data1.bug", numChains = 1, parametersToSave=c("r", "beta", "tau.y", "ypred"), nBurnin = 4000, nIter = 6000, nThin = 10,  DIC = TRUE, working.directory = NULL, digits = 3)
    
    thing2=BRugsFit("model2.bug", "Data1.bug", numChains = 1, parametersToSave=c("r",  "K", "tau.y", "ypred"), nBurnin = 4000, nIter = 6000, nThin = 4,    DIC = TRUE, working.directory = NULL, digits = 3)
    
    
    ## Predictions ypred1 from the stocahstic Gompertz model and  ypred1 from the stochastic Ricker model
    
    ypred1=thing1$Stats[4:102,1]
    
    ypred2=thing2$Stats[4:102,1]
    
    
    ## Posterior estimates (sigma2y1) of environmental variance from the stocahstic Gompertz (sigma2y2) from the stochastic Ricker model
    
    sigma2y1[i]=1/thing1$Stats[3,1]
    
    sigma2y2[i]=1/thing2$Stats[3,1]
    
    
    ## Root Mean Squared Error (RMSE1) under Gompertz model and (RMSE2) under the stochastic Ricker model 
    
    RMSE1[i]<-sqrt(mean( log(simData[202:m])-ypred1)^2)
    
    RMSE2[i]<-sqrt(mean( log(simData[202:m])-ypred2)^2)
    
    ## Deviance Information Criteria (DIC1) under Gompertz model  and (DIC2) under the Ricker model 
    
    DIC1[i]=thing1$DIC[3,3]
    
    DIC2[i]=thing2$DIC[3,3]
    
    ## Write the results to the file (here the file res_sim4.txt on my Desktop)
    
    sink("C:\\Users\\cmutshinda\\Desktop\\res_sim4.txt", append = TRUE)
    
    cat(vs[i], sigma2y1[i], sigma2y2[i], RMSE1[i], RMSE2[i], DIC1[i], DIC2[i], "\n", sep=" ")
    
    sink()
    
  }
  
}

