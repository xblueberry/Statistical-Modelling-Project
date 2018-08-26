#Statistical Modelling: exam project
#Author: Daniel Izquierdo Juncas
#June 2017

#Part A
#install.packages("http://stat.wharton.upenn.edu/~buja/PAPERS/PoSI_1.0.tar.gz", repos=NULL, type="source") 
#install.packages("mnormt")

library("PoSI")
library("mnormt") # to generate from a multivariate normal distribution: function rmnorm

datagen=function(seed){
  
  theta= c(3,3,3,0,0)
  N=200         #Sample size
  q=length(theta)    # length of parameter vector
  p0=1       #number of parameters which are in all models. Here it is 1 which corresponds to having an intercept in all models 
  pselected=0 #Stores the selected model in each data generation process 
  sig2=1        #Variance of the error (standard normal distribution)
  
  
  mun=as.vector(rep(0,(q-1)))  # Mean vector of your covariates. 
  rho=  0.75 #0.25  # Off diagonal elements of the varance covariance matrix of covariates. 
  Xsigma=rho*(matrix(1,q-1,q-1)-diag(q-1))+diag(q-1) #Variance covariance matrix of covariates
  
  set.seed(seed)      # Seed number such that you track your simulation. 
  
  X = cbind(rep(1,N),rmnorm(N,mean=mun,varcov=Xsigma))   # generated covariates in a matrix
  y = X%*%theta + rnorm(N,mean=0,sd=sqrt(sig2))     # Response variable in a vector
  outdgf=list(X,y)
  return(outdgf)
}

#this function generates the set of models, each row corresponds to a model. 1 means the presence of a covariate and 0 means absence.
combinations = function(n){
  comb = NULL
  for( i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
  return(comb)
}
Pi1=combinations(4)
p0=1
Pi=cbind(matrix(1,nrow=nrow(Pi1),ncol=p0),Pi1)   # Intercept is present in all models.
pfocus=2    # The model that you should consider for A.2. Change this for A.3. 
M=nrow(Pi)  # Number of candidate models in which you do the selection
pselected=0

conddatagen=function(seed,pfocus=2){
  while(TRUE){
    
    if(pselected==pfocus) break()
    
    else{
      
      set.seed(seed)      # Seed number such that you track your simulation. 
      
      #aic=rep(NA,M)       
      aic<-c()  # Generate a vector to store AIC (BIC) scores for each model
      datasim=datagen(seed)
      X=datasim[[1]]
      y=datasim[[2]]

      # Calculate AIC for the candidate models   
      for (i in 1:M){
        X2 <- c()
        for (j in 1:5){
          if(Pi[i,j]==1) {
            X2<-cbind(X2,X[,j])
          }
        }
        fit<-lm(y~X2+0)
        aic<-cbind(aic,AIC(fit))
      }
      #aic

      pselected=which.min(aic)# Which model is selected by AIC?
      
    }
    # If the above pselected is equal to the considered model, pfocus, you actually have a data set in which the selected model is equal to the focus model. 
    # If the above pselected is not equal to the focus model, you change the seed and generate another data set. Here, you should add the current seed number with a number bigger than 5000.
    seed=seed+5678
    
  }
  outcdgf=list(X,y)
  return(outcdgf)
}

# This while loop helps you to generate data such that the selected model based on AIC (BIC) is equal to the considered model, psel. 

simulationCI=function(seed,pfocus=2){
  
  conddata=conddatagen(seed,pfocus=pfocus)
  X=conddata[[1]]
  y=conddata[[2]]
  for (i in 1:M){
    X. <- c()
    for (j in 1:5){
      if(Pi[i,j]==1) {
        X.<-cbind(X.,X[,j])
      }
    }
    if (pfocus==i){regsel=lm(y~X.+0)} # Fit the selected model on the data
  }
  thetaest= regsel$coefficients	     # Get the estimates from the fitted model	
  naiveCI= cbind(matrix(thetaest-qnorm(0.975)*(sqrt(diag(vcov(regsel)))),ncol=1),matrix(thetaest+qnorm(0.975)*(sqrt(diag(vcov(regsel)))),ncol=1))  #Naive CI
    
  Posiconstant= summary(PoSI(conddata[[1]][,-1]))[1]
  
  posiCI= cbind(matrix(thetaest-Posiconstant*(sqrt(diag(vcov(regsel)))),ncol=1),matrix(thetaest+Posiconstant*(sqrt(diag(vcov(regsel)))),ncol=1)) #PoSI CI
  
  write.table(t(c(thetaest,c(t(naiveCI)),c(t(posiCI)))),file=paste("CI",".txt",sep=""),append=T,col.names=FALSE,row.names=FALSE)

# check for a single simulation run which number is which in the output file! You need this later.
}

#for(seed in 1:1000){
#  simulationCI(seed)
#}

myseed=654210 
i=0
while (i<1000){
  simulationCI(myseed,2)
  myseed=myseed+5000
  i=i+1
}


#Work with results:
results=read.table(file= "Daniel Izquierdo (A2,rho=0.75).txt" ) #same for A3,rho=0.75
Thetaests=matrix(NA,ncol=4,nrow=1000)
Thetaests=results[,1:4]

# empirical confidence intervals:
empiricalCI<-rbind(c(quantile(Thetaests[,1],p=0.025),quantile(Thetaests[,1],p=0.975)),
      c(quantile(Thetaests[,2],p=0.025),quantile(Thetaests[,2],p=0.975)),
      c(quantile(Thetaests[,3],p=0.025),quantile(Thetaests[,3],p=0.975)),
      c(quantile(Thetaests[,4],p=0.025),quantile(Thetaests[,4],p=0.975)))
rownames(empiricalCI)<-c("Parameter 1","Parameter 2","Parameter 3","Parameter 4")
empiricalCI

# count how many confidence intervals contain the true parameter (i.e., lower bound of interval smaller than true value and upper bound of interval larger than true value).
# Do this for naive intervals, Posi intervals, and for each parameter.
i=0
NaiveCI.param1.count=0
for (i in 1:1000){
  if(results[i,5]<3 & results[i,6]>3 ){NaiveCI.param1.count=NaiveCI.param1.count+1}
  i=i+1
}

i=0
NaiveCI.param2.count=0
for (i in 1:1000){
  if(results[i,7]<3 & results[i,8]>3 ){NaiveCI.param2.count=NaiveCI.param2.count+1}
  i=i+1
}

i=0
NaiveCI.param3.count=0
for (i in 1:1000){
  if(results[i,9]<3 & results[i,10]>3 ){NaiveCI.param3.count=NaiveCI.param3.count+1}
  i=i+1
}

i=0
NaiveCI.param4.count=0
for (i in 1:1000){
  if(results[i,11]<0 & results[i,12]>0 ){NaiveCI.param4.count=NaiveCI.param4.count+1}
  i=i+1
}

i=0
posiCI.param1.count=0
for (i in 1:1000){
  if(results[i,13]<3 & results[i,14]>3 ){posiCI.param1.count=posiCI.param1.count+1}
  i=i+1
}

i=0
posiCI.param2.count=0
for (i in 1:1000){
  if(results[i,15]<3 & results[i,16]>3 ){posiCI.param2.count=posiCI.param2.count+1}
  i=i+1
}

i=0
posiCI.param3.count=0
for (i in 1:1000){
  if(results[i,17]<3 & results[i,18]>3 ){posiCI.param3.count=posiCI.param3.count+1}
  i=i+1
}

i=0
posiCI.param4.count=0
for (i in 1:1000){
  if(results[i,19]<0 & results[i,20]>0 ){posiCI.param4.count=posiCI.param4.count+1}
  i=i+1
}

coverage.prob<-matrix(c(NaiveCI.param1.count,posiCI.param1.count,NaiveCI.param2.count,
                        posiCI.param2.count,NaiveCI.param3.count,posiCI.param3.count,
                        NaiveCI.param4.count,posiCI.param4.count),nrow=2,ncol=4)
colnames(coverage.prob)<-c("Parameter 1","Parameter 2","Parameter 3","Parameter 4")
rownames(coverage.prob)<-c("Naive CI","PoSI CI")
coverage.prob/1000

# construct average confidence interval: colMeans(...)
NaiveCI.param1=matrix(NA,ncol=2,nrow=1000)
NaiveCI.param1=results[,5:6]
colMeans(NaiveCI.param1)
# now similar for Posi and for other parameters.
NaiveCI.param2=matrix(NA,ncol=2,nrow=1000)
NaiveCI.param2=results[,7:8]

NaiveCI.param3=matrix(NA,ncol=2,nrow=1000)
NaiveCI.param3=results[,9:10]

NaiveCI.param4=matrix(NA,ncol=2,nrow=1000)
NaiveCI.param4=results[,11:12]

posiCI.param1=matrix(NA,ncol=2,nrow=1000)
posiCI.param1=results[,13:14]

posiCI.param2=matrix(NA,ncol=2,nrow=1000)
posiCI.param2=results[,15:16]

posiCI.param3=matrix(NA,ncol=2,nrow=1000)
posiCI.param3=results[,17:18]

posiCI.param4=matrix(NA,ncol=2,nrow=1000)
posiCI.param4=results[,19:20]


averageCI<-rbind(colMeans(NaiveCI.param1),colMeans(NaiveCI.param2),colMeans(NaiveCI.param3),
                 colMeans(NaiveCI.param4),colMeans(posiCI.param1),colMeans(posiCI.param2),
                 colMeans(posiCI.param3),colMeans(posiCI.param4))
colnames(averageCI)<-c("2.5%","97.5%")
rownames(averageCI)<-c("Naive: Parameter 1","Naive: Parameter 2","Naive: Parameter 3","Naive: Parameter 4",
                       "PoSI: Parameter 1","PoSI: Parameter 2","PoSI: Parameter 3","PoSI: Parameter 4")
averageCI
