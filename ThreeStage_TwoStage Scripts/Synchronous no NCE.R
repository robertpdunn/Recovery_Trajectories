##Robert Dunn
## Three stage urchin model with MPA/fisheries closure simulator
## Estimated large urchin carrying capacity by setting recruitment facilitation (sigma)
## =0 and with no lobsters present. KUl is ~40 large urchins. (just as before)
## Baseline values for FU=0.1, FL=0.25
## MPA establishment means no fishing, ie FU=0,FL=0

##Can use this model to test the effects of fishing intensity and size-structured predation on recovery trajectories
##as described in Samhouri et al. 2017, Nature Eco Evo
##Factorial cross of 3 levels of fishing intensity for both pred & prey: Low, FMSY, Overfished (already conducted in Supp Info of Samhouri et al)
##For each of those 9 combinations, initiate each of 3 different "restoration" scenarios: Pred-first, Prey-first, or Synchronous
##Measure the volatility and return time (as in Samhouri et al. Fig 3) of the transient 

##The above analysis should be the baseline to verify that the model works as in Samhouri et al. Then,
##Could run the same analysis with various degrees of size-structured predation & recruitment facilitation to advance the theory

## This also differs from Samhouri et al. because they include an alternate prey item for the predator- relatively straightfwd for me to do that too

rm(list=ls())
require(deSolve)

#################################
###Kelp forest initial state#####
#################################

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestagetwostage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul))*A        #kelp dynamics (logistic growth, type I func resp)
    
    dUs = ((aM*deltaUm*Um+aL*deltaUl*Ul)*A)*(1-sigma+sigma*(Ul/Kul))-
      (gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+ (Ls* deltaLss/(1+taus*deltaLss*Us+taums*deltaLms*Um)) +MUs)*Us #small urchin dynamics
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality)
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + (Ls*deltaLms/(1+taus*deltaLss*Us+taums*deltaLms*Um)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dLs = (b*(deltaLss*Us+deltaLms*Um)/(1+taus*deltaLss*Us+taums*deltaLms*Um)- gammaL - MLs)*Ls
    #^^conversion of small lobsters to urchins via type II func resp, loss via growth & natural mortality)
    
    dL = (gammaL*Ls) + (b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)-ML-FL)*L  #spiny lobster dynamics     
    #^^conversion of lobsters to urchins via type II func resp, loss via fishing & natural mortality)
    ##Including the alternate prey availability for lobsters makes the transient much longer following fishery closure.
    return(list(c(dA, dUs, dUm, dUl, dLs, dL))) # return dn/dt as a list with each state variable as a column
  })  
}
# recruitment facilitation term:   *(1-sigma+sigma*(Ul/Kul))
#To estimate carrying capacity for large urchins, remove the recruitment facilitation term
# and run the model with no lobsters present. Equilibrium value of Ul is the new 
# value of Kul (large urchin carrying capacity).
#What about if I change other parameters- do I need to recalculate Kul each time?
###----------parameters----------------------###
state<-c(A=2459.86, Us=8.209, Um=0.684, Ul=0.108, Ls=0.0001, L=10.658)   #Unexploited equilibrium conditions (simulated elsewhere)
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters (both size classes)
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
gammaL= 0.1   #growth of small lobsters to large lobsters
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
deltaLss= 0.1   #rate of small urchin consumption by small lobsters
deltaLms= 0.075   #rate of medium urchin consumption by small lobsters
MUl=0.1     #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
MUm= 0.1    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.0     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
taums= 1e-06 #handling time of medium urchins by small lobsters
#phi= 1     #recruitment of small lobsters each time step   (removed for now)
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
MLs = 0.35   #natural mortality of small lobsters (2x that of large lobsters) 0.35 [baseline]

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run simulation to set up exploited equilibrium: A=1843.17, Us=17.33, Um=1.65, Ul=0.293, Ls=3.78, L=7.27
out=as.data.frame(lsoda(y=state, times=times, func=threestagetwostage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.00         #close urchin fishery
FL=0.0    #close lobster fishery

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestagetwostage, parms=parms))  
tail(out2,7)

#Don't need to run these
#calculate new abundance scaled to fished equilibrium
#kelp<-out2$A/out$A[tf]
#smallurchin<-out2$Us/out$Us[tf]  
#medurchin<-out2$Um/out$Um[tf]
#largeurchin<-out2$Ul/out$Ul[tf]
#smalllobster<-out2$Ls/out$Ls[tf]
#lobster<-out2$L/out$L[tf]


#time series of raw biomass
plot(out2$time,out2$A, type="l", xlab="t", ylab="N", col="forestgreen", ylim=c(0,3000), lwd=2)
plot(out2$time, out2$Us, type="l", col="purple", ylim=c(0,20), lwd=2)
lines(out2$time,out2$Um, type="l", col="blue", lwd=2)
lines(out2$time,out2$Ul, type="l", col="darkred", lwd=2)
lines(out2$time, out2$Ls, type="l", col="red", lwd=2)
lines(out2$time,out2$L, type="l", col="indianred2", lwd=2)

y<-out2
y$time<-seq(1,200)

#Calculate community volatility##
#First, sum biomass of urchins and lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$Ls+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  #11.685   #with NCE:
fishedvol<-(maxfished-y[200,9])/y[200,9] #0.203  with NCE: 
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[200,8])/y[200,8] #,  with NCE: 

#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria

################################
###Overharvested initial state##
################################
state<-c(A=2545.075, Us=6.9208, Um=0.5640, Ul=0.0971, Ls=2.6556, L=9.6092)   #Unexploited equilibrium conditions (simulated elsewhere)
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters (both size classes)
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
gammaL= 0.1   #growth of small lobsters to large lobsters
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
deltaLss= 0.1   #rate of small urchin consumption by small lobsters
deltaLms= 0.075   #rate of medium urchin consumption by small lobsters
MUl=0.1     #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
MUm= 0.1    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.250     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
taums= 1e-06 #handling time of medium urchins by small lobsters
phi= 1     #recruitment of small lobsters each time step
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.8   # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
MLs = 0.35   #natural mortality of small lobsters (2x that of large lobsters)

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run simulation to set up exploited equilibrium (overfishing) A=953.7047 Us=29.47757 Um=3.818175 Ul=0.8091771 Ls=7.899637 L=2.43718
out=as.data.frame(lsoda(y=state, times=times, func=threestagetwostage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.00         #close urchin fishery
FL=0.0    #close lobster fishery

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestagetwostage, parms=parms))  
tail(out2,7)

#time series of raw biomass
plot(out2$time,out2$A, type="l", xlab="t", ylab="N", col="forestgreen", ylim=c(0,3000), lwd=2)
plot(out2$time, out2$Us, type="l", col="purple", ylim=c(0,20), lwd=2)
lines(out2$time,out2$Um, type="l", col="blue", lwd=2)
lines(out2$time,out2$Ul, type="l", col="darkred", lwd=2)
lines(out2$time, out2$Ls, type="l", col="red", lwd=2)
lines(out2$time,out2$L, type="l", col="indianred2", lwd=2)

y<-out2
y$time<-seq(1,200)

#Calculate community volatility##
#First, sum biomass of urchins and lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$Ls+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  #11.685   #overfished: 12.463
fishedvol<-(maxfished-y[200,9])/y[200,9] #0.203     overfished:0.284 
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[200,8])/y[200,8] #,  with NCE: 

#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria