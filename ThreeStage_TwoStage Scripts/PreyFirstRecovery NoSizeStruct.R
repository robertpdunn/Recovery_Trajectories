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

#############################
###Kelp forest initial state#
#############################

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
    
    dLs = phi + (b*(deltaLss*Us+deltaLms*Um)/(1+taus*deltaLss*Us+taums*deltaLms*Um)- gammaL - MLs)*Ls
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
#state<-c()  #Unexploited equilibrium conditions no NCE
state<-c(A=2315.215, Us= 10.52901, Um= 0.8327827, Ul= 0.05129015, Ls= 2.658525, L= 11.45617)   #Unexploited equilibrium, NoSizeStructure
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2  #rate of kelp consumption by small urchins baseline:0.2
deltaUm= 0.2   #rate  of kelp consumption by med. urchins baseline 0.2
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters (both size classes)
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
gammaL= 0.1   #growth of small lobsters to large lobsters
deltaLs= 0.133  #rate of small urchin consumption by lobsters, 0.2
deltaLm=0.133  #rate of med urchin consumption by lobsters, 0.15
deltaLl=0.133   #rate of large urchin consumption by lobsters, 0.05
deltaLss= 0.133   #rate of small urchin consumption by small lobsters, 0.1
deltaLms= 0.133   #rate of medium urchin consumption by small lobsters, 0.075
MUl=0.1     #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
MUm= 0.1    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
taums= 1e-06 #handling time of medium urchins by small lobsters
phi= 1     #recruitment of small lobsters each time step
ML= 0.175      # natural mortality of lobsters: 0.175 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
MLs = 0.35   #natural mortality of small lobsters (2x that of large lobsters)

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run simulation for baseline fishing mortality to both species
out=as.data.frame(lsoda(y=state, times=times, func=threestagetwostage, parms=parms))  #Multi-species exploited equilibrium
tail(out,7)
#basic plot
#plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen", main="3 urchin stages", ylim=c(0,100))
#lines(out$time, out$Us, type="l", col="purple")
#lines(out$time,out$Um, type="l", col="blue")
#lines(out$time,out$Ul, type="l", col="darkred")
#lines(out$time, out$Ls, type="l", col="red")
#lines(out$time,out$L, type="l", col="indianred2")

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.0         #close the urchin fishery

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestagetwostage, parms=parms))  
tail(out2,7)

#calculate new abundance scaled to fished equilibrium
#kelp<-out2$A/out$A[tf]
#smallurchin<-out2$Us/out$Us[tf]  
#medurchin<-out2$Um/out$Um[tf]
#largeurchin<-out2$Ul/out$Ul[tf]
#smalllobster<-out2$Ls/out$Ls[tf]
#lobster<-out2$L/out$L[tf]

##Time series of abundance relative to pre-MPA equilibria N values
#plot(out2$time,kelp, type="l", lwd=2,xlab="t (years since fishery closure)",
 #    ylab="Relative Abundance (Nt/N0)", col="forestgreen",ylim=c(0,2),
  #   main="Post-fishery closure dynamics: \n Pred Recovery no NCE", bty="L",las=1)
#lines(out2$time,smallurchin, type="l", lwd=2, col="purple")
#lines(out2$time,medurchin, type="l", lwd=2, col="blue")
#lines(out2$time,largeurchin, type="l", lwd=2, col="darkred")
#lines(out2$time,smalllobster, type="l", lwd=2, col="red")
#lines(out2$time,lobster, type="l", lwd=2, col="indianred2")

#Now stop fishing for lobsters
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200],Ls = out2$Ls[200], L=out2$L[200])
FL=0.0         #close the lobster fishery

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestagetwostage, parms=parms))  
tail(out3,7)

#calculate new abundance scaled to fished equilibrium
#kelp2<-out3$A/out$A[tf]
#smallurchin2<-out3$Us/out$Us[tf]  
#medurchin2<-out3$Um/out$Um[tf]
#largeurchin2<-out3$Ul/out$Ul[tf]
#smalllobster2<-out3$Ls/out$Ls[tf]
#lobster2<-out3$L/out$L[tf]

y<-rbind(out2,out3)
y$time<-seq(1,400)
plot(y$time,y$A, ylim=c(0,3000), type="l", lwd=2, col="forestgreen")
plot(y$time, y$Us, ylim=c(0,30), type="l", lwd=2, col="purple")
lines(y$Um, type="l", lwd=2, col="blue")
lines(y$Ul, type="l", lwd=2, col="darkred")
lines(y$Ls, type="l", lwd=2, col="red")
lines(y$L, type="l", lwd=2, col="indianred2")

#g1<-cbind(kelp,smallurchin,medurchin,largeurchin,smalllobster,lobster)
#g2<-cbind(kelp2,smallurchin2,medurchin2,largeurchin2,smalllobster2,lobster2)
#g<-rbind(g1,g2)

##Time series of abundance relative to pre-MPA equilibria N values
#pdf("PreyNoSizeStruc.pdf")
#plot(y$time,g[,1], type="l", lwd=2,xlab="t (years since fishery closure)",
 #    ylab="Relative Abundance (Nt/N0)", col="forestgreen",ylim=c(0,2),
  #   main="Prey first recovery: with NCE", bty="L",las=1)
#lines(y$time,g[,2], type="l", lwd=2, col="purple")
#lines(y$time,g[,3], type="l", lwd=2, col="blue")
#lines(y$time,g[,4], type="l", lwd=2, col="darkred")
#lines(y$time,g[,5], type="l", lwd=2, col="red")
#lines(y$time,g[,6], type="l", lwd=2, col="indianred2")
#dev.off()

#Extract final abundances (relative to fished equilibrium), record in excel, is same for all 3 recovery trajectory scenarios
#tail(g)

 
#Calculate community volatility##
#First, sum biomass of large urchins and large lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$Ls+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens) #9.815
fishedvol<-maxfished-y[400,9] # 0.03
preyNoSize<-y$FishedDens
maxcomm<-max(y$CommDens) #33.229
commvol<-maxcomm-y[400,8] #9.606

################
#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria




###############################
##Urchin barren initial state##    #These were nearly identical to kelp forest initial conditions
###############################

state<-c(A=1, Us=100, Um=100, Ul=100, Ls=5, L=5)   #Initial state of the system
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
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
taums= 1e-06 #handling time of medium urchins by small lobsters
phi= 1     #recruitment of small lobsters each time step
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
MLs = 0.35   #natural mortality of small lobsters (2x that of large lobsters)

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run simulation for baseline fishing mortality to both species
out=as.data.frame(lsoda(y=state, times=times, func=threestagetwostage, parms=parms))  
tail(out,7)
#basic plot
#plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen", main="3 urchin stages", ylim=c(0,100))
#lines(out$time, out$Us, type="l", col="purple")
#lines(out$time,out$Um, type="l", col="blue")
#lines(out$time,out$Ul, type="l", col="darkred")
#lines(out$time, out$Ls, type="l", col="red")
#lines(out$time,out$L, type="l", col="indianred2")
#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.0         #close the urchin fishery

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestagetwostage, parms=parms))  
tail(out2,7)

#calculate new abundance scaled to fished equilibrium
kelp<-out2$A/out$A[tf]
smallurchin<-out2$Us/out$Us[tf]  
medurchin<-out2$Um/out$Um[tf]
largeurchin<-out2$Ul/out$Ul[tf]
smalllobster<-out2$Ls/out$Ls[tf]
lobster<-out2$L/out$L[tf]

##Time series of abundance relative to pre-MPA equilibria N values
#plot(out2$time,kelp, type="l", lwd=2,xlab="t (years since fishery closure)",
#    ylab="Relative Abundance (Nt/N0)", col="forestgreen",ylim=c(0,2),
#   main="Post-fishery closure dynamics: \n Pred Recovery no NCE", bty="L",las=1)
#lines(out2$time,smallurchin, type="l", lwd=2, col="purple")
#lines(out2$time,medurchin, type="l", lwd=2, col="blue")
#lines(out2$time,largeurchin, type="l", lwd=2, col="darkred")
#lines(out2$time,smalllobster, type="l", lwd=2, col="red")
#lines(out2$time,lobster, type="l", lwd=2, col="indianred2")

#Now stop fishing for lobsters, institute non-consumptive effect on urchin grazing.
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200],Ls = out2$Ls[200], L=out2$L[200])
FL=0.0         #close the lobster fishery
deltaUs= 0.1   #reduced rate of kelp consumption by small urchins due to more lobsters, ie non-consumptive effect (tmi), originally 0.2
deltaUm= 0.1   #reduced rate  of kelp consumption by med. urchins, originally 0.2

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, gammaL=gammaL,
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, taums=taums, 
        phi=phi, ML=ML, FL=FL, MLs=MLs)  #run model again, same parameters except no fishing, with new initial conditions
out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestagetwostage, parms=parms))  
tail(out3,7)

#calculate new abundance scaled to fished equilibrium
kelp2<-out3$A/out$A[tf]
smallurchin2<-out3$Us/out$Us[tf]  
medurchin2<-out3$Um/out$Um[tf]
largeurchin2<-out3$Ul/out$Ul[tf]
smalllobster2<-out3$Ls/out$Ls[tf]
lobster2<-out3$L/out$L[tf]

y<-rbind(out2,out3)
y$time<-seq(1,400)
plot(y$time,y$A, ylim=c(0,3000), type="l", lwd=2, col="forestgreen")
plot(y$time, y$Us, ylim=c(0,20), type="l", lwd=2, col="purple")
lines(y$Um, type="l", lwd=2, col="blue")
lines(y$Ul, type="l", lwd=2, col="darkred")
lines(y$Ls, type="l", lwd=2, col="red")
lines(y$L, type="l", lwd=2, col="indianred2")

g1<-cbind(kelp,smallurchin,medurchin,largeurchin,smalllobster,lobster)
g2<-cbind(kelp2,smallurchin2,medurchin2,largeurchin2,smalllobster2,lobster2)
g<-rbind(g1,g2)

##Time series of abundance relative to pre-MPA equilibria N values
plot(y$time,g[,1], type="l", lwd=2,xlab="t (years since fishery closure)",
     ylab="Relative Abundance (Nt/N0)", col="forestgreen",ylim=c(0,2),
     main="Post-fishery closure dynamics: \n Prey Recovery with NCE", bty="L",las=1)
lines(y$time,g[,2], type="l", lwd=2, col="purple")
lines(y$time,g[,3], type="l", lwd=2, col="blue")
lines(y$time,g[,4], type="l", lwd=2, col="darkred")
lines(y$time,g[,5], type="l", lwd=2, col="red")
lines(y$time,g[,6], type="l", lwd=2, col="indianred2")