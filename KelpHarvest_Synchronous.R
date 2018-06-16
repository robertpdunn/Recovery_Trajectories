##Robert Dunn
## Three stage urchin model
##With KELP HARVEST

require(deSolve)
require(lattice)
rm(list=ls())

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
kelpharvest=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul)-Fk)*A        #kelp dynamics (logistic growth, type I func resp, harvest)
    
    dUs = ((aM*deltaUm*Um+aL*deltaUl*Ul)*A)*(1-sigma+sigma*(Ul/Kul))-(gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+MUs)*Us #small urch dyn
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality)
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dL = ((b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))-ML-FL)*L          #spiny lobster dynamics     
    #^^conversion of lobsters to urchins via type II func resp, loss via fishing & natural mortality)
    return(list(c(dA, dUs, dUm, dUl, dL))) # return dn/dt as a list with each state variable as a column
  })  
}
#############################################
###SYNCHRONOUS RECOVERY######################
#############################################
###----------parameters----------------------###
state<-c(A=2459.864, Us=8.209526, Um=0.6845794, Ul=0.1081592, L=10.65874) #Unexploited equilibrium (no fishing), inits for trajectory sims
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Fk=0.1        #kelp harvest mortality
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2  ##rate of small urchin consumption by lobsters 
deltaLm=0.15 ##rate of med urchin consumption by lobsters  
deltaLl=0.05 ##rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment   #0.75 for barrens inits, #0.6
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model to set up harvested equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=kelpharvest, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use harvested equilibrium population values from above
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
Fk=0.0        #close kelp harvest
FU=0.00       #close urchin fishery
FL=0.0        #close lobster fishery

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out2=as.data.frame(lsoda(y=newstate, times=times, func=kelpharvest, parms=parms))  
tail(out2,7)

#time series of raw biomass
plot(out2$time,out2$A, type="l", xlab="t", ylab="N", col="forestgreen", ylim=c(0,3000), lwd=2)
plot(out2$time, out2$Us, type="l", col="purple", ylim=c(0,20), lwd=2)
lines(out2$time,out2$Um, type="l", col="blue", lwd=2)
lines(out2$time,out2$Ul, type="l", col="darkred", lwd=2)
lines(out2$time,out2$L, type="l", col="indianred2", lwd=2)

y<-out2
y$time<-seq(1,200)

kelp<-(y$A/out$A[200])   #scale kelp, large urchins & lobsters to harvested equilibrium
urch<-(y$Ul/out$Ul[200])
lob<-(y$L/out$L[200])
tot<-(kelp+urch+lob)/3   #take average of relative biomass across 3 harvested groups
tot<-as.data.frame(tot)
maxtot<-max(tot)   #calculate maximum of avg relative biomass
totvol<-(maxtot-tot[200,1])/tot[200,1]   #calculate volatility of avg relative biomass
print(totvol)  ##volatility=0.1603

y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  # 15.537   #with NCE: 13.521
fishedvol<-(maxfished-y[200,7])/y[200,7] 
print(fishedvol)
maxkelp<-max(y$A)
kelpvol<-(maxkelp-y[200,2])/y[200,2]
print(kelpvol)
#RETURN TIME###
#Calculate manually by looking at tot dataframe, add up duration of transient periods until within 10% of equilibria
#6 yrs


#############################################
###SEQUENTIAL followed by Synchronous########
###Need to adjust order depending on ########
#####Simulation##############################
#############################################
###----------parameters----------------------###
state<-c(A=2459.864, Us=8.209526, Um=0.6845794, Ul=0.1081592, L=10.65874) #Unexploited equilibrium
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Fk=0.1        #kelp harvest mortality
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2  ##rate of small urchin consumption by lobsters 
deltaLm=0.15 ##rate of med urchin consumption by lobsters  
deltaLl=0.05 ##rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment   #0.75 for barrens inits, #0.6
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model to set up harvested equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=kelpharvest, parms=parms))  
tail(out,7)

#Now, close the first fishery, use harvested equilibrium population values from above
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
Fk=0.0        #close kelp harvest

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out2=as.data.frame(lsoda(y=newstate, times=times, func=kelpharvest, parms=parms))  
tail(out2,7)

#now stop fishing for both other species
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200], L=out2$L[200])
FU=0.00       #close urchin fishery
FL=0.0        #close lobster fishery

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out3=as.data.frame(lsoda(y=newerstate, times=times, func=kelpharvest, parms=parms))
tail(out3,7)

y<-rbind(out2,out3)
y$time<-seq(1,400)

kelp<-(y$A/out$A[200])   #scale kelp, large urchins & lobsters to harvested equilibrium
urch<-(y$Ul/out$Ul[200])
lob<-(y$L/out$L[200])
tot<-(kelp+urch+lob)/3   #take average of relative biomass across 3 harvested groups
maxtot<-max(tot)   #calculate maximum of avg relative biomass
totvol<-(maxtot-tot[400])/tot[400]   #calculate volatility of avg relative biomass
print(totvol)  ##Basal First: 0.1581    Prey first: 0.172     Pred First: 0.1215

y$FishedDens<-y$Ul+y$L
maxfished<-max(y$FishedDens)  
fishedvol<-(maxfished-y[400,7])/y[400,7] 
print(fishedvol)
maxkelp<-max(y$A)
kelpvol<-(maxkelp-y[400,2])/y[400,2]
print(kelpvol)


#RETURN TIME###
#Calculate manually by looking at tot dataframe, add up duration of transient periods until within 10% of equilibria
plot(tot, type="l")
###Basal First: 0+6=6    Prey first:0+6=6      Pred First: 6+0=6
head(tot,15)
tail(tot,1)
print(tot[200:220])


#############################################
###SEQUENTIAL ###############################
###Need to adjust order depending on ########
#####Simulation##############################
#############################################
###----------parameters----------------------###
state<-c(A=2459.864, Us=8.209526, Um=0.6845794, Ul=0.1081592, L=10.65874) #Unexploited equilibrium
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Fk=0.1        #kelp harvest mortality
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2  ##rate of small urchin consumption by lobsters 
deltaLm=0.15 ##rate of med urchin consumption by lobsters  
deltaLl=0.05 ##rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment   #0.75 for barrens inits, #0.6
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model to set up harvested equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=kelpharvest, parms=parms))  
tail(out,7)

#Now, close the first fishery, use harvested equilibrium population values from above
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.0        #close urchin harvest

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out2=as.data.frame(lsoda(y=newstate, times=times, func=kelpharvest, parms=parms))  
tail(out2,7)

#now stop fishing for next species
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200], L=out2$L[200])
Fk=0.00       #close kelp harvest

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out3=as.data.frame(lsoda(y=newerstate, times=times, func=kelpharvest, parms=parms))
tail(out3,7)

#now close the final fishery
neweststate<-c(A=out3$A[200],Us=out3$Us[200],Um=out3$Um[200],Ul=out3$Ul[200], L=out3$L[200])
FL=0.0        #close lobster fishery

parms=c(r=r, Ka=Ka, Fk=Fk, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out4=as.data.frame(lsoda(y=neweststate, times=times, func=kelpharvest, parms=parms))
tail(out4,7)

y<-rbind(out2,out3,out4)
y$time<-seq(1,600)

kelp<-(y$A/out$A[200])   #scale kelp, large urchins & lobsters to harvested equilibrium
urch<-(y$Ul/out$Ul[200])
lob<-(y$L/out$L[200])
tot<-(kelp+urch+lob)/3   #take average of relative biomass across 3 harvested groups
maxtot<-max(tot)   #calculate maximum of avg relative biomass
totvol<-(maxtot-tot[600])/tot[600]   #calculate volatility of avg relative biomass
print(totvol) #Basal-Pred-Prey:0.1287 Basal-Prey-Pred:0.1706, Prey-Basal-Pred:0.1706 Prey-Pred-Basal:0.1648, Pred-Prey-Basal:0.1215 Pred-Basal-Prey:0.1215 
plot(tot, type="l")

y$FishedDens<-y$Ul+y$L
maxfished<-max(y$FishedDens)  
fishedvol<-(maxfished-y[600,7])/y[600,7] 
print(fishedvol)
maxkelp<-max(y$A)
kelpvol<-(maxkelp-y[600,2])/y[600,2]
print(kelpvol)

#RETURN TIME###
#Calculate manually by looking at tot dataframe, add up duration of transient periods until within 10% of equilibria
###All = 6 *Something probably off with thse because of relativizing
head(tot,15)
tail(tot,1)
print(tot[400:420])
print(tot[400])

