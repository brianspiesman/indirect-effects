rm(list=ls())
#library(MASS);library(deSolve);library(rootSolve)


Parameters = matrix(0,1,6)
Equalib = matrix(0,1,5)
Steady = matrix(0,1,1)
g = seq(0,.08, by=.001)
reps = length(g)

CommMod <- function (Time, State, Pars) {
	with(as.list(c(State, Pars)), {
		dN1 = r1*N1*(1-(N1+(a21*N2))/(k1+((g1*M1*s1*N1)/(1+s1*h1*N1))))
		dN2 = r2*N2*(1-(N2+(a12*N1))/(k2+((g2*M2*s2*N2)/(1+s2*h2*N2))))
		dM1 = ((b1*M1*s1*N1)/(1+s1*h1*N1))-(d1*M1*(1 + e1*M1))-(a1*P*M1/(1+a1*f*M1))
		dM2 = ((b2*M2*s2*N2)/(1+s2*h2*N2))-(d2*M2*(1 + e2*M2))-(a2*P*M2/(1+a2*f*M2))
		dP = c*P*((a1*M1 + a2*M2)/(1+a1*f*M1+a2*f*M2)) - m*P

		return(list(c(dN1, dN2, dM1, dM2, dP)))
	})
}

pb = winProgressBar(title = "progress bar", min = 0, max = reps, width = 300) #Progress bar

for (i in 1:reps){
Sys.sleep(0.0001)
setWinProgressBar(pb, i, title=paste( round(i/reps*100, 1), "% complete"))

Pars <- c(
a12 = .5, 		#round(runif(1,min=0,max=1.2),4),				
a21 = .5,		#round(runif(1,min=0,max=1.2),4),					
g1 = g[i]+.01, 	#(round(runif(1,min=0,max=0.2),4)),
g2 = g[i], 		#(round(runif(1,min=0,max=0.2),4)),
a1 = .09, 		#round(runif(1,min=0,max=.15),4),
a2 = .095, 		#round(runif(1,min=0,max=.15),4),
r1=1, r2=1,
k1=10, k2=10,
b1=.2, b2=.2,
d1=.1, d2=.1,
e1=.1, e2=.1,
s1=.3, s2=.3,
h1=.2, h2=.2,
c=.1,
m=.1,
f=.1)

State <- c(N1=10, N2=10, M1=5, M2=5, P=5)

out <- runsteady(func = CommMod, y = State, parms = Pars, times = c(0,Inf), maxsteps = 100000)

Parameters = rbind(Parameters,Pars[1:6])
Equalib = rbind(Equalib,t(as.matrix(out$y)))
Steady = rbind(Steady,attr(out,"steady"))
}

Parameters = Parameters[2:(reps+1),]
Equalib = Equalib[2:(reps+1),]
	for(i in 1:reps){
		for (j in 1:5){
			if (Equalib[i,j]<.0001){Equalib[i,j]=0}
		}
	}
Steady = Steady[2:(reps+1),]

#matplot(out[,-1], type = "l", xlab = "time", ylab = "population")
#legend("topright", c("N1", "N2", "M1", "M2", "P"), lty = c(1,2,3,4,5), col = c(1,2,3,4,5), box.lwd = 0)

output = as.data.frame(cbind(Equalib,Parameters));rownames(output) <- seq(1:reps)
N1e = output[,1]
N2e = output[,2]
M1e = output[,3]
M2e = output[,4]
Pe = output[,5]
a12e = output[,6]
a21e = output[,7]
g1e = output[,8]
g2e = output[,9]
a1e = output[,10]
a2e = output[,11]

apcom12 = seq(1:reps)			#empty vector to fill with interaction coefficient
apcom21 = seq(1:reps)
netpred1 = seq(1:reps)
netpred2 = seq(1:reps)
partial.pred12 = seq(1:reps)
partial.comp12 = seq(1:reps)
partial.pred21 = seq(1:reps)
partial.comp21 = seq(1:reps)

r1=1; r2=1							#N intrinsic rate of growth
k1=10; k2=10 						#N carrying capacity
b1=.2; b2=.2   						#M intrinsic rate of growth
d1=.1; d2=.1    						#M density independent mortality
e1=.1; e2=.1    						#M density dependent mortality
s1=.3; s2=.3    						#M encounter rate
h1=.2; h2=.2    						#M processing time
c=.1								#conversion efficiency
m=.1								#P mortality
f=.1								#P processing time

for (i in 1:reps){
State2 <- c(N1=N1e[i], N2=N2e[i], M1=M1e[i], M2=M2e[i], P=Pe[i])
a12=a12e[i]; a21=a21e[i]; g1=g1e[i]; g2=g2e[i]; a1=a1e[i]; a2=a2e[i]
commat = jacobian.full(y = State2, func = CommMod)
icommat = -ginv(t(commat))

apcom12[i] = icommat[3,4]								#net effect of M1 on M2 for each replicate
apcom21[i] = icommat[4,3]								#net effect of M2 on M1 for each replicate
netpred1[i] = icommat[5,3]
netpred2[i] = icommat[5,4]

partial.pred12[i] = (icommat[5,4]*icommat[3,5])/icommat[5,5]
partial.pred21[i] = (icommat[5,3]*icommat[4,5])/icommat[5,5]
partial.comp12[i] = apcom12[i]-partial.pred12[i]
partial.comp21[i] = apcom21[i]-partial.pred21[i]

}
output[output==0] <- NA
output1 = cbind(apcom21,apcom12,partial.pred12,partial.comp12,partial.pred21,partial.comp21,netpred1,netpred2,output,Steady)
output2 = output1[complete.cases(output),]
#write.table(output2,file="testode.txt")
#write.table(output1,file="testode full.txt")
close(pb)
output1
plot(g,partial.pred12)

