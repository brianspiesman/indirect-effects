rm(list=ls())
#library(MASS);library(lattice);library(rootSolve)

Equalib = matrix(0,1,5)
Steady = matrix(0,1,1)

gmin = 0
gmax = .1
res = .01
giter = 1+((gmax-gmin)/res)
gam1=seq(gmin,gmax, by=res) 				#Benefit provided by pollinator
gam2=seq(gmin,gmax, by=res) 				#Benefit provided by pollinator


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

pb = winProgressBar(title = "progress bar", min = 0, max = giter, width = 300) #Progress bar
for (i in 1:giter){
for (j in 1:giter){
Sys.sleep(0.0001)
setWinProgressBar(pb, i, title=paste( round(i/(giter)*100, 1), "% complete"))

Pars <- c(
a12 = .5,				
a21 = .5,					
g1 = gam1[i],
g2 = gam2[j],
a1 = .09,
a2 = .09,
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

out <- runsteady(func = CommMod, y = State, parms = Pars, times = c(0,Inf), maxsteps = 250000)

Equalib = rbind(Equalib,t(as.matrix(out$y)))
Steady = rbind(Steady,attr(out,"steady"))
}
}

Equalib = Equalib[2:((giter^2)+1),]
	for(i in 1:(giter^2)){
		for (j in 1:5){
			if (Equalib[i,j]<.00001){Equalib[i,j]=0}
		}
	}
Steady = Steady[2:((giter^2)+1),]

N1e = Equalib[,1]
N2e = Equalib[,2]
M1e = Equalib[,3]
M2e = Equalib[,4]
Pe = Equalib[,5]

g1 = as.matrix(rep(seq(gmin,gmax, by=res),each=giter))
g2 = as.matrix(rep(seq(gmin,gmax, by=res),giter))

output = cbind(N1e,N2e,M1e,M2e,Pe,g1,g2)
#write.table(output,file="coex_p1_a8.txt")
close(pb)

x = seq(1:giter^2)
for (i in 1:giter^2){
	x[i] = N1e[i]*N2e[i]
	if (x[i]>0){x[i]=1}
}
levelplot(M1e~g1*g2)








