import numpy as np
import matplotlib.pyplot as plt
import math

pi=math.pi
e=math.e

file = "gruneisen_originals.txt"
data = np.loadtxt(file)
radiusData = data[:,0]
gruneisenData = data[:,1]

rVals=np.zeros(348001)
gruneisenVals=np.zeros(348001)
counter=0
    
for eachVal in rVals:
        
    rVals[counter]=counter/100
    gruneisenIndex=(np.abs(radiusData-rVals[counter])).argmin()
    gruneisenVals[counter]=gruneisenData[gruneisenIndex]
    counter=counter+1

def supercoolingCalc(T0Fe, χ0, Cp, α0, ρ0, γ, G,):
    
    counter2=0
    
    ΔT0=np.zeros(348001)
    
    for eachVal in rVals:
        
        r=rVals[counter2]*1000
        #γ=gruneisenVals[counter2]
        
        Lp = ((3*Cp)/(2*pi*α0*ρ0*G))**0.5
        
        term1 = 1-((e**((-2)*(1-(1/(3*γ)))*((r**2)/(Lp**2))))/(e**(-((r**2)/(Lp**2)))))
        
        term2 = ((1)/(e**(-((r**2)/(Lp**2)))))-1
        
        ΔT0[counter2] = (T0Fe*term1)+(χ0*term2)
        
        counter2=counter2+1
    
    return ΔT0

Ric=1221
RicStr=str(Ric)+" km"

#fig,ax=plt.subplots()

plt.plot([0, 2500],[Ric, Ric], color="royalblue", linestyle="--", label="Inner Core Radius")
#ax.text(230, Ric+50, RicStr, size=8)
plt.text(390, Ric+50, RicStr, size=8)

ΔT0min=supercoolingCalc(6000, 500, 774, 1.02e-5, 12000, 1.3, 6.6743e-11)
ΔT0max=supercoolingCalc(7000, 1000, 946, 1.95e-5, 13000, 1.7, 6.6743e-11)

intersect = (np.abs(rVals-Ric)).argmin()
ΔTmin=ΔT0min[intersect]
ΔTmax=ΔT0max[intersect]

print("Minimum Initial Supercooling of the Inner Core: "+str(ΔTmin)+" K")
print("Maximum Initial Supercooling of the Inner Core: "+str(ΔTmax)+" K")

#ax.plot([ΔTmin, ΔTmin], [0, Ric], color='black', linewidth=0.25)
#ax.plot([ΔTmax, ΔTmax], [0, Ric], color='black', linewidth=0.25)

#ax.fill_betweenx([0, Ric], [ΔTmin, ΔTmin], [ΔTmax, ΔTmax], color="royalblue", label="Estimated Initial Supercooling of the Inner Core")

plt.plot(ΔT0min, rVals, color='black', linewidth=0.5)
plt.plot(ΔT0max, rVals, color='black', linewidth=0.5)

plt.fill_betweenx(rVals, ΔT0min, ΔT0max, color='lightcoral', label="Modelled Saturation Surface of a Frozen Iron Mixture")
ΔTminStr=str(int(round(ΔTmin, 0)))
ΔTmaxStr=str(int(round(ΔTmax, 0)))
#ax.text(180, 50, ΔTminStr+" - "+ΔTmaxStr+" K", size=8)
plt.text(330, 50, ΔTminStr+" - "+ΔTmaxStr+" K", size=8)

plt.plot([ΔTmin, ΔTmin], [0, Ric], color='royalblue', linestyle="-.", label="Centre Undercooling for the Inner Core Radius")
plt.plot([ΔTmax, ΔTmax], [0, Ric], color='royalblue', linestyle="-.")

plt.plot(ΔTmin, 0, marker='o', ms=7, mec="black", color='royalblue')
plt.plot(ΔTmax, 0, marker='o', ms=7, mec="black", color='royalblue')

plt.title("Saturation Surface for a given Earth Centre Undercooling")
plt.xlabel("Undercooling at Earth's centre, ΔT₀ (K)")
plt.ylabel("Radius, r (km)")

#ax.set_xlim(0, 1500)
plt.xlim(0, 2500)
plt.ylim(0, 3480)

#ax2=ax.twiny()
#ax2.plot(gruneisenData, radiusData, color="darkorange", label="Gruneisen Parameter Value", linestyle=":")
##ax2.set_xlabel("Gruneisen Parameter, γ")
#ax2.set_xlim(min(gruneisenData), max(gruneisenData))

plt.legend(fontsize=8, loc=2)
#ax2.legend(fontsize=8, loc=4)

filename="initial_supercooling_estimate.png"
#plt.savefig(filename,dpi=300)

plt.show()


    
    
    
    
    
    