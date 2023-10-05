# http://codepad.org/lOGkuR1X

import math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import leastsq
import warnings

e = math.e
warnings.filterwarnings('ignore')

file3800 = "3800K_EOS_data.txt"
data3800 = np.loadtxt(file3800)
volumeData3800 = data3800[:,0]
avePressureData3800 = data3800[:,1]
avePressureErrorData3800 = data3800[:,2]

file4000 = "4000K_EOS_data.txt"
data4000 = np.loadtxt(file4000)
volumeData4000 = data4000[:,0]
avePressureData4000 = data4000[:,1]
avePressureErrorData4000 = data4000[:,2]

file4200 = "4200K_EOS_data.txt"
data4200 = np.loadtxt(file4200)
volumeData4200 = data4200[:,0]
avePressureData4200 = data4200[:,1]
avePressureErrorData4200 = data4200[:,2]

Vplot = np.linspace(2.0, 10.0, 800000)
    
V0initial = 6.661
B0initial = 217.0
B0dashinitial = 4.39

def birchmurnaghan(parameters, volumes):
    V0 = parameters[0]
    B0 = parameters[1]
    B0dash = parameters[2]
    Term1 = (3*B0)/2
    Term2 = (V0/volumes)**(7/3)
    Term3 = 1+((0.75*(B0dash-4))*(((V0/volumes)**(2/3))-1))
    P = Term1*Term2*Term3
    return P

def error(parametervals, pressvals, volvals):
    misfits = pressvals - birchmurnaghan(parametervals, volvals)
    #sqrMisfits = misfits**2
    #meanSqrMisfit = np.sum(sqrMisfits)/len(sqrMisfits)
    #rmsMisfit = (meanSqrMisfit)**0.5
    return misfits

plt.xlabel('Simulation Volume, V (Å³/atom)')
plt.ylabel("Pressure, P (GPa)")
plt.title("Fe-2.7 wt% Si Birch-Murnaghan EOS")

initialguesses3800 = [V0initial, B0initial, B0dashinitial]
paramsguess3800, ier = leastsq(error, initialguesses3800, args = (avePressureData3800,volumeData3800))
finalfitPs3800 = birchmurnaghan(paramsguess3800, Vplot)
print("3800K EOS Parameters [V0, B0, B0']: "+str(paramsguess3800[0])+", "+str(paramsguess3800[1])+", "+str(paramsguess3800[2]))
plt.errorbar(volumeData3800, avePressureData3800, yerr=avePressureErrorData3800, color='limegreen', fmt='*')
#plt.plot(Vplot, finalfitPs3800, color='black', linewidth=2)
plt.plot(Vplot, finalfitPs3800, color='limegreen', label="3800K Fe-2.7 wt% Si (Pressure Data Fit)", linewidth=2, zorder=70)

initialguesses4000 = [V0initial, B0initial, B0dashinitial]
paramsguess4000, ier = leastsq(error, initialguesses4000, args = (avePressureData4000,volumeData4000))
finalfitPs4000 = birchmurnaghan(paramsguess4000, Vplot)
print("4000K EOS Parameters [V0, B0, B0']: "+str(paramsguess4000[0])+", "+str(paramsguess4000[1])+", "+str(paramsguess4000[2]))
plt.errorbar(volumeData4000, avePressureData4000, yerr=avePressureErrorData4000, color='firebrick', fmt='*')
#plt.plot(Vplot, finalfitPs4000, color='black', linewidth=2)
plt.plot(Vplot, finalfitPs4000, color='firebrick', label="4000K Fe-2.7 wt% Si (Pressure Data Fit)", linewidth=2, zorder=70)

initialguesses4200 = [V0initial, B0initial, B0dashinitial]
paramsguess4200, ier = leastsq(error, initialguesses4200, args = (avePressureData4200,volumeData4200))
finalfitPs4200 = birchmurnaghan(paramsguess4200, Vplot)
print("4200K EOS Parameters [V0, B0, B0']: "+str(paramsguess4200[0])+", "+str(paramsguess4200[1])+", "+str(paramsguess4200[2]))
plt.errorbar(volumeData4200, avePressureData4200, yerr=avePressureErrorData4200, color='darkturquoise', fmt='*')
#plt.plot(Vplot, finalfitPs4200, color='black', linewidth=2)
plt.plot(Vplot, finalfitPs4200, color='darkturquoise', label="4200K Fe-2.7 wt% Si (Pressure Data Fit)", linewidth=2, zorder=70)

print("-------------------------")

fischerparams = [7.203, 129.1, 5.29]
fischer2params = [6.661, 217.0, 4.39]
hiraoparams = [6.71, 198.0, 4.7]
linparams = [6.882, 141.0, 5.7]

fischerV0 = fischerparams[0]
fischer2V0 = fischer2params[0]
hiraoV0 = hiraoparams[0]
linV0 = linparams[0]

fischerPs = birchmurnaghan(fischerparams, Vplot)
fischer2Ps = birchmurnaghan(fischer2params, Vplot)
hiraoPs = birchmurnaghan(hiraoparams, Vplot)
linPs = birchmurnaghan(linparams, Vplot)

def extrap(initialα, minTExtrap, maxTExtrap, ExtrapParams, ExtrapPs):
    
    lenTExtrap = int(((maxTExtrap - minTExtrap)/10)+1)
    counter = 0
    TExtrapIincrease = 0.0
    
    TExtrap = np.zeros(lenTExtrap)
    VExtrap = np.zeros(lenTExtrap)
    
    for eachvalue in TExtrap:
        TExtrap[counter] = minTExtrap + TExtrapIincrease
        TExtrapIincrease = TExtrapIincrease + 10.0
        counter = counter + 1
     
    initialVExtrap = ExtrapParams[0]
    #deltaVExtrap = minTExtrap*initialVExtrap*α
    #counter2 = 0
    counter2=1
    VExtrap[0]=initialVExtrap

    finalα=8.72712422581875E-06
    αincrease=(initialα-finalα)/lenTExtrap
    counter3=0
    αVals=np.zeros(lenTExtrap)
    α=finalα
    
    for eachValue in αVals:
        αVals[counter3]=α
        α=α+αincrease
        counter3=counter3+1
    
    for eachvalue in TExtrap:
        #VExtrap[counter2] = initialVExtrap+deltaVExtrap
        #deltaVExtrap = (TExtrap[counter2]-minTExtrap)*VExtrap[counter2]*α
        #counter2=counter2+1
        exponent=(TExtrap[counter2]*αVals[counter2])-((TExtrap[counter2]-10.0)*αVals[counter2])
        term4=(e**(exponent))-1
        deltaVExtrap=VExtrap[counter2-1]*term4
        VExtrap[counter2]=VExtrap[counter2-1]+deltaVExtrap
        #α=α*(e**(α*TExtrap[counter2]))
        counter2=counter2+1
        if counter2>(len(TExtrap)-1):
            break
    
    finalVExtrap=(VExtrap[-1])
    
    initialVExtrapIntersect = (np.abs(Vplot-initialVExtrap)).argmin()
    finalVExtrapIntersect = (np.abs(Vplot-finalVExtrap)).argmin()

    VExtrapTransformx = finalVExtrap - initialVExtrap
    VExtrapTransformy =  finalfitPs3800[finalVExtrapIntersect] - ExtrapPs[initialVExtrapIntersect]
    
    VExtraplinex = [initialVExtrap, finalVExtrap]
    VExtrapliney = [ExtrapPs[initialVExtrapIntersect], finalfitPs3800[finalVExtrapIntersect]]
    
    VplotTransformed = Vplot+VExtrapTransformx
    PsTransformed = ExtrapPs+VExtrapTransformy
    
    return(VplotTransformed, PsTransformed, VExtraplinex, VExtrapliney, finalVExtrap)
    
fischerPsExtrap = extrap(4.8987832e-5, 300.0, 3800.0, fischerparams, fischerPs)
fischer2PsExtrap = extrap(4.7880112e-5, 0.0, 3800.0, fischer2params, fischer2Ps)
hiraoPsExtrap = extrap(4.8261660e-5, 300.0, 3800.0, hiraoparams, hiraoPs)
linPsExtrap = extrap(4.949246e-5, 300.0, 3800.0, linparams, linPs)

fischerparams = [fischerPsExtrap[4], 129.1, 5.29]
fischer2params = [fischer2PsExtrap[4], 217.0, 4.39]
hiraoparams = [hiraoPsExtrap[4], 198.0, 4.7]
linparams = [linPsExtrap[4], 141.0, 5.7]

print("3800K [300K] Fe-9 wt% Si (Fischer et al., 2014) EOS Parameters [V0, V0exp, B0, B0']: "+str(fischerV0)+", "+str(fischerparams[0])+", "+str(fischerparams[1])+", "+str(fischerparams[2]))
print("3800K [0K] Fe-3 wt% Si (Fischer et al., 2014) EOS Parameters [V0, V0exp, B0, B0']: "+str(fischer2V0)+", "+str(fischer2params[0])+", "+str(fischer2params[1])+", "+str(fischer2params[2]))
print("3800K [300K] Fe-8.7 wt% Si (Hirao et al., 2004) EOS Parameters [V0, V0exp, B0, B0']: "+str(hiraoV0)+", "+str(hiraoparams[0])+", "+str(hiraoparams[1])+", "+str(hiraoparams[2]))
print("3800K [300K] Fe-8 wt% Si (Lin et al., 2003) EOS Parameters [V0, V0exp, B0, B0']: "+str(linV0)+", "+str(linparams[0])+", "+str(linparams[1])+", "+str(linparams[2]))

print("-------------------------")

fischerPs = birchmurnaghan(fischerparams, Vplot)
fischer2Ps = birchmurnaghan(fischer2params, Vplot)
hiraoPs = birchmurnaghan(hiraoparams, Vplot)
linPs = birchmurnaghan(linparams, Vplot)

fischerPsIntersect=(np.abs(fischerPs-323.0)).argmin()
fischer2PsIntersect=(np.abs(fischer2Ps-323.0)).argmin()
hiraoPsIntersect=(np.abs(hiraoPs-323.0)).argmin()
linPsIntersect=(np.abs(linPs-323.0)).argmin()

Vaverage=(Vplot[fischerPsIntersect] + Vplot[fischer2PsIntersect] + Vplot[hiraoPsIntersect] + Vplot[linPsIntersect])/4
VaverageIntersect=(np.abs(Vplot-Vaverage)).argmin()

print("3800K [300K] Fe-9 wt% Si (Fischer et al., 2014) Volume at 323 GPa: "+str(Vplot[fischerPsIntersect]))
print("3800K [0K] Fe-3 wt% Si (Fischer et al., 2014) Volume at 323 GPa: "+str(Vplot[fischer2PsIntersect]))
print("3800K [300K] Fe-8.7 wt% Si (Hirao et al., 2004) Volume at 323 GPa: "+str(Vplot[hiraoPsIntersect]))
print("3800K [300K] Fe-8 wt% Si (Lin et al., 2003) Volume at 323 GPa: "+str(Vplot[linPsIntersect]))
print("Average of Literature Volumes at 323 GPa: "+str(Vplot[VaverageIntersect]))

plt.scatter(Vplot[fischerPsIntersect], fischerPs[fischerPsIntersect], color='lightcoral', label="3800K [300K] Fe-9 wt% Si (Fischer et al., 2014)", s=50, edgecolors="black", zorder=100)
plt.scatter(Vplot[fischer2PsIntersect], fischer2Ps[fischer2PsIntersect], color='darkorange', label="3800K [0K] Fe-3 wt% Si (Fischer et al., 2014)", s=50, edgecolors="black", zorder=100)
plt.scatter(Vplot[hiraoPsIntersect], hiraoPs[hiraoPsIntersect], color='royalblue', label="3800K [300K] Fe-8.7 wt% Si (Hirao et al., 2004)", s=50, edgecolors="black", zorder=100)
plt.scatter(Vplot[linPsIntersect], linPs[linPsIntersect], color='gold', label="3800K [300K] Fe-8 wt% Si (Lin et al., 2003)", s=50, edgecolors="black", zorder=100)
plt.scatter(Vplot[VaverageIntersect], 323.0, color='none', label="Average of Literature Volumes at 323 GPa", marker="o", s=150, edgecolors="black", zorder=100)

#plt.plot(fischerPsExtrap[2], fischerPsExtrap[3], color='lightcoral', linestyle="--", alpha=0.7)
#plt.plot(fischer2PsExtrap[2], fischer2PsExtrap[3], color='darkorange', linestyle="--", alpha=0.7)
#plt.plot(hiraoPsExtrap[2], hiraoPsExtrap[3], color='royalblue', linestyle="--", alpha=0.7)
#plt.plot(linPsExtrap[2], linPsExtrap[3], color='gold', linestyle="--", alpha=0.7)

#plt.plot(fischerPsExtrap[0], fischerPsExtrap[1], color='lightcoral', label="300K Fe-9 wt% Si (Fischer et al., 2014)")
#plt.plot(fischer2PsExtrap[0], fischer2PsExtrap[1], color='darkorange', label="0K Fe-3 wt% Si (Fischer et al., 2014)")
#plt.plot(hiraoPsExtrap[0], hiraoPsExtrap[1], color='royalblue', label="300K Fe-8.7 wt% Si (Hirao et al., 2004)")
#plt.plot(linPsExtrap[0], linPsExtrap[1], color='gold', label="300K Fe-8 wt% Si (Lin et al., 2003)")

gpalinex = (3.0, 10.0)
gpaliney = (323.0, 323.0)

plt.plot(gpalinex, gpaliney, linestyle = "-.", alpha=0.7, color='darkviolet', zorder=1)

plt.text(7.475, 325, '323 GPa', fontsize=9)

print("-------------------------")

intersect323_3800 = (np.abs(finalfitPs3800-323.0)).argmin()
intersectlinex323_3800 = (Vplot[intersect323_3800], Vplot[intersect323_3800])
intersectliney323_3800 = (0.0, 323.0)
idealVolume3800=Vplot[intersect323_3800]
idealVolumeClean3800 = str(round(idealVolume3800,3))
print("Ideal Volume to reach 323 GPa at 3800K: "+str(idealVolume3800)+" Å³/atom")
plt.plot(intersectlinex323_3800, intersectliney323_3800, linestyle = ":", alpha=0.7, color='darkviolet', zorder=1)

intersect323_4000 = (np.abs(finalfitPs4000-323.0)).argmin()
intersectlinex323_4000 = (Vplot[intersect323_4000], Vplot[intersect323_4000])
intersectliney323_4000 = (0.0, 323.0)
idealVolume4000=Vplot[intersect323_4000]
idealVolumeClean4000 = str(round(idealVolume4000,3))
print("Ideal Volume to reach 323 GPa at 4000K: "+str(idealVolume4000)+" Å³/atom")
plt.plot(intersectlinex323_4000, intersectliney323_4000, linestyle = ":", alpha=0.7, color='darkviolet', zorder=1)

intersect323_4200 = (np.abs(finalfitPs4200-323.0)).argmin()
intersectlinex323_4200 = (Vplot[intersect323_4200], Vplot[intersect323_4200])
intersectliney323_4200 = (0.0, 323.0)
idealVolume4200=Vplot[intersect323_4200]
idealVolumeClean4200 = str(round(idealVolume4200,3))
print("Ideal Volume to reach 323 GPa at 4200K: "+str(idealVolume4200)+" Å³/atom")
plt.plot(intersectlinex323_4200, intersectliney323_4200, linestyle = ":", alpha=0.7, color='darkviolet', zorder=1)

plt.xlim(6.6, 7.6)
plt.ylim(200, 340)

plt.legend(loc=3, prop={'size': 8})

plt.tick_params(axis="y", direction='inout')
plt.tick_params(axis="x", direction='inout')

filename="equations_of_state.png"
#plt.savefig(filename,dpi=300)

plt.show()