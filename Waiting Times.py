import math
import numpy as np
import matplotlib.pyplot as plt
import warnings
import scipy.optimize as sco

def equation6(kB, τvic, undercoolings, multiplier, γ, hf, Tm, I0, δTc, γ2, hf2, Tm2, I02, δTc2, γ3, hf3, Tm3, I03, δTc3, γ4, hf4, Tm4, I04, δTc4, plottype):
    
    if plottype=="final":
        print("Iron-Silicon Melting Temperature: "+str(iron_sil_tm)+" K")
    
    warnings.filterwarnings('ignore')
    pi = math.pi
    
    #WAITINGTIMES=[(1.54917714870026E-36), (2.58177528262359E-35), (4.89212582346998E-36)]
    #UNDERCOOLINGS=[(3800-Tm3), (4000-Tm3), (4200-Tm3)]
    
    WAITINGTIMES=[(1.54917714870026E-36), (2.58177528262359E-35)]
    UNDERCOOLINGS=[(3800-Tm3), (4000-Tm3)]
    
    wtError1=(np.array(WAITINGTIMES))/(50**0.5)
    #wtError2=(1.09682063809614E-36, 2.34646115128324E-36, 2.38724111795338E-36)
    wtError2=(1.09682063809614E-36, 2.34646115128324E-36)
    wtError=wtError1+wtError2
    ucError=75.16671764
    
    WAITINGTIMESMAX = np.zeros(len(WAITINGTIMES))
    WAITINGTIMESMIN = np.zeros(len(WAITINGTIMES))
    
    WAITINGTIMESMAX = WAITINGTIMES + wtError
    WAITINGTIMESMIN = WAITINGTIMES - wtError
    
    def waitingtimes(γx, hfx, Tmx, I0x, δTcx):
    
        arraylen=(undercoolings*multiplier)+1
        τv =  np.zeros(arraylen)
        T = np.zeros(arraylen)
        
        δTx = np.linspace(0, -undercoolings, arraylen)   
    
        for counter1, eachValue in enumerate(T):
            T[counter1] = Tmx + δTx[counter1]
        
        for counter2, eachValue in enumerate(τv):
            hc = 1.0-(7.05e-5*δTx[counter2])
            a = 16.0*pi*(γx**3.0)*(Tmx**2.0)
            b = 3.0*kB*(T[counter2])*((δTx[counter2])**2.0)*(hfx**2.0)*(hc**2.0)
            τv[counter2] = (I0x/2.0)*np.exp(a/b)
            
        τvConvertx=(np.where(np.isinf(τv),1e+260,τv))
            
        intersectx=(np.abs(τvConvertx-τvic)).argmin()
        
        if δTcx == δTc3:
            δTcx=δTx[intersectx]
        
        critUndercoolingPos=(np.argwhere(δTx==δTcx))
        critWaitingTime=str(τv[critUndercoolingPos]).replace('[','').replace(']','').replace('\'','').replace('\"','')
        if plottype=="final":    
            print("Waiting Time For Critical Undercooling of "+str(δTcx)+" K: "+critWaitingTime+" sm³")
        
        return δTx, τvConvertx, intersectx
    
    def waitingtimesfit(params, undercoolingVals):
        
        γy=params[0]
        I0y=params[1]
        
        arrayleny=len(undercoolingVals)
        τvy =  np.zeros(arrayleny)
        Ty = np.zeros(arrayleny)
    
        for counter3, eachValue in enumerate(Ty):
            Ty[counter3] = Tm3 + undercoolingVals[counter3]
        
        for counter4, eachValue in enumerate(τvy):
            hcy = 1.0-(7.05e-5*undercoolingVals[counter4])
            ay = 16.0*pi*(γy**3.0)*(Tm3**2.0)
            by = 3.0*kB*(Ty[counter4])*((undercoolingVals[counter4])**2.0)*(hf3**2.0)*(hcy**2.0)
            τvy[counter4] = (I0y/2)*np.exp(ay/by)
        
        τvConverty=(np.where(np.isinf(τvy),1e+260,τvy))
        
        return τvConverty
        
    pureiron=waitingtimes(γ, hf, Tm, I0, δTc)
    ironoxygen=waitingtimes(γ2, hf2, Tm2, I02, δTc2)
    idealcurve=waitingtimes(γ4, hf4, Tm4, I04, δTc4)
    
    def error(freeparams, Waitingtimes, Undercoolings):
        if freeparams[0]>0 and freeparams[1]>0:
            misfits = Waitingtimes - waitingtimesfit(freeparams, Undercoolings)
            return misfits
        else:
            return [1e6, 1e6]
            #return [1e6, 1e6, 1e6]
    
    initialguesses=[γ3, I03]
    
    paramsguess, ier = sco.leastsq(error, initialguesses, args = (WAITINGTIMES, UNDERCOOLINGS))
    #paramsguess, ier = sco.leastsq(error, initialguesses, args = (WAITINGTIMESMAX, UNDERCOOLINGS))
    #paramsguess, ier = sco.leastsq(error, initialguesses, args = (WAITINGTIMESMIN, UNDERCOOLINGS))
    
    γ3 = paramsguess[0]
    I03 = paramsguess[1]
    
    ironsilicon=waitingtimes(γ3, hf3, Tm3, I03, δTc3)
    
    if plottype=="final":
        print("γ: "+str(γ3)+", I0: "+str(I03))
        
    if plottype=="final":
        print("Waiting Times (Calculated): "+str(WAITINGTIMES)[1:-1])
        #print("Waiting Times (Maximum): "+str(WAITINGTIMESMAX)[1:-1])
        #print("Waiting Times (Minimum): "+str(WAITINGTIMESMIN)[1:-1])
        
        print("Undercoolings: "+str(UNDERCOOLINGS)[1:-1])
    
    δT=pureiron[0]
    τvConvert=pureiron[1]
    intersect=pureiron[2]
    δT2=ironoxygen[0]
    τv2Convert=ironoxygen[1]
    intersect2=ironoxygen[2]
    δT3=ironsilicon[0]
    τv3Convert=ironsilicon[1]
    intersect3=ironsilicon[2]
    δT4=idealcurve[0]
    τv4Convert=idealcurve[1]
    intersect4=idealcurve[2]
    
    if plottype=="final" or plottype=="initial":
        plt.plot(pureiron[0], pureiron[1], label='Pure Iron System', color='limegreen')
        plt.plot(δT2, τv2Convert, label='Iron-Oxygen System', color='darkorange')
    if plottype=="final" or plottype=="zoomed":
        plt.plot(δT3, τv3Convert, label='Iron-Silicon System', color='firebrick')
    if plottype=="initial":
        plt.plot(δT4, τv4Convert, label='Modelled Secular Cooling', color='darkviolet')
    
    estimatelinex = [0, δT[-1]]
    estimateliney = [τvic, τvic]
    if plottype=="final" or plottype=="initial":
        plt.plot(estimatelinex, estimateliney, '-.', label="τv Estimate from Earth's Secular Cooling", color='darkturquoise')
    
    intersectlinex=[δT[intersect], δT[intersect]]
    intersectliney=[τvConvert[-1], τvic]
    if plottype=="final" or plottype=="initial":
        plt.plot(intersectlinex, intersectliney, '--', color='limegreen', alpha=0.6)
        plt.plot(δT[intersect],10e-40,marker='o', ms=7, mec="black", color='limegreen', label='Pure Iron Undercooling')
    
    intersect2linex=[δT2[intersect2], δT2[intersect2]]
    intersect2liney=[τv2Convert[-1], τvic]
    if plottype=="final" or plottype=="initial":
        plt.plot(intersect2linex, intersect2liney, '--', color='darkorange', alpha=0.6)
        plt.plot(δT2[intersect2],10e-40,marker='o', ms=7, mec="black", color='darkorange', label='Iron-Oxygen Undercooling')
    
    intersect3linex=[δT3[intersect3], δT3[intersect3]]
    intersect3liney=[τv3Convert[-1], τvic]
    if plottype=="final" or plottype=="zoomed":
        plt.plot(intersect3linex, intersect3liney, '--', color='firebrick', alpha=0.6)
        plt.plot(δT3[intersect3],10e-40,marker='o', ms=7, mec="black", color='firebrick', label='Iron-Silicon Undercooling')
    
    intersect4linex=[δT4[intersect4], δT4[intersect4]]
    intersect4liney=[τv4Convert[-1], τvic]
    if plottype=="initial":
        plt.plot(intersect4linex, intersect4liney, '--', color='darkviolet', alpha=0.6)
        plt.plot(δT4[intersect4],10e-40,marker='o', ms=7, mec="black", color='darkviolet', label='Secular Undercooling')
    
    if plottype=="initial":
        plt.text(-1600, τvic+1e37, '2.3 x 10³⁵ sm³', fontsize=10.5)
    if plottype=="final":
        #plt.text(-75, τvic+1e37, '2.3 x 10³⁵ sm³', fontsize=10.5)
        plt.text(-25, τvic+1e37, '2.3 x 10³⁵ sm³', fontsize=10.5)
    if plottype=="initial":
        plt.text(δT[intersect]-20, 10e-37, str(int(δT[intersect]))+' K', fontsize=8.5, rotation=90)
        plt.text(δT2[intersect2]-20, 10e-37, str(int(δT2[intersect2]))+' K', fontsize=8.5, rotation=90)
    if plottype=="final":
        plt.text(δT3[intersect3]-20, 10e-37, str(int(δT3[intersect3]))+' K', fontsize=8.5, rotation=90)
    if plottype=="initial":
        plt.text(δT4[intersect4]-20, 10e-37, str(int(δT4[intersect4]))+' K', fontsize=8.5, rotation=90)
    
    if plottype=="final":
        plt.errorbar(UNDERCOOLINGS, WAITINGTIMES, xerr=ucError, yerr=wtError, fmt="*", color="firebrick", elinewidth = 2, capsize=3, ecolor="royalblue")
        
    if plottype=="zoomed":
        plt.errorbar(UNDERCOOLINGS, WAITINGTIMES, xerr=ucError, yerr=wtError, fmt="*", color="firebrick", elinewidth = 2, capsize=3, ecolor="royalblue")
        
    plt.xlabel('Undercooling, ΔT (K)')
    plt.ylabel("Waiting Time to Observe a Freezing Event, τv (sm³)")
    plt.title("Plot of Waiting Time for a Given Undercooling (CNT Modelled)")
    if plottype=="initial":
        plt.legend(loc=4, prop={'size': 8})
    if plottype=="zoomed" or plottype=="final":
        plt.legend(loc=1, prop={'size': 8})
    plt.yscale('log')
    if plottype=="initial":
        plt.xlim(-2100, 0)
        plt.ylim(10e-40, 10e60)
    if plottype=="final":
        #plt.xlim(-2450, -550)
        plt.xlim(-2450, 0)
        plt.ylim(10e-40, 10e60)
    if plottype=="zoomed":
        plt.xlim(-2450, -1800)
        plt.ylim(10e-38, 10e-35)
    plt.gca().invert_xaxis()
    plt.tick_params(axis="y", direction='inout')
    plt.tick_params(axis="x", direction='inout')
    filename="waitingtimes_"+plottype+".png"
    #plt.savefig(filename,dpi=300)
    plt.show()
 
#iron_sil_tm = 6139.07535553353
#iron_sil_tm = 6063.90863789632
iron_sil_tm = 6214.24207317073

equation6(1.3806e-23, 2.3978e35, 3850, 100, 1.08, 0.98e10, 6215.0, 0.71e-48, -730.0, 1.02, 0.98e10, 5600.0, 0.79e-45, -675.0, 1.02, 0.98e10, iron_sil_tm , 0.79e-45, -675.0, 0.592, 0.98e10, 6215.0, 0.71e-48, -291.0, "initial")
equation6(1.3806e-23, 2.3978e35, 3850, 100, 1.08, 0.98e10, 6215.0, 0.71e-48, -730.0, 1.02, 0.98e10, 5600.0, 0.79e-45, -675.0, 1.02, 0.98e10, iron_sil_tm , 0.79e-45, -675.0, 0.592, 0.98e10, 6215.0, 0.71e-48, -291.0, "zoomed")
equation6(1.3806e-23, 2.3978e35, 3850, 100, 1.08, 0.98e10, 6215.0, 0.71e-48, -730.0, 1.02, 0.98e10, 5600.0, 0.79e-45, -675.0, 1.02, 0.98e10, iron_sil_tm , 0.79e-45, -675.0, 0.592, 0.98e10, 6215.0, 0.71e-48, -291.0, "final")
    

