import math
import numpy as np
import matplotlib.pyplot as plt

def equation1(minr, maxr, multiplier, γ, hf, Tm, δTc, γ2, hf2, Tm2, δTc2, γ3, hf3, Tm3, δTc3):
    
    pi = math.pi
    hc = 1.0-(7.05e-5*δTc)
    hc2 = 1.0-(7.05e-5*δTc2)
    hc3 = 1.0-(7.05e-5*δTc3)
    gsl = hf*(δTc/Tm)*hc
    gsl2 = hf2*(δTc2/Tm2)*hc2
    gsl3 = hf3*(δTc3/Tm3)*hc3
    arraylen = ((maxr-minr)*multiplier)+1
    deltaG = np.zeros(arraylen)
    deltaGa = np.zeros(arraylen)
    deltaGb = np.zeros(arraylen)
    deltaG2 = np.zeros(arraylen)
    deltaG3 = np.zeros(arraylen)
    rc = ((-2.0*γ)/gsl)/1.0e-10
    rc2 = ((-2.0*γ2)/gsl2)/1.0e-10
    rc3 = ((-2.0*γ3)/gsl3)/1.0e-10
    
    r = np.linspace(minr,maxr,arraylen)
    r = r*1.0e-10
    r2 = np.linspace(minr,maxr,arraylen)
    r2 = r2*1.0e-10
    r3 = np.linspace(minr,maxr,arraylen)
    r3 = r3*1.0e-10
    
    for counter, eachValue in enumerate(deltaG):
        a = ((4.0/3.0)*pi*(r[counter]**3.0)*gsl)
        a2 = ((4.0/3.0)*pi*(r2[counter]**3.0)*gsl2)
        a3 = ((4.0/3.0)*pi*(r3[counter]**3.0)*gsl3)
        b = (4.0*pi*(r[counter]**2.0)*γ)
        b2 = (4.0*pi*(r2[counter]**2.0)*γ2)
        b3 = (4.0*pi*(r3[counter]**2.0)*γ3)
        deltaG[counter] = (a + b)*6.242e18
        deltaGa[counter] = a*6.242e18
        deltaGb[counter] = b*6.242e18
        deltaG2[counter] = (a2 + b2)*6.242e18
        deltaG3[counter] = (a3 + b3)*6.242e18
        
    r = r/1.0e-10
    r2 = r2/1.0e-10
    r3 = r3/1.0e-10

    plt.plot(r, deltaG, label='Pure Iron System', color='limegreen')
    plt.plot(r, deltaGa, label='Free Energy Component', color='limegreen',linestyle='-.', alpha=0.7)
    plt.plot(r, deltaGb, label='Interfacial Energy Component', color='limegreen', linestyle='--', alpha=0.7)
    plt.plot(r2, deltaG2, label='Iron-Oxygen System', color='darkorange')
    #plt.plot(r3, deltaG3, label='Iron-Silicon System', color='gold')
    plt.plot(rc, np.amax(deltaG), marker='o', ms=7, mec="black", label='Critical Radius for Pure Iron of '+str(round(rc,3))+'Å (5sf)', color='limegreen')
    plt.plot(rc2, np.amax(deltaG2), marker='o', ms=7, mec="black", label='Critical Radius for Iron-Oxygen of '+str(round(rc2,3))+'Å (5sf)', color='darkorange')
    #plt.plot(rc3, np.amax(deltaG3), marker='o', ms=7, mec="black", label='Critical Radius for Iron-Silicon of '+str(round(rc3,3))+'Å (5sf)', color='darkviolet')

    print('Critical Radius for Pure Iron of '+str(rc)+'Å')
    print('Critical Radius for Iron-Oxygen of '+str(rc2)+'Å')

    plt.xlabel('Nucleus Radius, r (Å)')
    plt.ylabel("Gibb's Free Energy Change, ΔG (eV)")
    plt.title("Plot of Gibb's Free Energy Change against Nucleus Radius")
    plt.xlim(0, 35)
    plt.ylim(-250, 250)
    plt.tick_params(axis="y", direction='inout')
    plt.tick_params(axis="x", direction='inout')
    plt.legend(loc=3, prop={'size': 8}, )
    
    filename="gibbs_free_energies.png"
    #plt.savefig(filename,dpi=300)
    
    plt.show()
    
equation1(0, 50, 1000, 1.08, 0.98e10, 6215.0, -730.0, 1.02, 0.98e10, 5600.0, -675.0, 1.02, 0.98e10, 5897.75, -675.0)