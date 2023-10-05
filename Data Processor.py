import numpy as np
import matplotlib.pyplot as plt

def simulationdata(file, plottype):
    run = file[0]+file[1]+file[2]
    
    data = np.loadtxt(file, max_rows=1015000)
    
    openfile = open(file, "r")
    siminfo=openfile.readlines()
    
    startTime=str(siminfo[1015000]).strip()
    endTime=str(siminfo[1015001]).strip()
    fullInitialTemp=str(siminfo[1015002]).strip()
    volume=str(siminfo[1015003]).strip()
    ranSeed=str(siminfo[1015004]).strip()
    simLength=str(siminfo[1015005]).strip()
    
    startTimeClean=startTime[8]+startTime[9]+"/"+startTime[5]+startTime[6]+"/"+startTime[2]+startTime[3]+" at "+startTime[11]+startTime[12]+startTime[13]+startTime[14]+startTime[15]+startTime[16]+startTime[17]+startTime[18]
    endTimeClean=endTime[8]+endTime[9]+"/"+endTime[5]+endTime[6]+"/"+endTime[2]+endTime[3]+" at "+endTime[11]+endTime[12]+endTime[13]+endTime[14]+endTime[15]+endTime[16]+endTime[17]+endTime[18]
    initialTemp=fullInitialTemp[0]+fullInitialTemp[1]+fullInitialTemp[2]+fullInitialTemp[3]
    
    initialTempNum=int(float(siminfo[1015002]))
    
    fullTimeVals = data[:,0]
    fullTempVals = data[:,1]
    fullPotVals = data[:,2]
    fullKinVals = data[:,3]
    fullEnergyVals = data[:,4]
    fullPressVals = data[:,5]
    
    value=0
    value1=0
    value2=0
    value3=0
    value4=0
    value5=0
    value6=0
    value7=0
    
    arrayLen = len(fullTimeVals) - 14999
    
    mainTimeVals = np.empty(arrayLen)
    mainTempVals = np.empty(arrayLen)
    mainPotVals = np.empty(arrayLen)
    mainKinVals = np.empty(arrayLen)
    mainEnergyVals = np.empty(arrayLen)
    mainPressVals = np.empty(arrayLen)
    
    plotLen=(len(mainTimeVals)*1.1)*0.000001
    
    eosTimeVals = np.empty(1000)
    eosTempVals = np.empty(1000)
    eosPotVals = np.empty(1000)
    eosKinVals = np.empty(1000)
    eosEnergyVals = np.empty(1000)
    eosPressVals = np.empty(1000)
    
    equilTimeVals = np.empty(10000)
    equilTempVals = np.empty(10000)
    equilPotVals = np.empty(10000)
    equilKinVals = np.empty(10000)
    equilEnergyVals = np.empty(10000)
    equilPressVals = np.empty(10000)
    
    eightyArrayLen = round((arrayLen*0.8),0)
    twentyPercent = 14.999+((arrayLen-eightyArrayLen)*0.001)
    
    averageTimeVals = np.empty(int(eightyArrayLen))
    averageTempVals = np.empty(int(eightyArrayLen))
    averagePotVals = np.empty(int(eightyArrayLen))
    averageKinVals = np.empty(int(eightyArrayLen))
    averageEnergyVals = np.empty(int(eightyArrayLen))
    averagePressVals = np.empty(int(eightyArrayLen))
    
    for eachValue in fullTimeVals:
        if fullTimeVals[value]>14.999:
            mainTimeVals[value1]=fullTimeVals[value]
            mainTempVals[value1]=fullTempVals[value]
            mainPotVals[value1]=fullPotVals[value]
            mainKinVals[value1]=fullKinVals[value]
            mainEnergyVals[value1]=fullEnergyVals[value]
            mainPressVals[value1]=fullPressVals[value]
            value1=value1+1
        value=value+1
        
    for eachValue in fullTimeVals:
        if fullTimeVals[value2]>13.999 and fullTimeVals[value2]<15.0:
            eosTimeVals[value3]=fullTimeVals[value2]
            eosTempVals[value3]=fullTempVals[value2]
            eosPotVals[value3]=fullPotVals[value2]
            eosKinVals[value3]=fullKinVals[value2]
            eosEnergyVals[value3]=fullEnergyVals[value2]
            eosPressVals[value3]=fullPressVals[value2]
            value3=value3+1
        value2=value2+1
        
    for eachValue in fullTimeVals:
        if fullTimeVals[value4]>4.999 and fullTimeVals[value4]<15.0:
            equilTimeVals[value5]=fullTimeVals[value4]
            equilTempVals[value5]=fullTempVals[value4]
            equilPotVals[value5]=fullPotVals[value4]
            equilKinVals[value5]=fullKinVals[value4]
            equilEnergyVals[value5]=fullEnergyVals[value4]
            equilPressVals[value5]=fullPressVals[value4]
            value5=value5+1
        value4=value4+1
        
    for eachValue in fullTimeVals:
        if fullTimeVals[value6]>twentyPercent:
            averageTimeVals[value7]=fullTimeVals[value6]
            averageTempVals[value7]=fullTempVals[value6]
            averagePotVals[value7]=fullPotVals[value6]
            averageKinVals[value7]=fullKinVals[value6]
            averageEnergyVals[value7]=fullEnergyVals[value6]
            averagePressVals[value7]=fullPressVals[value6]
            value7=value7+1
        value6=value6+1
    
    fullPressVals = fullPressVals*0.0001
    mainPressVals = mainPressVals*0.0001
    eosPressVals = eosPressVals*0.0001
    equilPressVals = equilPressVals*0.0001
    averagePressVals = averagePressVals*0.0001
    
    eosPressStanDev = np.std(eosPressVals)
        
    eosTemp = (sum(eosTempVals))/(len(eosTempVals))
    eosPot = (sum(eosPotVals))/(len(eosPotVals))
    eosKin = (sum(eosKinVals))/(len(eosKinVals))
    eosEnergy = (sum(eosEnergyVals))/(len(eosEnergyVals))
    eosPress = (sum(eosPressVals))/(len(eosPressVals))
    averageTemp = (sum(averageTempVals))/(len(averageTempVals))
    averagePot = (sum(averagePotVals))/(len(averagePotVals))
    averageKin = (sum(averageKinVals))/(len(averageKinVals))
    averageEnergy = (sum(averageEnergyVals))/(len(averageEnergyVals))
    averagePress = (sum(averagePressVals))/(len(averagePressVals))
    
    fullTimeVals = (fullTimeVals-15.0)*0.001
    mainTimeVals = (mainTimeVals-15.0)*0.001
    equilTimeVals = (equilTimeVals-15.0)*0.001
    
    first50=np.zeros(50)
    last50=np.zeros(50)
    
    counter1=0
    counter2=0
    counter3=-51
    
    for eachVal in first50:
        first50[counter1]=mainTempVals[counter1]
        counter1=counter1+1
        
    for eachVal in last50:
        last50[counter2]=mainTempVals[counter3]
        counter2=counter2+1
        counter3=counter3-1
    
    first50average=sum(first50)/len(first50)
    last50average=sum(last50)/len(last50)
    
    tempchange=last50average-first50average
    
    if tempchange<500.0:
        waitingtime=1.0
        
    if tempchange>500.0:
        approxhalfwaytemp=first50average+(tempchange/2)
        halfwaytemp=mainTempVals[(np.abs(mainTempVals-approxhalfwaytemp)).argmin()]
        waitingtime=mainTimeVals[(np.abs(mainTempVals-approxhalfwaytemp)).argmin()]
        if plottype == "mainTempPlot":
            plt.plot([waitingtime, waitingtime],[0, halfwaytemp], color='firebrick', linestyle="--")
            plt.plot([0, waitingtime],[halfwaytemp, halfwaytemp], color='firebrick', linestyle="--", label="Waiting Time of "+str(round(waitingtime,3))+" ns")
            
    print("Run #"+run)
    print("-----")
    print("Waiting Time to Freeze: "+str(waitingtime)+" ns")
    print("-----")
    print("Start: "+startTimeClean)
    print("End: "+endTimeClean)
    print("Goal Initial Temperature: "+fullInitialTemp+" K")
    print("Box Volume Multiplier: "+volume+"")
    print("Random Seed #"+ranSeed)
    print("Length of Simulation: "+simLength+" Timesteps")
    print("-----")
    print('EOS Average Temperature: '+str(eosTemp)+' K')
    print('EOS Average Potential Energy: '+str(eosPot)+' eV')
    print('EOS Average Kinetic Energy: '+str(eosKin)+' eV')
    print('EOS Average Total Energy: '+str(eosEnergy)+' eV')
    print('EOS Average Pressure: '+str(eosPress)+' GPa')
    print('EOS Average Pressure Error: '+str(eosPressStanDev)+' GPa')
    print("-----")
    print('Simulation Average Temperature: '+str(averageTemp)+' K')
    print('Simulation Average Potential Energy: '+str(averagePot)+' eV')
    print('Simulation Average Kinetic Energy: '+str(averageKin)+' eV')
    print('Simulation Average Total Energy: '+str(averageEnergy)+' eV')
    print('Simulation Average Pressure: '+str(averagePress)+' GPa')
    
    if plottype == "fullTempPlot":
        plt.plot(fullTimeVals, fullTempVals, color='royalblue')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Temperature of Simulation, T (K)")
        plt.title("Plot of Temperature against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.015,1)
        plt.ylim((initialTempNum-100.0), round((max(fullTempVals)+100.0),-2))
        plt.show()
        
    if plottype == "mainTempPlot":
        plt.plot(mainTimeVals, mainTempVals, color='royalblue', label="Temperature Data")
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Temperature of Simulation, T (K)")
        plt.title("Plot of Temperature against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(0, 1)
        plt.ylim((initialTempNum-100.0), round((max(mainTempVals)+100.0),-2))
        plt.legend(loc=4)
        
        #filename=run+"_"+str(initialTempNum)+"K_"+plottype+".png"
        #plt.savefig(filename,dpi=300)
        
        plt.show()
        
    if plottype == "equilTempPlot":
        plt.plot(equilTimeVals, equilTempVals, color='royalblue')
        plt.xlabel('Time Before Simulation Start, t (ns)')
        plt.ylabel("Temperature of Simulation, T (K)")
        plt.title("Plot of Temperature against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.01, 0)
        plt.ylim((initialTempNum-100.0), round((max(equilTempVals)+100.0),-2))
        
        #filename=run+"_"+str(initialTempNum)+"K_"+plottype+".png"
        #plt.savefig(filename,dpi=300)
        
        plt.show()
        
    if plottype == "fullPotPlot":
        plt.plot(fullTimeVals, fullPotVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Potential Energy of Simulation, PE (eV)")
        plt.title("Plot of Potential Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.015,1)
        plt.show()
        
    if plottype == "mainPotPlot":
        plt.plot(mainTimeVals, mainPotVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Potential Energy of Simulation, PE (eV)")
        plt.title("Plot of Potential Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(0, 1)
        plt.show()
        
    if plottype == "equilPotPlot":
        plt.plot(equilTimeVals, equilPotVals, color='darkorange')
        plt.xlabel('Time Before Simulation Start, t (ns)')
        plt.ylabel("Potential Energy of Simulation, PE (eV)")
        plt.title("Plot of Potential Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.01, 0)
        plt.show()
        
    if plottype == "fullKinPlot":
        plt.plot(fullTimeVals, fullKinVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Kinetic Energy of Simulation, KE (eV)")
        plt.title("Plot of Kinetic Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.015,1)
        plt.show()
        
    if plottype == "mainKinPlot":
        plt.plot(mainTimeVals, mainKinVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Kinetic Energy of Simulation, KE (eV)")
        plt.title("Plot of Kinetic Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(0, 1)
        plt.show()
        
    if plottype == "equilKinPlot":
        plt.plot(equilTimeVals, equilKinVals, color='darkorange')
        plt.xlabel('Time Before Simulation Start, t (ns)')
        plt.ylabel("Kinetic Energy of Simulation, KE (eV)")
        plt.title("Plot of Kinetic Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.01, 0)
        plt.show()
        
    if plottype == "fullEnergyPlot":
        plt.plot(fullTimeVals, fullEnergyVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Total Energy of Simulation, E (eV)")
        plt.title("Plot of Total Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.015,1)
        
        #filename=run+"_"+str(initialTempNum)+"K_"+plottype+".png"
        #plt.savefig(filename,dpi=300)
        
        plt.show()
        
    if plottype == "mainEnergyPlot":
        plt.plot(mainTimeVals, mainEnergyVals, color='darkorange')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Total Energy of Simulation, E (eV)")
        plt.title("Plot of Total Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(0, 1)
        plt.show()
        
    if plottype == "equilEnergyPlot":
        plt.plot(equilTimeVals, equilEnergyVals, color='darkorange')
        plt.xlabel('Time Before Simulation Start, t (ns)')
        plt.ylabel("Total Energy of Simulation, E (eV)")
        plt.title("Plot of Total Energy against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.01, 0)
        plt.show()
        
    if plottype == "fullPressPlot":
        plt.plot(fullTimeVals, fullPressVals, color='firebrick')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Pressure of Simulation, P (GPa)")
        plt.title("Plot of Pressure against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.015,1)
        plt.show()
        
    if plottype == "mainPressPlot":
        plt.plot(mainTimeVals, mainPressVals, color='firebrick')
        plt.xlabel('Time Since Simulation Start, t (ns)')
        plt.ylabel("Pressure of Simulation, P (GPa)")
        plt.title("Plot of Pressure against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(0, 1)
        
        #filename=run+"_"+str(initialTempNum)+"K_"+plottype+".png"
        #plt.savefig(filename,dpi=300)
        
        plt.show()
        
    if plottype == "equilPressPlot":
        plt.plot(equilTimeVals, equilPressVals, color='firebrick')
        plt.xlabel('Time Before Simulation Start, t (ns)')
        plt.ylabel("Pressure of Simulation, P (GPa)")
        plt.title("Plot of Pressure against Time for Run #"+run+" ("+initialTemp+" K, "+str(round((eosPress),1))+" GPa)")
        plt.xlim(-0.01, 0)
        plt.show()    

simulationdata("001.txt", "equilTempPlot")





