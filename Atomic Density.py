import numpy as np
import matplotlib.pyplot as plt

file = "001step2posfinal.txt"
position="start"
data = np.loadtxt(file)
atomNumData = data[:,0]
atomTypeData = data[:,1]
xData = data[:,2]
yData = data[:,3]
zData = data[:,4]

binsNum=250
dataSet=zData
rangeMin=min(dataSet)
rangeMax=max(dataSet)

#plt.hist(dataSet, bins=binsNum, range=(rangeMin,rangeMax), color="royalblue", histtype='step')
histo, binedges = np.histogram(dataSet, bins=binsNum, range=(rangeMin,rangeMax))
multipliedhisto=histo/68732.6312980447
plt.hist(binedges[:-1], binedges, weights=multipliedhisto, histtype='step', color="royalblue")
if position=="start":
    plt.title("Atomic Density at Beginning Timestep of Simulation")
if position=="end":
    plt.title("Atomic Density at Final Timestep of Simulation")
plt.xlabel("Z-Axis Position, z (Å)")
plt.ylabel("Density, ρ (Å⁻³)")
plt.xlim(rangeMin,rangeMax)
plt.ylim(0, max(multipliedhisto))

plt.gcf().subplots_adjust(left=0.15)

filename="atomicdensity_"+position+".png"
#plt.savefig(filename,dpi=300)
    
plt.show()