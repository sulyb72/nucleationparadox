import numpy as np
import matplotlib.pyplot as plt

nvtend = np.loadtxt('001rdf2.txt', skiprows=4)
nveend = np.loadtxt('001rdf3.txt', skiprows=49603)

nvtenddistance = nvtend[:,1]
nvtendgr1 = nvtend[:,2]
nvtendgr2 = nvtend[:,4]
nvtendgr3 = nvtend[:,6]

nveenddistance = nveend[:,1]
nveendgr1 = nveend[:,2]
nveendgr2 = nveend[:,4]
nveendgr3 = nveend[:,6]

plt.plot(nvtenddistance, nvtendgr1, color='limegreen', label='iron-iron')
plt.plot(nvtenddistance, nvtendgr2, color='darkviolet', label='silicon-silicon')
plt.plot(nvtenddistance, nvtendgr3, color='darkorange', label='iron-silicon')
plt.xlim(0, 6)
plt.ylim(0, 30)
plt.xlabel('Distance, r (å)')
plt.ylabel('Correlation function, g(r)')
plt.legend()

filename="radial_distribution_functions_1.png"
plt.savefig(filename,dpi=300)

plt.show()

plt.plot(nveenddistance, nveendgr1, color='limegreen', label='iron-iron')
plt.plot(nveenddistance, nveendgr2, color='darkviolet', label='silicon-silicon')
plt.plot(nveenddistance, nveendgr3, color='darkorange', label='iron-silicon')
plt.xlim(0, 6)
plt.ylim(0, 30)
plt.xlabel('Distance, r (å)')
plt.ylabel('Correlation function, g(r)')
plt.legend()

filename="radial_distribution_functions_2.png"
#plt.savefig(filename,dpi=300)

plt.show()