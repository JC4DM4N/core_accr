import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('densities.dat')

au = 1.496e13

msun = 2e33
mgas = 0 
mdust = 0
for i in range(1,len(data)):
    mgas += np.pi*(data[i,0]**2 - data[i-1,0]**2)*data[i,1]
    mdust += np.pi*(data[i,0]**2 - data[i-1,0]**2)*data[i,2]

print('mgas = %.3f , mdust = %.3f , dust-to-gas = %.3f' %(mgas/msun, mdust/msun, mdust/mgas))

plt.plot(data[:,0]/au, data[:,2], label='Dust')
plt.plot(data[:,0]/au, data[:,1], label='Gas')
plt.grid(alpha=0.25)
plt.ylabel(r'Planetesimal surface density, gcm$^{-2}$')
plt.xlabel(r'Radius, AU')
plt.xlim([0,500])
plt.ylim([0,25])

plt.show()
