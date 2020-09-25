import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('mvt.dat')
datamig = np.genfromtxt('migration.dat')

plt.figure()
#plot mcore
plt.plot(data[:,0],data[:,1])
#plot mplanet
plt.plot(data[:,0],data[:,2],'--')
#plt.xlim([0,4e6])
#plt.ylim([0,100])
plt.ylabel(r'Planet Mass, M$_{\oplus}$')
plt.xlabel('Time, yrs')
plt.legend()

plt.figure()
plt.plot(datamig[:,0],datamig[:,1])
plt.ylabel('Semi-major axis, AU')
plt.xlabel('Time, yrs')
plt.show()
