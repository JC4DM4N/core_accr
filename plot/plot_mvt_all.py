import numpy as np
import matplotlib.pyplot as plt

colors=['blue','green','red','orange']

for k,i in enumerate([5,10,20,30]):
	plt.figure(k)
	for l,j in enumerate([7.5,10,15]):
		data = np.genfromtxt('mvt_r%i_sig%s.dat' %(i,str(j)))
		#anything above isolation mass=np.nan
		data[:,1][data[:,1]==data[:,1].max()] = np.nan
		data[:,2][data[:,2]==data[:,2].max()] = np.nan
		#plot mcore
		plt.plot(data[:,0],data[:,1],
				 label=r'$\Sigma_{p, \rm 5 AU} = %s$gcm$^{-2}$' %str(j),
				 c=colors[l])
		#plot mplanet
		plt.plot(data[:,0],data[:,2],'--',c=colors[l])
		try:
			trunaway = data[data[:,2]>2*data[:,1],0][0]
			trunaway = data[data[:,2]>1270][0,0]
			print('Runaway growth time for %i AU %.3f gcm-2: %.3f' %(i, j,
																	trunaway))
		except:
			print('Planet at %i AU %.3f gcm-2 does not reach ruanway' %(i,j))																	
	plt.xlim([0,10e6])
	if i==50:
		plt.ylim([0,1])
		plt.text(1e5,0.9,r'$R=$%s AU' %str(i), fontsize=20)
	else:
		plt.ylim([0,100])
		plt.text(1e5,90,r'$R=$%s AU' %str(i), fontsize=20)
	plt.xticks([0,0.2e7,0.4e7,0.6e7,0.8e7,1e7],labels=['0','2','4','6','8','10'])
	plt.ylabel(r'Planet Mass, M$_{\oplus}$',fontsize=12)
	plt.xlabel('Time, Myrs',fontsize=12)
	#plt.title(r'$M_{\rm core,init}=0.1$M$_{\oplus}$, $\Sigma \propto R^{-1}$' +
	#		   ' \n Disc profile and evolution is taken from visag runs')
	plt.legend(loc='upper right')
	plt.grid()
plt.show()
