import numpy as np 
import matplotlib.pyplot as plt 

mlms = np.array([0.091, 0.345,1.83,7.24,44.7,147,794,3259,14728])
fft = np.array([0.054,0.126,0.193,0.875,4.7,12.7,75.7,511,1852])
naive = np.array([0.136,2.19,36.6,589,9435,150970,2548218])
g = np.array([8,16,32,64,128,256,512,1024,2048])
g2 = np.array([8,16,32,64,128,256,512])
g = g**2
g2 = g2**2
fig = plt.gcf()
#fig.set_size_inches(4,3)
fig.set_dpi(300)
plt.plot(g2,g2*g2*1e-4, 'k--')
plt.plot(g2,naive)
plt.plot(g,mlms)
plt.plot(g,g*np.log(g)*1e-4, 'k:')
plt.yscale("log")
plt.xscale("log")
plt.ylabel("Time in ms", fontsize=16)
plt.xlabel("Gridsize", fontsize=16)
plt.plot(g,fft)
plt.legend([r'O$(N^2)$',"Naive","MLMS",r"O$(N \ln{(N)})$","FFT"])
plt.tight_layout()
plt.savefig('log_fft_mlms.pdf')
plt.show()
