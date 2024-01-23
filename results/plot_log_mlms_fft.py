import numpy as np 
import matplotlib.pyplot as plt 
# import matplotlib
# print(matplotlib.get_cachedir())

plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = 'Helvetica'
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

mlms = np.array([0.091, 0.345,1.83,7.24,44.7,147,794,3259,14728])
fft = np.array([0.054,0.126,0.193,0.875,4.7,12.7,75.7,511,1852])
naive = np.array([0.136,2.19,36.6,589,9435,150970,2548218])
g = np.array([8,16,32,64,128,256,512,1024,2048])
g2 = np.array([8,16,32,64,128,256,512])
g = g**2
g2 = g2**2
fig = plt.gcf()
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#fig.set_size_inches(4,3)
fig.set_dpi(300)
plt.plot(g2,g2*g2*1e-4, linestyle='--')
plt.plot(g2,naive, linestyle='-.' )
plt.plot(g,mlms, dashes=[8, 2])
plt.plot(g,g*np.log(g)*1e-4, 'k:')
plt.yscale("log")
plt.xscale("log")
plt.ylabel("Time ($ms$)", fontname="Arial", fontsize=16)
plt.xlabel("Gridsize $x \cdot y$ ", fontname="Arial", fontsize=16)
plt.plot(g,fft)
plt.legend([r'O$(N^2)$',"Naive","MLMS",r"O$(N \ln{(N)})$","FFT"])
plt.tight_layout()
plt.savefig('log_fft_mlms.pdf')
plt.show()
