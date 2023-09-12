import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def conf_plt(logx=False,logy=False):
  fig, ax = plt.subplots()
  fig.set_size_inches(3,4)
  #fig = plt.gcf()
  fig.set_dpi(300)
  if(logy):
    ax.set_yscale("log")
  if logx:
    ax.set_xscale('log')
  return ax, fig

#arr = np.fromfile("pressure_4", sep='\t')
#arr = np.fromfile("initial_press", sep='\t')
#arr = np.fromfile("part_res_93", sep='\t')
arr = np.fromfile("res_mlms2048", sep='\t')
n = int(np.sqrt(arr.size))
arr2 = arr.reshape(n,n)

ax,fig = conf_plt(0,0)
surf = ax.imshow(arr2,cmap='viridis')
#surf = ax.imshow(arr2,cmap='Greys')

ax.invert_yaxis()
plt.colorbar(surf, shrink=0.6, aspect=6)
plt.tight_layout()
#plt.savefig("coarse_press.pdf",bbox_inches='tight')
#plt.savefig("init_press.pdf",bbox_inches='tight')
#plt.savefig("part_res_coarse.pdf",bbox_inches='tight')
plt.savefig("res_mlms2k.pdf",bbox_inches='tight')
plt.show()
