import numpy as np 
import matplotlib.pyplot as plt 

means = np.loadtxt("res_l2", delimiter="\t")
t_min2 = means[:,0]
t_max2 = means[:,1]
t_mean2 = means[:,2]

'''means2 = np.loadtxt("res_l2_corr", delimiter="\t")

t_min3 = means[:,0]
t_max3 = means[:,1]
t_mean3 = means[:,2]
print(t_mean2 == t_mean3)'''
plt.yscale("log")
plt.xscale("log")

fig = plt.gcf()
#fig.set_size_inches(4,3)
fig.set_dpi(300)

t = np.array([8,16,32,64,128,256,512,1024,2048])# ,4096])
t = 1/np.sqrt(t)
#plt.plot(t, t_mean3)
plt.title("L2 norm")
plt.plot(t, t_mean2,'b')
plt.plot(t, t_mean2,'b^')
plt.ylabel("L2-norm", fontsize=16)
plt.xlabel(r'$\frac{1}{\sqrt{N}}$', fontsize=16)
plt.legend(["L2"])
plt.tight_layout()
plt.savefig('l2_err.pdf')  
plt.show()
