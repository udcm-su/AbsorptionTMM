"""
Created on Mon Feb 11 11:20:06 2019
An example where we consider the layers 
Pt___Co_______Cr_____MgO_______________
We want to show that Co gets uniformly heated
@author: lUKAS
"""
import numpy as np
from matplotlib import pyplot as plt
import TMM_abs as atmm 

#layers:  Air   Pt          Co                  Cr          MgO     Air
n_list = [1, 1.0433+3.0855j,1.0454+3.2169j,2.0150+2.8488j,1.766, 1] 
d_list = [np.inf,3,15,5,2000, np.inf] #in nm 


th0 = np.pi/4
lam0 = 400
pol = 'p'
plotpoints = 100000



[absorp,grid] = atmm.absorption(th0,lam0,n_list,d_list,pol,plotpoints)

plt.figure()
plt.suptitle('Local absorbtion profile of multi layer thin film', fontsize=12)
plt.title(r"$\lambda=400$nm, $\theta_0=\frac{\pi}{4}$, polarization =$p$",fontsize=10)
plt.xlabel("Depth of material in nm"); plt.ylabel("Normalized Absorbtion profile")
plt.xlim(0,30)
plt.plot(grid,1*absorp)
#plt.legend()
plt.show()

[M,M_n,t,r,T,R,A,theta] = atmm.TM(th0,lam0,n_list,d_list,pol) 


