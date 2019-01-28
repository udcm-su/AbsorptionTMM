
"""
Created on Mon Jan 28 11:49:27 2019
Usint Transfer matrix representation to calculate the local absorbtion a(z)
in a multi layer thin film. 
Each of the main functions correspond to one of the 5 steps described in 
the README provided at GitHub. 
@author: lUKAS
"""

from matplotlib import pyplot as plt
import numpy as np
import numericalunits as u
u.reset_units('SI')



lamb = 800#in nm
theta = 0 #in degree with respect to normal of suface 
#refractiv index
n_800nm = np.array([1, 3.896 + 0.109j, 0.18 + 4.7j, 1.111 + 3.480j, 1.54, 3.681 + 0.005j, 1])
thickness = np.array([0, 1,2,5,8,6, 0])



def degr2rad(degr):
    return(degr*np.pi/180)

def fresnel(theta_in,n_in,n_out,polarization=0):
    """
    0 is unpolarized, 1 is p-polarized, 2 is s-polarized
    """
    #map theta to radiants
    theta_in = degr2rad(complex(theta_in))
    n_in = complex(n_in); n_out = complex(n_out)
    #Senllius law
    theta_out = complex(np.arcsin(n_in*np.sin(theta_in)/n_out))
    #here maybe change - to +
    rp = -(n_out*np.cos(theta_in)-n_in*np.cos(theta_out))/(n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
    tp = 2*n_in*np.cos(theta_in)/(n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
    rs = (n_in*np.cos(theta_in)-n_out*np.cos(theta_out))/(n_in*np.cos(theta_in)+n_out*np.cos(theta_out))
    ts = 2*n_in*np.cos(theta_in)/(n_in*np.cos(theta_in)+n_out*np.cos(theta_out))
    #Theta in RAD
    if polarization==0:
    	return (180*theta_out/np.pi, 0.5*(rs+rp), 0.5*(ts+tp))
    elif polarization==1:
    	return (180*theta_out/np.pi, rp, tp)
    elif polarization==2:
    	return (180*theta_out/np.pi, rs, ts)
    
    
def TM(theta_in,lambda0,n_vec,d_vec,polarization):
    #Create space for variables
    theta   = np.zeros(len(n_vec),dtype=complex); theta[0] = theta_in
    phi     = np.zeros(len(n_vec)-1,dtype=complex)
    rn      = np.zeros(len(n_vec)-1,dtype=complex)
    tn      = np.zeros(len(n_vec)-1,dtype=complex)
    Tr      = np.ones((len(n_vec)-1,2,2),dtype=complex)
    Pr      = np.zeros((len(n_vec)-1,2,2),dtype=complex)
    Mr      =  np.zeros((len(n_vec)-1,2,2),dtype=complex)
    for i in range(0,len(n_vec)-1):
        [theta[i+1],rn[i],tn[i]] = fresnel(theta[i],n_vec[i],n_vec[i+1],polarization)
        phi[i] = (2*np.pi*n_vec[i+1]*d_vec[i+1]*np.cos(degr2rad(theta[i+1])))/lambda0
        Tr[i][0][1] = rn[i]; Tr[i][1][0] = rn[i]; Tr[i]/= tn[i]
        Pr[i][0][0] = np.exp(-1j*phi[i]); Pr[i][1][1] = np.exp(1j*phi[i])
        #M0 = Tr0*Pr0; M1 = Tr1*Pr1 .... =>M
        "Here the order matters! Tr*Pr vs Pr*Tr"
        Mr[i] = np.dot(Tr[i],Pr[i]) 
    #Mr[0] = np.array([[1,rn[0]],[rn[0],1]])/tn[0]
    return(rn,tn,phi,theta,Tr,Pr,Mr)
    
# =============================================================================
# #Multiplication from right to left. I.e.M = M0*M1*M2.... <-- 
# def TRA(theta_in,lambda0,n_vec,d_vec,polarization):
#     [rn, tn,phi,Tr,Pr,Mr] = TM(theta_in,lambda0,n_vec,d_vec,polarization)
#     #Multiply all matrices M0 = Pr0*Tr0 .... =>M = M0*M1*M2...
#     M =  np.dot(Mr[-2],Mr[-1])
#     l = np.shape(Mr)[0]
#     print(l)
#     for i in range(0,l-2):
#         M = np.dot(Mr[-3-i],M)
#     t = 1/M[0,0]; r = M[1,0]/M[0,0]
#     return(t,r)
# =============================================================================


def TRA(theta_in,lambda0,n_vec,d_vec,polarization):
    [rn, tn,phi,theta,Tr,Pr,Mr] = TM(theta_in,lambda0,n_vec,d_vec,polarization)
    M = np.eye(2)
    #Multiplication from left to right. I.e. starting with M = --> M0*M1*M2....
    for i in range(0,np.shape(Mr)[0]): 
        M = np.dot(M,Mr[i])
    #Transmission Amplitude
    t = 1/M[0,0]; r = M[1,0]/M[0,0]
    #Fraction of power transmitted
    if polarization == 2: #s-polarized
        T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(degr2rad(theta[-1])))/\
            np.real(n_vec[0]*np.cos(degr2rad(theta[0])))
    elif polarization == 1: #p-polarized
        T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(np.conj(degr2rad(theta[-1]))))/\
            np.real(n_vec[0]*np.cos(np.conj(degr2rad(theta[0]))))
    #Fraction of power reflected
    R = np.abs(r)**2
    A = 1.-T-R   
    return(M,t,r,T,R,A)
    
def layerAmpl(theta_in,lambda0,n_vec,d_vec,polarization): 
    [rn, tn,phi,theta,Tr,Pr,Mr] = TM(theta_in,lambda0,n_vec,d_vec,polarization)
    [M,t,r,T,R,A] = TRA(theta_in,lambda0,n_vec,d_vec,polarization)
    vw_n = np.zeros((2,len(n_vec)),dtype=complex)
    vw_n[:,0] = np.array([1,r]).T
    #After getting the overal reflection r and transmission t
    #Solve a system of equations to get v and w (forward/backward)
    #coefficients for every layer
    for i in range(1,len(n_vec)-1): 
        vw_n[:,i] = np.linalg.solve(Mr[i-1],vw_n[:,i-1])
    return(vw_n,theta)
    
def absorbtion(theta_in,lambda0,n_vec,d_vec,polarization,grid): 
    [vw_n,theta] = layerAmpl(theta_in,lambda0,n_vec,d_vec,polarization)
    accum = np.cumsum(d_vec)
    a = np.zeros(len(grid),dtype=complex)
    for i in range(0,len(n_vec)):
        kz = 2*np.pi*n_vec[i]*np.cos(degr2rad(theta[i]))/lambda0        
        if polarization == 2:#s-polarized
            normal = np.real(n_vec[0]*np.cos(degr2rad(theta[0])))
            A1 = np.abs(vw_n[1,i])**2*np.imag(n_vec[i]*np.cos(degr2rad(theta[i]))*kz)/normal #w^2
            A2 = np.abs(vw_n[0,i])**2*np.imag(n_vec[i]*np.cos(degr2rad(theta[i]))*kz)/normal#v^2
            A3 = vw_n[0,i]*np.conj(vw_n[1,i])*np.imag(n_vec[i]*np.cos(degr2rad(theta[i]))*kz)/normal#vw*
        if polarization == 1:#p-polarized
            normal = np.real(n_vec[0]*np.cos(np.conj(degr2rad(theta[0]))))
            A1 = np.abs(vw_n[1,i])**2*2*np.imag(kz)*np.real(n_vec[i]*np.cos(np.conj(degr2rad(theta[i]))))/normal
            A2 = np.abs(vw_n[0,i])**2*2*np.imag(kz)*np.real(n_vec[i]*np.cos(np.conj(degr2rad(theta[i]))))/normal
            A3 = vw_n[0,i]*np.conj(vw_n[1,i])*2*np.real(kz)*np.imag(n_vec[i]*np.cos(np.conj(degr2rad(theta[i]))))/normal
        left = np.where(grid >=accum[i])[0]; left = left[0]
        right = np.where(grid<=accum[i+1])[0]; right = right[-1]
        #absorbtion profile for each layer
        a[left:right] = A1*np.exp(2*grid[left:right]*np.imag(kz))+\
                    A2*np.exp(-2*grid[left:right]*np.imag(kz))+\
                    A3*np.exp(2*1j*grid[left:right]*np.real(kz))+\
                    np.conj(A3)*np.exp(-2*1j*grid[left:right]*np.real(kz))
    a = np.real(a)
    return(a)
    

#vw_n = layerAmpl(theta,lamb,n_800nm,thickness,2)   
[M,t,r,T,R,A]=TRA(theta,lamb,n_800nm,thickness,2)
#grid = np.linspace(0,np.cumsum(thickness)[-1],100)
#absorb = absorbtion(theta,lamb,n_800nm,thickness,2,grid)


n_vector = [2.3+0.1j,3.3+0.1j]
d_test = [0,100,200]
grid_test = np.linspace(0,np.cumsum(d_test)[-1],100)
lambda0 = 800#nm
theta = 0;
absorb = absorbtion(theta,lambda0,n_vector,d_test,2,grid_test)
plt.figure()
plt.plot(grid_test,absorb)

#Test, if M*[t,0].T = [1,r].T
#print(np.dot(M,np.array([t,0]).T))



    
    
    
    
    