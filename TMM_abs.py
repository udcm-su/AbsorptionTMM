"""
Created on Thu Jan 31 16:29:36 2019
Transfer matrix representation to calculate local absorbtion. 
Based on and compared with the work of: 
Steven J. Byrnes. Multilayer optical calculations.arXiv:1603.02720v3,2018
@author: lUKAS
"""

import numpy as np


def fresnel(theta_in,n_in,n_out,pol): 
    n_in = complex(n_in); n_out = complex(n_out)
    theta_out = np.arcsin(n_in*np.sin(theta_in)/n_out)
    if pol == 's': 
        rs = (n_in*np.cos(theta_in) - n_out*np.cos(theta_out))/\
             (n_in*np.cos(theta_in) + n_out*np.cos(theta_out))
        ts = 2*n_in*np.cos(theta_in)/(n_in*np.cos(theta_in)+n_out*np.cos(theta_out))
        return(theta_out,rs,ts)
    if pol == 'p':
        rp = (n_out*np.cos(theta_in)-n_in*np.cos(theta_out))/\
             (n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
        tp = 2* n_in*np.cos(theta_in)/(n_out*np.cos(theta_in)+n_in*np.cos(theta_out))
        return(theta_out,rp,tp) 


def TM(theta_in,lambda0,n_vec,d_vec,pol): 
    #create complex arrays for variables
    theta   = np.zeros(len(n_vec), dtype = complex); theta[0] = theta_in
    phi     = np.zeros(len(n_vec),dtype = complex)
    rn      = np.zeros(len(n_vec)-1, dtype = complex) 
    tn      = np.zeros_like(rn,dtype = complex)
    M_n     = np.zeros((len(n_vec),2,2), dtype = complex)
    M       = np.eye(2,dtype = complex)
    for i in range(len(n_vec)-1): # to obtian all angels/rn/tn for each layer
        [theta[i+1],rn[i],tn[i]] = fresnel(theta[i],n_vec[i],n_vec[i+1],pol)
    #M = M0*M1*M2*M4*....
    for k in range(1,len(n_vec)-1):#loop over all interfaces except 1st
        phi[k]  = 2*np.pi*n_vec[k]*np.cos(theta[k])*d_vec[k]/lambda0
        Tn      = np.array([[np.exp(-1j*phi[k]),0],[0,np.exp(1j*phi[k])]],dtype = complex)/tn[k]
        Pn      = np.array([[1,rn[k]],[rn[k],1]],dtype = complex)
        M_n[k]     = np.dot(Tn,Pn)
        M = np.dot(M,M_n[k])
    #compute for the first interface: 
    trans0 = np.array([[1,rn[0]],[rn[0],1]],dtype= complex)/tn[0]
    M = np.dot(trans0,M)
    #Complex transmission/reflection amplitude
    t = 1/M[0,0]
    r = M[1,0]/M[0,0]
    #Fraction of power transmitted
    if pol == 's': #s-polarized
        T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(theta[-1]))/\
            np.real(n_vec[0]*np.cos(theta[0]))
    elif pol == 'p': #p-polarized
        T = np.abs(t)**2*np.real(n_vec[-1]*np.cos(np.conj(theta[-1])))/\
            np.real(n_vec[0]*np.cos(np.conj(theta[0])))
    #Fraction of power reflected
    R = np.abs(r)**2
    A = 1.-T-R
    return(M,M_n,t,r,T,R,A,theta)

def layerAmpl(theta_in,lambda0,n_vec,d_vec,pol): 
    """
    After r & t have been calculated and all the respective matrices M_n
    for each layer are known, we can go 'backwards', i.e. from the last to the
    first layer, and compute all the amplituedes for the forward v_n and 
    backward w_n traveling wave. -> [v_n,w_n].T = M_n @ [v_{n+1},w_{n+1}].T
    """
    [M,M_n,t,r,T,R,A,theta] = TM(theta_in,lambda0,n_vec,d_vec,pol)
    vw_list = np.zeros((len(n_vec),2),dtype = complex)
    vw =np.array([[t],[0]])
    vw_list[-1,:] = vw.T
    for i in range(len(n_vec)-2,0,-1):
        vw = np.dot(M_n[i],vw)
        vw_list[i,:] = vw.T
    return(vw_list,theta)
    
def absorption(theta_in,lambda0,n_vec,d_vec,pol,points):
    #reload the forward and backward wave coefficients for every layer
    [vw_n,theta]  = layerAmpl(theta_in,lambda0,n_vec,d_vec,pol)
    total_len = np.sum(d_vec[1:-1])
    pointcount = 0
    #a is an array where the normalized absorbtion for the entire grid is stored 
    a = []
    for i in range(1,len(n_vec)-1):
        kz      = 2*np.pi*n_vec[i]*np.cos(theta[i])/lambda0
        #How many points, out of the total 'points' does each layer get, with respect to 
        #the length of the layer
        points_per_layer = int(np.round(points*d_vec[i]/total_len))
        #for every layer, the space grid is reset. I.e. every layer starts at 0
        layer   = np.linspace(0,d_vec[i],points_per_layer)
        v = vw_n[i,0]; w = vw_n[i,1];#complex wave amplitudes for every layer
        Ef = v * np.exp(1j * kz * layer)#forward traveling wave 
        Eb = w * np.exp(-1j * kz *layer)#backward traveling wave
        if pol == 'p':#p-polarized
            a_layer = (n_vec[i]*np.conj(np.cos(theta[i]))*(kz*np.abs(Ef-Eb)**2-np.conj(kz)*np.abs(Ef+Eb)**2)).imag /\
                    (n_vec[0]*np.conj(np.cos(theta[0]))).real
        if pol == 's': 
            a_layer = (np.abs(Ef+Eb)**2 *np.imag(kz*n_vec[i]*np.cos(theta[i])))/\
                    (np.real(n_vec[0]*np.cos(theta[0])))
        #for every layer calculate the absorbtion grid and append it to the total 
        a   = np.append(a,a_layer)
        #to get the right amount of points considered in the grid, since we round.
        pointcount += points_per_layer 
    a = np.real(a)
    grid = np.linspace(0,total_len,pointcount)
    return(a,grid)
    
   
#the analytic absorbtion should give the same result as the absorption function
#this is true for s- polarized light but for some reason there is a phase shift
# when the p- polarized light is considered.... 
def absorption_analytic(theta_in,lambda0,n_vec,d_vec,pol,points):
    [vw_n,theta]  = layerAmpl(theta_in,lambda0,n_vec,d_vec,pol)
    total_len = np.sum(d_vec[1:-1])
    grid            = np.linspace(0,total_len,points) 
    #create a grid along z and space for absorbtion vector
    a               = []
    for i in range(1,len(n_vec)-1):
        kz      = 2*np.pi*n_vec[i]*np.cos(theta[i])/lambda0  
        if pol == 's':#s-polarized
            normal = np.real(n_vec[0]*np.cos(theta[0]))
            A1 = np.abs(vw_n[i,1])**2*np.imag(n_vec[i]*np.cos(theta[i])*kz)/normal #w^2
            A2 = np.abs(vw_n[i,0])**2*np.imag(n_vec[i]*np.cos(theta[i])*kz)/normal#v^2
            A3 = vw_n[i,0]*np.conj(vw_n[i,1])*np.imag(n_vec[i]*np.cos(theta[i])*kz)/normal#vw*
        if pol == 'p':#p-polarized
            normal = np.real(n_vec[0]*np.conj(np.cos(theta[0])))
            A1 = (np.abs(vw_n[i,1])**2*2*np.imag(kz)*np.real(n_vec[i]*np.cos(np.conj(theta[i]))))/normal
            A2 = (np.abs(vw_n[i,0])**2*2*np.imag(kz)*np.real(n_vec[i]*np.cos(np.conj(theta[i]))))/normal
            A3 = vw_n[i,0]*np.conj(vw_n[i,1])*2*np.real(kz)*np.imag(n_vec[i]*np.cos(np.conj(theta[i])))/normal

        layer   = np.linspace(0,d_vec[i],int(np.round(points*d_vec[i]/total_len)))
        a_layer = A1*np.exp(2*layer*np.imag(kz))+\
                    A2*np.exp(-2*layer*np.imag(kz))+\
                    A3*np.exp(2*1j*layer*np.real(kz))+\
                    np.conj(A3)*np.exp(-2*1j*layer*np.real(kz))
        a       = np.append(a,a_layer)
    a = np.real(a)
    return(a,grid)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    