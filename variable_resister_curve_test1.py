import numpy as np
import matplotlib.pyplot as plt

mattheta=np.linspace(0,2*np.pi,100)

Rc=1e6
Ro=100e3
a=Rc/Ro*0.3
matB=np.logspace(np.log10(0.1),np.log10(1.6),5)

matRv=[]
matinvG=[]
for theta in mattheta:
    #Rv=Rc/((a-Rc/Ro)/(2*np.pi)*theta+Rc/Ro)-Ro
    Rv=Ro*(1/(1-theta/(2*np.pi)+a/(2*np.pi)*Ro/Rc*theta)-1)
    invG=1+Rc/(Rv+Ro)

    matRv.append(Rv)
    matinvG.append(invG)

mattheta=np.array(mattheta)
fig, ax = plt.subplots(nrows=2, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
ax[0,0].plot(mattheta*180/np.pi,matRv,"k-")
ax[1,0].plot(mattheta*180/np.pi,matinvG,"k-")

for B in matB:
    matRv_approx=[]
    matinvG=[]
    for theta in mattheta:
        Rv_approx=Ro*(Rc/(a*Ro)-1)*(np.exp(theta*B)-1)/(np.exp(2*np.pi*B)-1)
        invG=1+Rc/(Rv_approx+Ro)
        matRv_approx.append(Rv_approx)
        matinvG.append(invG)
    ax[0,0].plot(mattheta*180/np.pi,matRv_approx,"r-")
    ax[1,0].plot(mattheta*180/np.pi,matinvG,"r-")

ax[1,0].set_xlabel("theta [deg]")
ax[0,0].set_ylabel("Rv [Ohm]")
ax[1,0].set_ylabel("1/G [-]")
plt.show()
