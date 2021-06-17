import numpy as np
import regrid
import remap
import matplotlib.pyplot as plt#import matplotlib.pyplot as plt

VERBOSE=True
H=100.;n0=10;z0=np.linspace(0,-H,n0+1)
T_s=2;T_b=-2
S_s=34.2;S_b=35.
T=np.linspace(T_s,T_b,n0)
S=np.linspace(S_s,S_b,n0)
ps=np.array([1.e4])
fs=np.array([0.])
zbot=np.array([H])

#ZSTAR
# print('====== ZSTAR (UNIFORM) ======')
# #Uniform Resolution
# CoordRes=z0[:-1]-z0[1:]
# #Perturb the Surface
# z0[0]=z0[0]+0.1
# z_in = z0[np.newaxis,np.newaxis,:]
# T_in = T[np.newaxis,np.newaxis,:]
# S_in = S[np.newaxis,np.newaxis,:]
# z_out = np.zeros(z_in.shape)
# ps_in = ps[np.newaxis,:]
# fs_in = fs[np.newaxis,:]
# zbot_in=zbot[np.newaxis,:]
# z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'Z*',CoordRes,'PCM')
# print('Initial Interface Positions=',np.squeeze(z_in))
# print('Final Interface Positions=',np.squeeze(z_out))

print('====== RHO (Linear EOS) ======')
RHO_T0_S0=1000.
dRdT=-0.2
dRdS=0.8
rho=RHO_T0_S0+dRdT*T + dRdS*S
CoordRes=rho
print(CoordRes)
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'RHO',CoordRes,'PCM')
diff=z_out-z_in
print('RMS motion= ',np.std(diff))

if VERBOSE:
    print('Initial Interface Positions=',np.squeeze(z_in))
    print('Final Interface Positions=',np.squeeze(z_out))
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.show()
