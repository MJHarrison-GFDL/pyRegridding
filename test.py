import numpy as np
import regrid
import remap
#import matplotlib.pyplot as plt

H=100.;n0=10;z0=np.linspace(0,-H,n0+1)
T_s=-2;T_b=2
S_s=34.2;S_b=35.
T=np.linspace(T_s,T_b,n0)
S=np.linspace(S_s,S_b,n0)
ps=np.array([0.])
fs=np.array([0.])
zbot=np.array([H])
CoordRes=z0[:-1]-z0[1:]

z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]

regrid.regrid_mod.regrid(z_in,T_in,S_in,zbot,ps_in,fs_in,'ZSTAR',CoordRes,'PCM',z_out)

print('z_in=',z_in)
print('z_out=',z_out)
