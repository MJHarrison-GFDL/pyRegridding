import numpy as np
import regrid
import remap
import matplotlib.pyplot as plt#import matplotlib.pyplot as plt


def write_MOM_PF(ale_coord_config="UNIFORM",ale_units="m",eos='LINEAR',interpolation_scheme='PLM'):

    file=open('MOM_input','w')
    txt="REGRIDDING_COORDINATE_UNITS = "+ale_units+" \n"
    print(txt)
    file.write(txt)
    txt='ALE_COORDINATE_CONFIG = '+ale_coord_config+" \n"
    file.write(txt)
    txt='EQN_OF_STATE = '+eos+" \n"
    file.write(txt)
    txt='INTERPOLATION_SCHEME = '+interpolation_scheme+" \n"
    file.write(txt)
    file.close()

def destroy_MOM_PF():
    import os
    if os.path.isfile('MOM_input'):
        os.remove('MOM_input')

VERBOSE=False
H=100.;n0=10;z0=np.linspace(0,-H,n0+1)
T_s=2;T_b=-2
S_s=34.2;S_b=35.
T=np.linspace(T_s,T_b,n0)
S=np.linspace(S_s,S_b,n0)
ps=np.array([1.e4])
fs=np.array([0.])
zbot=np.array([H])

print('====== ZSTAR (UNIFORM) ======')
#Uniform Resolution
CoordRes=z0[:-1]-z0[1:]
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF()
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'Z*',CoordRes,'PCM')
#print('Initial Interface Positions=',np.squeeze(z_in))
#print('Final Interface Positions=',np.squeeze(z_out))
diff=z_out-z_in
print('RMS motion= ',np.std(diff))
nfig=1
if VERBOSE:
#    print('Initial Interface Positions=',np.squeeze(z_in))
#    print('Final Interface Positions=',np.squeeze(z_out))
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')

print('====== ZSTAR (FNC1) ======')
#Uniform Resolution
CoordRes=z0[:-1]-z0[1:]
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF(ale_coord_config="FNC1:1.0,100.,0,1.e-6",ale_units="m")
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'Z*',CoordRes,'PCM')
#print('Initial Interface Positions=',np.squeeze(z_in))
#print('Final Interface Positions=',np.squeeze(z_out))
diff=z_out-z_in
print('RMS motion= ',np.std(diff))

if VERBOSE:
#    print('Initial Interface Positions=',np.squeeze(z_in))
#    print('Final Interface Positions=',np.squeeze(z_out))
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')


print('====== RHO (Linear EOS;RFNC1) ======')
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
destroy_MOM_PF()
write_MOM_PF(ale_coord_config="\"RFNC1:10,1026.96,1027.12,1027.28,1.065,1028.4,0.001,1\"",ale_units="kg m^-3",eos='LINEAR',interpolation_scheme='PLM')
#write_MOM_PF(ale_coord_config="\"RFNC1:10,1026.96,1027.12,1027.28,1.12,1028.4,0.001,0\"",ale_units="kg m^-3",eos='LINEAR',interpolation_scheme='PLM')
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'RHO',CoordRes,'PCM')
diff=z_out-z_in
print('RMS motion= ',np.std(diff))

if VERBOSE:
#    print('Initial Interface Positions=',np.squeeze(z_in))
#    print('Final Interface Positions=',np.squeeze(z_out))
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')
    plt.show()
