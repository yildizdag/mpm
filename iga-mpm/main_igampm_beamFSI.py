import sys
import time
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
import matplotlib.pyplot as plt

from patch import *
from regular2Dmesh import *

comm = PETSc.COMM_WORLD
size = comm.getSize()
rank = comm.getRank()

# BACKGROUND MESH
mesh = patch('igampm_10x5_')
print('mesh done')

# PARTICLE INITIALIZATION
L1,h1p,L2,h2p,order = 9,0.25,2.5,0.25,1
p_n,p_c,_,_,h1p,h2p = regular2Dmesh(L1,h1p,L2,h2p,order)
r_p = np.zeros((len(p_c),2))
r_p[:,0] = 0.5*(p_n[p_c[:,1],0]+p_n[p_c[:,0],0])
r_p[:,1] = 0.5*(p_n[p_c[:,2],1]+p_n[p_c[:,1],1])

#SOLID and FLUID PARTS
r_p_solid = []
r_p_fluid = []
for i in range(0,len(r_p)):
    if ((r_p[i,0]<8.) and (r_p[i,1]>1.) and (r_p[i,1]<1.5)):
        r_p_solid.append([r_p[i][0],r_p[i][1]])
    else:
        r_p_fluid.append([r_p[i][0],r_p[i][1]])

r_p_solid = np.array(r_p_solid)
r_p_fluid = np.array(r_p_fluid)

#plt.scatter(r_p_solid[:,0],r_p_solid[:,1])
#plt.scatter(r_p_fluid[:,0],r_p_fluid[:,1])
#plt.axis('equal')
#plt.show()

Yp = 2E11
nup = 0.3
C_mat = (Yp/((1.+nup)*(1.-2.*nup)))*np.array([[1-nup,nup,0],[nup,1-nup,0],[0,0,0.5*(1-2*nup)]])

rhop_solid = 7800.
rhop_fluid = 1000.
c_fluid = 50.
K_fluid = 1.4E6
mu_fluid = 0.001
TOL = 1E-10
maxF = 100
load_time = 0.005
#
V_p_solid = h1p**2*np.ones(len(r_p_solid))
V_p_fluid = h1p**2*np.ones(len(r_p_fluid))
M_p_solid = rhop_solid*V_p_solid
M_p_fluid = rhop_fluid*V_p_fluid
v_p_solid = np.zeros(length(r_p_solid),2); %Velocity Solid
v_p_fluid = np.zeros(length(r_p_fluid),2); %Velocity Fluid
