import petsc4py
#import slepc4py
import sys
import numpy as np
import math
import time
#from mpi4py import MPI
from igafem import *
from readPatch2D import *
from bc import *
from nurbs import *
from model import *

from petsc4py import PETSc
#from slepc4py import SLEPc

petsc4py.init(sys.argv)

#comm = PETSc.COMM_WORLD
##############
# PATCH FAMILY
##############
fem = igafem('lindholm_200_')
print("##################")
print("IGAFEM PATCH INFO:")
print("Patch Orders:",fem.p[0])
print("Number of Elements:",fem.nel[0])
print("Number of Control Points:",fem.nnp[0])
print("Local Stiffness Size:",fem.nen[0])
print("##################")
fem.zeroDirichlet(0,'Y','GT',0.99)
fem.fem_connectivity(3)
############################
x_gp, w_gp = fem.quad2D(3)
K = PETSc.Mat().create(comm=None)
K.setType('seqaij')
K.setSizes([fem.dof[0],fem.dof[0]])
K.setFromOptions()
K.setPreallocationNNZ((150,150))
start_time = time.time()
for i in range(0,fem.nel[0]):
    print(i)
    iu = fem.INC[0][fem.IEN[0][i,0],0]
    iv = fem.INC[0][fem.IEN[0][i,0],1]
    el_CP, el_w = fem.getCP(0,i)
    el_LM = np.reshape(np.flipud(fem.LM[0][i,:,:]),3*fem.nen[0])
    K_loc = fem.local_stiffness(iu,iv,el_CP,el_w,x_gp,w_gp)
    for j in range(0,3*fem.nen[0]):
        if el_LM[j] != -1:
            for k in range(0,3*fem.nen[0]):
                if el_LM[k] != -1:
                    K.setValues(el_LM[j],el_LM[k],K_loc[j,k],addv=True)
K.assemblyBegin()
K.assemblyEnd()
end_time = time.time()
print('Total Time=',end_time-start_time)



#K = fem.global_stiffness(3)
#M = fem.global_mass(3)

#E = SLEPc.EPS().create(comm=comm)
#E.setProblemType(SLEPc.EPS.ProblemType.GHEP)
#E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
#E.setDimensions(24)
#E.setOperators(K,M)
#E.setFromOptions()
#E.setTolerances(tol=1e-8, max_it=5000)
#E.solve()
#print(E.getConverged())
#print(E.getEigenvalue(0))
#print(E.getEigenvalue(1))
#print(E.getEigenvalue(2))
