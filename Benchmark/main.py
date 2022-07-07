import petsc4py
import slepc4py
import sys
import numpy as np
import os.path
import math
import matplotlib.pyplot as plt
from mpi4py import MPI
from readPatch2D import *
from bc import *
from igafem import *

petsc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

comm = MPI.COMM_WORLD
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
fem.fem_connectivity(2)
print(fem.LM[0][0,:,:])
#print(fem.LM[0][-1])
