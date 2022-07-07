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
fem = igafem('igabem_DR05_4th_768_')
print("##################")
print("PATCH FAMILY INFO:")
print("Number of Patches = ",fem.np)
fem.zeroDirichlet(0,'X','LT',-0.64199)
fem.zeroDirichlet(0,'X','GT',0.64199)
fem.fem_connectivity(2)
print(fem.ID[0])
print(fem.LM[0][-1])
