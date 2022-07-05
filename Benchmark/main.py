import petsc4py
import slepc4py
import sys
import numpy as np
import os.path
import math
import matplotlib.pyplot as plt
from mpi4py import MPI
from patch import *

petsc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

comm = MPI.COMM_WORLD
##############
# PATCH FAMILY
##############
pf = patchFamily2D('igabem_DR05_4th_768_',2)
print("##################")
print("PATCH FAMILY INFO:")
print("Number of Patches = ",pf.np)
