import numpy as np
from readPatch2D import *

class bc(readPatch2D):
    def __init__(self,file_ID):
        super().__init__(file_ID)
        self.bcZero_no = []
        self.bcZero = []
    def zeroDirichlet(self,p_no,axis,cond,coord):
        self.bcZero_no.append(p_no)
        if axis == 'X':
            if cond == 'LT':
                ind = np.where(self.CP[p_no][:,0]<coord)
            elif cond == 'GT':
                ind = np.where(self.CP[p_no][:,0]>coord)
        if axis == 'Y':
            if cond == 'LT':
                ind = np.where(self.CP[p_no][:,1]<coord)
            elif cond == 'GT':
                ind = np.where(self.CP[p_no][:,1]>coord)
        if axis == 'Z':
            if cond == 'LT':
                ind = np.where(self.CP[p_no][:,2]<coord)
            elif cond == 'GT':
                ind = np.where(self.CP[p_no][:,2]>coord)
        self.bcZero.append(ind[0])
