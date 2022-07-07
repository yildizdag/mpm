import numpy as np
from readPatch2D import *
from bc import *

class igafem(bc,readPatch2D):
    def __init__(self,file_ID):
        super().__init__(file_ID)

    def fem_connectivity(self,ndof):
        count = 0
        for i in range(0,self.np):
            nnp_p = self.nnp[i]
            ID_p = np.zeros((nnp_p,ndof)).astype(int)
            for j in range(0,nnp_p):
                dirichlet = False
                for p in range(0,len(self.bcZero_no)):
                    if i == self.bcZero_no[p]:
                        if j not in self.bcZero[p]:
                            dirichlet = False
                        else:
                            dirichlet = True
                    if dirichlet == True:
                        break
                for k in range(0,ndof):
                    if dirichlet == False:
                        ID_p[j,k] = count
                        count += 1
                    elif dirichlet == True:
                        ID_p[j,k] = -1
            self.ID.append(ID_p)
        # LM - FEM CONNECTIVITY
        for i in range(0,self.np):
            LM_p = np.zeros(((self.nel[i],self.nen[i],ndof))).astype(int)
            for j in range(0,self.nel[i]):
                for k in range(0,self.nen[i]):
                    LM_p[j,k,:] = self.ID[i][self.IEN[i][j,k],:]
                    self.LM.append(LM_p)
