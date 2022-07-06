import numpy as np
from patch import *

class igafem(readPatch2D,bc):
    def igafem_connectivity(self,ndof):
        # ID ARRAY
        count = 0
        for i in range(0,self.np):
            nnp_p = self.nnp[i]
            ID_p = np.zeros((nnp_p,ndof)).astype(int)
            for j in range(0,nnp_p):
                for k in range(0,ndof):
                    ID_p[j,k] = count
                    count += 1
            self.ID.append(ID_p)
        # LM ARRAY
        for i in range(0,self.np):
            LM_p = np.zeros(((self.nel[i],self.nen[i],ndof))).astype(int)
            for j in range(0,self.nel[i]):
                for k in range(0,self.nen[i]):
                    LM_p[j,k,:] = self.ID[i][self.IEN[i][j,k],:]
            self.LM.append(LM_p)
