import numpy as np
from readPatch2D import *
from bc import *

class igafem(bc,readPatch2D):
    def __init__(self,file_ID):
        super().__init__(file_ID)
        self.dof = []
        self.ldof = []
    # IGAFEM CONNECTIVITY:
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
            self.dof.append(count*ndof)
            self.ldof.append((self.nen[0]*ndof)**2)
        # LM
        for i in range(0,self.np):
            LM_p = np.zeros(((self.nel[i],self.nen[i],ndof))).astype(int)
            for j in range(0,self.nel[i]):
                for k in range(0,self.nen[i]):
                    LM_p[j,k,:] = self.ID[i][self.IEN[i][j,k],:]
                    self.LM.append(LM_p)
        # FEM PROCEDURE
        def global_stiffness(self,ngp):

            x_gp, w_gp = self.quad2D(ngp)
            K = PETSc.Mat().create(comm=comm)
            K.setSizes((self.dof[0],self.dof[0]))
            K.setFromOptions()
            K.setPreallocationNNZ(self.ldof[0]*self.nen[0])
            for i in range(0,self.nel[0]):
                el_CP = self.getCP(0,i)
                el_LM = self.LM[0][i,:,:]
                k_local = self.local_stiffness()
            return K

        def quad1D(self,ngp):
            if ngp == 1:
                xgp = [0.0]
                wgp = [2.0]
            elif ngp == 2:
                xgp = [-1./sqrt(3.), 1./sqrt(3.)];
                wgp = [1., 1.];
            elif ngp ==3:
                xgp = [-sqrt(3./5), 0., sqrt(3./5)];
                wgp = [5./9, 8./9, 5./9];

            return xgp, wgp

        def quad2D(self,ngp):
            x_ksi, w_ksi = self.quad1D(ngp)
            x_eta, w_eta = self.quad1D(ngp)
            x_gp = np.zeros((len(x_ksi)*len(x_eta)),2)
            w_gp = np.zeros(len(x_ksi)*len(x_eta))
            for i in range(0,len(x_eta)):
                k = i*len(x_ksi)
                for j in range(0,len(x_ksi)):
                    x_gp[j+k,:] = [x_ksi[j], x_eta[i]]
                    w_gp[j+k] = w_ksi[j]*w_eta[i]
            return x_gp, w_gp

        def local_stiffness(self,iu,iv,el_CP,x_gp,w_gp):
