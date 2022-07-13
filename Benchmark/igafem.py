import numpy as np
from readPatch2D import *
from bc import *
from nurbs import *
from model import *

class igafem(readPatch2D,bc,nurbs,model):
    def __init__(self,file_ID):
        super().__init__(file_ID)
        self.dof = []
        self.ID = []
        self.LM = []
        self.bcZero_no = []
        self.bcZero = []
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
            self.dof.append(count)
        # LM
        for i in range(0,self.np):
            LM_p = np.zeros(((self.nel[i],self.nen[i],ndof))).astype(int)
            for j in range(0,self.nel[i]):
                for k in range(0,self.nen[i]):
                    LM_p[j,k,:] = self.ID[i][self.IEN[i][j,k],:]
                    self.LM.append(LM_p)
    # FEM PROCEDURE - STIFFNESS
    #def global_stiffness(self,ngp):
    #    x_gp, w_gp = self.quad2D(ngp)
    #    K = PETSc.Mat().create(comm=None)
    #    K.setSizes((self.dof[0],self.dof[0]))
    #    K.setFromOptions()
    #    K.setPreallocationNNZ((75,75))
    #    Istart,Iend = K.getOwnershipRange()
    #    for i in range(0,self.nel[0]):
    #        print(i)
    #        iu = self.INC[0][self.IEN[0][i,0],0]
    #        iv = self.INC[0][self.IEN[0][i,0],1]
    #        el_CP, el_w = self.getCP(0,i)
    #        el_LM = np.reshape(np.flipud(self.LM[0][i,:,:]),3*self.nen[0])
    #        K_loc = self.local_stiffness(iu,iv,el_CP,el_w,x_gp,w_gp)
    #        for j in range(0,3*self.nen[0]):
    #            if el_LM[j] != -1:
    #                for k in range(0,3*self.nen[0]):
    #                    if el_LM[k] != -1:
    #                        K.setValues(el_LM[j],el_LM[k],K_loc[j,k],addv=True)
    #    K.assemblyBegin()
    #    K.assemblyEnd()
    #    return K
    # FEM PROCEDURE - MASS
    def global_mass(self,ngp):
        x_gp, w_gp = self.quad2D(ngp)
        M = PETSc.Mat().create(comm=None)
        M.setSizes((self.dof[0],self.dof[0]))
        M.setFromOptions()
        #K.setPreallocationNNZ((self.ldof[0]*self.nen[0],self.ldof[0]*self.nen[0]))
        M.setPreallocationNNZ((75,75))
        Istart,Iend = M.getOwnershipRange()
        for i in range(0,self.nel[0]):
            print(i)
            iu = self.INC[0][self.IEN[0][i,0],0]
            iv = self.INC[0][self.IEN[0][i,0],1]
            el_CP, el_w = self.getCP(0,i)
            el_LM = np.reshape(np.flipud(self.LM[0][i,:,:]),3*self.nen[0])
            M_loc = self.local_mass(iu,iv,el_CP,el_w,x_gp,w_gp)
            for j in range(0,3*self.nen[0]):
                if el_LM[j] != -1:
                    for k in range(0,3*self.nen[0]):
                        if el_LM[k] != -1:
                            M.setValues(el_LM[j],el_LM[k],M_loc[j,k],addv=True)
        M.assemblyBegin()
        M.assemblyEnd()
        return M
    # GAUSSIAN QUADRATURE 1-D
    def quad1D(self,ngp):
        if ngp == 1:
            xgp = [0.0]
            wgp = [2.0]
        elif ngp == 2:
            xgp = [-1./np.sqrt(3.), 1./np.sqrt(3.)];
            wgp = [1., 1.];
        elif ngp == 3:
            xgp = [-np.sqrt(3./5), 0., np.sqrt(3./5)];
            wgp = [5./9, 8./9, 5./9];

        return xgp, wgp
    # GAUSSIAN QUADRATURE 2-D
    def quad2D(self,ngp):
        x_ksi, w_ksi = self.quad1D(ngp)
        x_eta, w_eta = self.quad1D(ngp)
        x_gp = np.zeros(((len(x_ksi)*len(x_eta)),2))
        w_gp = np.zeros(len(x_ksi)*len(x_eta))
        for i in range(0,len(x_eta)):
            k = i*len(x_ksi)
            for j in range(0,len(x_ksi)):
                x_gp[j+k,:] = [x_ksi[j], x_eta[i]]
                w_gp[j+k] = w_ksi[j]*w_eta[i]
        return x_gp, w_gp
