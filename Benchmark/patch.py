import numpy as np
import os.path
import math

class patchFamily2D:
    def __init__(self,file_ID,ndof):
        # Knot Vectors
        self.U = []
        # Control Points
        self.CP = []
        # Orders
        self.p = []
        # Number of Shape Functions
        self.nsf = []
        # Connectivity
        self.nel = []
        self.nnp = []
        self.nen = []
        self.INC = []
        self.IEN = []
        self.ID = []
        self.LM = []
        # READ PATCH FAMILY NAMED WITH "file_ID"
        file_exists = os.path.exists(file_ID+'1')
        fNo = 1
        while (file_exists == True):
            f = open(file_ID+str(fNo), "r")
            data = f.readlines()
            # Knot Vector - V
            lenV = int(data[0])
            V_p = np.zeros(lenV)
            dataV = data[1].split(" ")
            for i in range(0,lenV):
                V_p[i] = float(dataV[i])
            # Knot Vector - U
            lenU = int(data[2])
            U_p = np.zeros(lenU)
            dataU = data[3].split(" ")
            for i in range(0,lenU):
                U_p[i] = float(dataU[i])
            self.U.append([])
            self.U[fNo-1].append(U_p)
            self.U[fNo-1].append(V_p)
            # Control Points:
            lenCP = int(data[4])
            CP_p = np.zeros((lenCP,4))
            for i in range(0,lenCP):
                dataCP = data[i+5].split(" ")
                for j in range(0,4):
                    CP_p[i][j] = float(dataCP[j])
            self.CP.append(CP_p)
            # Order:
            p1 = 0
            p2 = 0
            while True:
                if (U_p[p1] == U_p[p1+1]):
                    p1 += 1
                else:
                    break
            while True:
                if (V_p[p2] == V_p[p2+1]):
                    p2 += 1
                else:
                    break
            self.p.append(np.array([p1,p2]))
            n = len(U_p)-p1-1
            m = len(V_p)-p2-1
            self.nsf.append(np.array([n,m]))
            # Connectivity:
            nel = (n-p1)*(m-p2)
            self.nel.append(nel)
            nnp = n*m
            self.nnp.append(nnp)
            nen = (p1+1)*(p2+1)
            self.nen.append(nen)
            INC = np.zeros((n*m,2)).astype(int)
            IEN = np.zeros((nel,nen)).astype(int)
            A = 0
            B = 0
            for j in range(1,m+1):
                for i in range(1,n+1):
                    INC[A][0] = (i-1)
                    INC[A][1] = (j-1)
                    if (i >= p1+1 and j >= p2+1):
                        for k in range(0,p2+1):
                            for l in range(0,p1+1):
                                C = int(k*(p1+1)+l)
                                IEN[B][C] = ((A-k*n-l))
                        B += 1
                    A += 1
            self.INC.append(INC)
            self.IEN.append(IEN)
            #
            fNo = fNo + 1
            file_exists = os.path.exists(file_ID+str(fNo))
        # Number of Patches:
        self.np = fNo-1
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
    ############
    # SPAN INDEX
    def findSpan(self,pN,pD,u):
        if (u == self.U[pN][pD][self.nsf[pN][pD]+1]):
            return n
        low = self.p[pN][pD]
        high = self.nsf[pN][pD]+1
        mid = int((low+high)/2)
        while (u<self.U[pN][pD][mid] or u>self.U[pN][pD][mid+1]):
            if (u<self.U[pN][pD][mid]):
                high = mid
            else:
                low = mid
            mid = int(np.floor((low+high)/2))
        return mid
    ############
    # DERIVATIVES
    def dersBasisFuns(self,pN,pD,u,n):
        i = self.findSpan(pN,pD,u)
        left = np.zeros(self.p[pN][pD]+1)
        right = np.zeros(self.p[pN][pD]+1)
        ndu = np.zeros((self.p[pN][pD]+1,self.p[pN][pD]+1))
        a = np.zeros((2,self.p[pN][pD]+1))
        ders = np.zeros((n+1,self.p[pN][pD]+1))
        #
        ndu[0][0] = 1.0
        for j in range(1,self.p[pN][pD]+1):
            left[j] = u-self.U[pN][pD][i+1-j]
            right[j] = self.U[pN][pD][i+j]-u
            saved = 0.0
            for r in range(0,j):
                ndu[j][r] = right[r+1]+left[j-r]
                temp = ndu[r][j-1]/ndu[j][r]
                ndu[r][j] = saved+right[r+1]*temp
                saved = left[j-r]*temp
            ndu[j][j] = saved
        for j in range(0,self.p[pN][pD]+1):
            ders[0][j] = ndu[j][self.p[pN][pD]]
        for r in range(0,self.p[pN][pD]+1):
            s1 = 0
            s2 = 1
            a[0][0] = 1.0
            for k in range(1,n+1):
                d =0.0
                rk = r-k
                pk = self.p[pN][pD]-k
                if (r >= k):
                    a[s2][0] = a[s1][0]/ndu[pk+1][rk]
                    d = a[s2][0]*ndu[rk][pk]
                if (rk >= -1):
                    j1 = 1
                else:
                    j1 = -rk
                if (r-1 <= pk):
                    j2 = k-1
                else:
                    j2 = self.p[pN][pD]-r
                for j in range(j1,j2+1):
                    a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j]
                    d += a[s2][j]*ndu[rk+j][pk]
                if (r <= pk):
                    a[s2][k] = -a[s1][k-1]/ndu[pk+1][r]
                    d += a[s2][k]*ndu[r][pk]
                ders[k][r] = d
                j = s1
                s1 = s2
                s2 = j
        r = self.p[pN][pD]
        for k in range(1,n+1):
            for j in range(0,self.p[pN][pD]+1):
                ders[k][j] *= r
            r *= (self.p[pN][pD]-k)
        return ders

    def getCP(self,pN,p_el):
        cp = np.zeros(((3,self.p[pN][0]+1,self.p[pN][1]+1)))
        wMat = np.zeros((self.p[pN][0]+1,self.p[pN][1]+1))
        iu = self.INC[pN][self.IEN[pN][p_el][0],0]
        iv = self.INC[pN][self.IEN[pN][p_el][0],1]
        ny = iv-self.p[pN][1]
        k = 0
        for j in range(0,self.p[pN][1]+1):
            nx = iu-self.p[pN][0]
            for i in range(0,self.p[pN][0]+1):
                cp[0][i,j] = self.CP[pN][self.nsf[pN][0]*ny+nx][0]
                cp[1][i,j] = self.CP[pN][self.nsf[pN][0]*ny+nx][1]
                cp[2][i,j] = self.CP[pN][self.nsf[pN][0]*ny+nx][2]
                wMat[i,j] = self.CP[pN][self.nsf[pN][0]*ny+nx][3]
                nx += 1
                k += 1
            ny += 1
        return cp, wMat

    def ratBasisFuns2D(self,pN,u,v,wMat,du,dv):
        N = self.dersBasisFuns(pN,0,u,du)
        M = self.dersBasisFuns(pN,1,v,dv)
        R = np.zeros((((self.p[pN][0]+1,self.p[pN][1]+1,du+1,dv+1))))
        wders = np.zeros((du+1,dv+1))
        for k in range(0,du+1):
            for l in range(0,dv+1):
                wders[k][l] = N[k][:].dot(wMat.dot(M[l][:]))
                temp1 = np.tensordot(N[k][:],(M[l][:]),0)*wMat
                for j in range(1,l+1):
                    temp1 -= (self.nchoosek(l,j)*wders[0][j])*R[:,:,k,l-j]
                for i in range(1,k+1):
                    temp1 -= (self.nchoosek(k,i)*wders[i][0])*R[:,:,k-i,l]
                    temp2 = np.zeros((self.p[pN][0]+1,self.p[pN][1]+1))
                    for j in range (1,l+1):
                        temp2 += (self.nchoosek(l,j)*wders[i][j])*R[:,:,k-i,l-j]
                    temp1 -= self.nchoosek(k,i)*temp2
                R[:,:,k,l] = temp1/wders[0][0]
        return R

    def nchoosek(self,n,r):
        return math.factorial(n)/math.factorial(r)/math.factorial(n-r)

    def surfaceDerivs(self,R,cp,du,dv):
        S = np.zeros(((3,du+1,dv+1)))
        for i in range(0,du+1):
            for j in range(0,dv+1):
                S[0][i,j] = np.tensordot(cp[0],R[:,:,i,j])
                S[1][i,j] = np.tensordot(cp[1],R[:,:,i,j])
                S[2][i,j] = np.tensordot(cp[2],R[:,:,i,j])
        return S
