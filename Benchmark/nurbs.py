import numpy as np
import math

class nurbs:
    ############
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
    ###############
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
    ################
    def surfaceDerivs(self,R,cp,du,dv):
        S = np.zeros(((3,du+1,dv+1)))
        for i in range(0,du+1):
            for j in range(0,dv+1):
                S[0][i,j] = np.tensordot(cp[0],R[:,:,i,j])
                S[1][i,j] = np.tensordot(cp[1],R[:,:,i,j])
                S[2][i,j] = np.tensordot(cp[2],R[:,:,i,j])
        return S
    #################
    # Control Points:
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
    #################
    def nchoosek(self,n,r):
        return math.factorial(n)/math.factorial(r)/math.factorial(n-r)
