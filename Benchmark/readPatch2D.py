import numpy as np
import os.path
import math

class readPatch2D:
    def __init__(self,file_ID):
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
            a = 0
            b = 0
            for j in range(1,m+1):
                for i in range(1,n+1):
                    INC[a][0] = (i-1)
                    INC[a][1] = (j-1)
                    if (i >= p1+1 and j >= p2+1):
                        for k in range(0,p2+1):
                            for l in range(0,p1+1):
                                c = int(k*(p1+1)+l)
                                IEN[b][c] = ((a-k*n-l))
                        b += 1
                    a += 1
            self.INC.append(INC)
            self.IEN.append(IEN)
            #
            fNo = fNo + 1
            file_exists = os.path.exists(file_ID+str(fNo))
        # Number of Patches:
        self.np = fNo-1
