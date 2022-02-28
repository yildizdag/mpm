import numpy as np
def linear2Dshapefun(xi,eta):
    N = np.zeros(4)
    N[0] = 0.25*(1.0-xi)*(1.0-eta)
    N[1] = 0.25*(1.0+xi)*(1.0-eta)
    N[2] = 0.25*(1.0+xi)*(1.0+eta)
    N[3] = 0.25*(1.0-xi)*(1.0+eta)
    dN = np.zeros((2,4))
    dN[0,0] = -0.25*(1.0-eta)
    dN[0,1] =  0.25*(1.0-eta)
    dN[0,2] =  0.25*(1.0+eta)
    dN[0,3] = -0.25*(1.0+eta)
    dN[1,0] = -0.25*(1.0-xi)
    dN[1,1] = -0.25*(1.0+xi)
    dN[1,2] =  0.25*(1.0+xi)
    dN[1,3] =  0.25*(1.0-xi)
    return N,dN
