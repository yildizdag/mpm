import numpy as np

def regular2Dmesh(L1,h1,L2,h2,p):
    if p == 1:
        el1 = int(np.ceil(L1/h1))
        el2 = int(np.ceil(L2/h2))
        h1 = L1/el1
        h2 = L2/el2
        size1 = np.linspace(0,L1,el1+1)
        size2 = np.linspace(0,L2,el2+1)
        nel = el1*el2
        X1,Y1 = np.meshgrid(size1,size2)
        nodes = np.zeros((len(size1)*len(size2),2))
        k = 0
        for i in range(0,len(size2)):
            for j in range(0,len(size1)):
                nodes[k,:] = X1[i,j],Y1[i,j]
                k += 1
        conn = np.zeros((nel,4)).astype(int)
        for j in range(0,el2):
            for i in range(0,el1):
                conn[i+el1*j,:] = i+(el1+1)*j, i+(el1+1)*j+1, i+(el1+1)*(j+1)+1, i+(el1+1)*(j+1)
    elif p != 1:
        print('p must be equal to 1!')
    return nodes,conn,el1,el2,h1,h2
