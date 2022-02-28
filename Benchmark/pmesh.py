import numpy as np
from scipy.spatial import Delaunay
from matplotlib import path
def pmesh(pv,hmax,nref):
    edge = len(pv)-1
    for i in range(0,edge):
        d = np.sqrt((pv[i+1,0]-pv[i,0])**2+(pv[i+1,1]-pv[i,1])**2)
        if hmax<d:
            n = int(np.ceil(d/hmax))
            nx = np.linspace(pv[i,0],pv[i+1,0],n+1)
            ny = np.linspace(pv[i,1],pv[i+1,1],n+1)
            for j in range(0,len(nx[1:-1])):
                pv = np.append(pv,[[nx[j+1],ny[j+1]]],axis=0)
    poly = pv[0:edge+1,:]
    pv = np.unique(pv,axis=0)
    a = hmax**2
    p = 0
    while (hmax**2/2.0)<a:
        tri = Delaunay(pv)
        conn = tri.simplices
        A = np.zeros(len(conn))
        for i in range(0,len(conn)):
            A[i] = np.abs(0.5*(pv[conn[i,0],0]*(pv[conn[i,1],1]-pv[conn[i,2],1])+pv[conn[i,1],0]*(pv[conn[i,2],1]-pv[conn[i,0],1])+pv[conn[i,2],0]*(pv[conn[i,0],1]-pv[conn[i,1],1])))
        C = np.zeros((len(conn),2))
        D = np.zeros((len(conn),2))
        for j in range(0,len(conn)):
            K=np.array([[2*(pv[conn[j,0],0]-pv[conn[j,1],0]), 2*(pv[conn[j,0],1]-pv[conn[j,1],1])] , [2*(pv[conn[j,1],0]-pv[conn[j,2],0]), 2*(pv[conn[j,1],1]-pv[conn[j,2],1])]])
            F=np.array([((pv[conn[j,0],0])**2-(pv[conn[j,1],0])**2)+((pv[conn[j,0],1])**2-(pv[conn[j,1],1])**2) , ((pv[conn[j,1],0])**2-(pv[conn[j,2],0])**2)+((pv[conn[j,1],1])**2-(pv[conn[j,2],1])**2)])
            c = np.linalg.inv(K).dot(F)
            C[j,:] = [c[0], c[1]]
            D[j,:] = [np.mean([(pv[conn[j,0],0]),(pv[conn[j,1],0]),(pv[conn[j,2],0])]) , np.mean([(pv[conn[j,0],1]),(pv[conn[j,1],1]),(pv[conn[j,2],1])])]
        pp = path.Path(poly)
        in_check = pp.contains_points(D)
        delete = np.where(in_check == False)
        A = np.delete(A,delete)
        C = np.delete(C,delete,axis=0)
        conn = np.delete(conn,delete,axis=0)
        if max(A)<(hmax**2/2):
            break
        else:
            ind = np.argmax(A)
            pv = np.append(pv,[[C[ind,0],C[ind,1]]],axis=0)
            pv = np.unique(pv,axis=0)
        p = p+1

    return pv, conn
