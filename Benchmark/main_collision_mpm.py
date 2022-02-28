import numpy as np
from scipy.spatial import Delaunay
from matplotlib import path
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from pmesh import *
from regular2Dmesh import *
from linear2Dshapefun import *

n = 60
phi = (2*np.pi/n)*np.linspace(0,n,n+1)
pv1 = np.transpose(np.array([0.2*np.cos(phi), 0.2*np.sin(phi)]))
points,tri = pmesh(pv1,0.032,0)
set1 = (1/3.0)*np.transpose(np.array([points[tri[:,0],0]+points[tri[:,1],0]+points[tri[:,2],0],
                                      points[tri[:,0],1]+points[tri[:,1],1]+points[tri[:,2],1]]))
set1 = set1 + 0.2*np.ones((len(tri),2))
set2 = set1 + 0.6*np.ones((len(tri),2))
p_set = np.append(set1,set2,axis=0)
L1 = 1.0
L2 = 1.0
h1 = 0.05
h2 = 0.05
nodes, conn, el1, el2, h1, h2 = regular2Dmesh(L1,h1,L2,h2)
#plt.figure(figsize=(12,12))
#for i in range(0,len(conn)):
#    plt.gca().add_patch(patches.Polygon([[nodes[conn[i,0],0],nodes[conn[i,0],1]],[nodes[conn[i,1],0],nodes[conn[i,1],1]],[nodes[conn[i,2],0],nodes[conn[i,2],1]],[nodes[conn[i,3],0],nodes[conn[i,3],1]]],facecolor='w',edgecolor='k',fill=False))
#plt.scatter(p_set[:,0],p_set[:,1])
#plt.show()
##########
E = 1000.0
nu = 0.3
C=(E/((1.+nu)*(1.-2.*nu)))*np.array([[1-nu,nu,0],
                                     [nu,1-nu,0],
                                     [0,0,(1-2*nu)/2]])
rho = 1000.0
v = 0.1
TOL = 1E-12
Vp1 = np.zeros(len(set1))
Mp1 = np.zeros(len(set1))
vp1 = np.zeros((len(set1),2))
Sp1 = np.zeros((len(set1),3))
Ep1 = np.zeros((len(set1),3))
Fp1 = np.zeros((len(set1),4))
#############################
# Initial State
for i in range(0,len(set1)):
    aMatrix = np.array([[points[tri[i,0],0],points[tri[i,0],1],1.0],
                        [points[tri[i,1],0],points[tri[i,1],1],1.0],
                        [points[tri[i,2],0],points[tri[i,2],1],1.0]])
    a = np.linalg.det(aMatrix)/2.0
    Vp1[i] = a
    Mp1[i] = a*rho
    vp1[i,:] = [v,v]
    Fp1[i,:] = [1.0,0.0,0.0,1.0]
#############################
# All Particles
Vp = np.append(Vp1,Vp1,axis=0)
Mp = np.append(Mp1,Mp1,axis=0)
vp = np.append(vp1,-1.0*vp1,axis=0)
rp = p_set
Sp = np.append(Sp1,Sp1,axis=0)
Ep = np.append(Ep1,Ep1,axis=0)
Fp = np.append(Fp1,Fp1,axis=0)
Vp0 = Vp
elNp = np.zeros(len(rp))
#############################
n_mass = np.zeros(len(nodes))
n_momentum = np.zeros((len(nodes),2))
n_fi = np.zeros((len(nodes),2))
#############################
# Simulation Parameters:
deltaT = 0.005
t = 0.0
time = 3.6
frame = 100
X = np.zeros((len(rp),frame))
Y = np.zeros((len(rp),frame))
X[:,0]=rp[:,0]
Y[:,0]=rp[:,1]
Tstep = time/frame
step = int(np.floor(Tstep/deltaT))
count = 0
istep = 1
while (t<=time):
    print(t)
    if ((count)%step == 0):
        if (istep<frame):
            X[:,istep]=rp[:,0]
            Y[:,istep]=rp[:,1]
            istep+=1
    count = count + 1
    # From Particles to Nodes:
    for i in range(0,len(rp)):
        xp = rp[i,0]
        yp = rp[i,1]
        elp1 = np.ceil(xp/h1)
        elp2 = np.floor(yp/h2)
        el = int(elp2*el1+elp1-1)
        elNp[i] = el
        el_nodes = nodes[conn[el,:],:]
        xi = (2.0*xp-(el_nodes[0,0]+el_nodes[1,0]))/h1
        eta = (2.0*yp-(el_nodes[1,1]+el_nodes[2,1]))/h2
        N,dN = linear2Dshapefun(xi,eta)
        Jmatrix = np.array([[dN[0,:].dot(el_nodes[:,0]),dN[0,:].dot(el_nodes[:,1])],
                            [dN[1,:].dot(el_nodes[:,0]),dN[1,:].dot(el_nodes[:,1])]])
        J = np.linalg.det(Jmatrix)
        dNxy = np.linalg.inv(Jmatrix).dot(dN)
        #dNxy = np.linalg.solve(Jmatrix,dN)
        for j in range(0,4):
            n_mass[conn[el,j]] += N[j]*Mp[i]
            n_momentum[conn[el,j],:] += (N[j]*Mp[i])*vp[i,:]
            n_fi[conn[el,j],:] += (-1.0*Vp[i])*np.array([(Sp[i,0]*dNxy[0,j]+Sp[i,2]*dNxy[1,j]),
                                                         (Sp[i,1]*dNxy[1,j]+Sp[i,2]*dNxy[0,j])])
    n_momentum += deltaT*n_fi
    #From Nodes to Particles
    for i in range(0,len(rp)):
        xp = rp[i,0]
        yp = rp[i,1]
        el = int(elNp[i])
        el_nodes = nodes[conn[el,:],:]
        xi = (2.0*xp-(el_nodes[0,0]+el_nodes[1,0]))/h1
        eta = (2.0*yp-(el_nodes[1,1]+el_nodes[2,1]))/h2
        N,dN = linear2Dshapefun(xi,eta)
        Jmatrix = np.array([[dN[0,:].dot(el_nodes[:,0]),dN[0,:].dot(el_nodes[:,1])],
                            [dN[1,:].dot(el_nodes[:,0]),dN[1,:].dot(el_nodes[:,1])]])
        #J = np.linalg.det(Jmatrix)
        dNxy = np.linalg.inv(Jmatrix).dot(dN)
        #dNxy = np.linalg.solve(Jmatrix,dN)
        Lp = np.zeros((2,2))
        for j in range(0,4):
            vl = np.zeros(2)
            if (n_mass[conn[el,j]]>TOL):
                vp[i,:] += (deltaT*N[j]/n_mass[conn[el,j]])*(n_fi[conn[el,j],:])
                rp[i,:] += (deltaT*N[j]/n_mass[conn[el,j]])*n_momentum[conn[el,j],:]
                vl = n_momentum[conn[el,j],:]/n_mass[conn[el,j]]
            Lp[0,0] += vl[0]*dNxy[0,j]
            Lp[0,1] += vl[0]*dNxy[1,j]
            Lp[1,0] += vl[1]*dNxy[0,j]
            Lp[1,1] += vl[1]*dNxy[1,j]
        F = (np.identity(2)+deltaT*Lp).dot(np.array([[Fp[i,0],Fp[i,2]],
                                                     [Fp[i,1],Fp[i,3]]]))
        Fp[i,:] = [F[0,0],F[1,0],F[0,1],F[1,1]]
        Vp[i] = np.linalg.det(F)*Vp0[i]
        dEps = (0.5*deltaT)*(Lp+np.transpose(Lp))
        dSigma = C.dot([dEps[0,0],dEps[1,1],2*dEps[0,1]])
        Sp[i,:] += [dSigma[0], dSigma[1], dSigma[2]]
        Ep[i,:] += [dEps[0,0],dEps[1,1],dEps[0,1]]
    t += deltaT
    n_mass = np.zeros(len(nodes))
    n_momentum = np.zeros((len(nodes),2))
    n_fi = np.zeros((len(nodes),2))

plt.figure(figsize=(8,8))
for i in range(0,len(conn)):
    plt.gca().add_patch(patches.Polygon([[nodes[conn[i,0],0],nodes[conn[i,0],1]],[nodes[conn[i,1],0],nodes[conn[i,1],1]],[nodes[conn[i,2],0],nodes[conn[i,2],1]],[nodes[conn[i,3],0],nodes[conn[i,3],1]]],facecolor='w',edgecolor='k',fill=False))
plt.scatter(X[:,95],Y[:,95])
plt.show()
