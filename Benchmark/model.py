import numpy as np

class model:
    def local_stiffness(self,iu,iv,el_CP,el_w,x_gp,w_gp):
        K_loc = np.zeros((3*self.nen[0],3*self.nen[0]))
        E = 2.1E11
        nu = 0.3
        t = 0.00484
        cu1 = self.U[0][0][iu] + self.U[0][0][iu+1]
        cu2 = self.U[0][0][iu+1] - self.U[0][0][iu]
        cv1 = self.U[0][1][iv] + self.U[0][1][iv+1]
        cv2 = self.U[0][1][iv+1] - self.U[0][1][iv]
        for i in range(0,len(w_gp)):
            u = 0.5*(cu1+cu2*x_gp[i,0])
            v = 0.5*(cv1+cv2*x_gp[i,1])
            R = self.ratBasisFuns2D(0,u,v,el_w,2,2)
            S = self.surfaceDerivs(R,el_CP,2,2)
            A1 = S[:,1,0]
            A2 = S[:,0,1]
            J1 = np.linalg.norm(np.cross(A1,A2))
            A3 = np.cross(A1,A2)/J1
            A1d1 = S[:,2,0]
            A2d2 = S[:,0,2]
            A1d2 = S[:,1,1]
            A2d1 = A1d2
            A3d1 = (1/J1)*(np.cross(A1d1,A2)+np.cross(A1,A2d1)-(np.dot(A1d1,np.cross(A2,A3))+np.dot(A2d1,np.cross(A3,A1)))*A3)
            A3d2 = (1/J1)*(np.cross(A1d2,A2)+np.cross(A1,A2d2)-(np.dot(A1d2,np.cross(A2,A3))+np.dot(A2d2,np.cross(A3,A1)))*A3)
            A11 = np.dot(A1,A1)
            A12 = np.dot(A1,A2)
            A22 = np.dot(A2,A2)
            Ac = np.linalg.inv(np.array([[A11, A12], [A12, A22]]))
            C = np.array([[Ac[0,0]**2, nu*Ac[0,0]*Ac[1,1]+(1-nu)*Ac[0,1]**2, Ac[0,0]*Ac[0,1]],
                         [nu*Ac[0,0]*Ac[1,1]+(1-nu)*Ac[0,1]**2, Ac[1,1]**2, Ac[1,1]*Ac[0,1]],
                         [Ac[0,0]*Ac[0,1], Ac[1,1]*Ac[0,1], 0.5*((1-nu)*Ac[0,0]*Ac[1,1]+(1+nu)*Ac[0,1]**2)]])
            J2 = np.array([[0.5*(self.U[0][0][iu+1]-self.U[0][0][iu]),0], [0, 0.5*(self.U[0][1][iv+1]-self.U[0][1][iv])]])
            J = J1*np.linalg.det(J2)
            # Second Derivatives
            dR2B = np.zeros((3,self.nen[0]))
            dR2B[0,:] = np.reshape(R[:,:,2,0],self.nen[0],order='F')
            dR2B[1,:] = np.reshape(R[:,:,0,2],self.nen[0],order='F')
            dR2B[2,:] = np.reshape(R[:,:,1,1],self.nen[0],order='F')
            # First Derivatives
            dR1B = np.zeros((2,self.nen[0]))
            dR1B[0,:] = np.reshape(R[:,:,1,0],self.nen[0],order='F')
            dR1B[1,:] = np.reshape(R[:,:,0,1],self.nen[0],order='F')
            # B1
            B1 = np.zeros((3,3*self.nen[0]))
            B1[0,0::3] = dR1B[0,:]*A1[0]
            B1[0,1::3] = dR1B[0,:]*A1[1]
            B1[0,2::3] = dR1B[0,:]*A1[2]
            B1[1,0::3] = dR1B[1,:]*A2[0]
            B1[1,1::3] = dR1B[1,:]*A2[1]
            B1[1,2::3] = dR1B[1,:]*A2[2]
            B1[2,0::3] = dR1B[0,:]*A2[0]+dR1B[1,:]*A1[0]
            B1[2,1::3] = dR1B[0,:]*A2[1]+dR1B[1,:]*A1[1]
            B1[2,2::3] = dR1B[0,:]*A2[2]+dR1B[1,:]*A1[2]
            #B2
            B2 = np.zeros((3,3*self.nen[0]))
            phi1x = (1./J1)*(A3[0]*dR1B[1,:])
            phi1y = (1./J1)*(A3[1]*dR1B[1,:])
            phi1z = (1./J1)*(A3[2]*dR1B[1,:])

            phi2x = (1./J1)*(-A3[0]*dR1B[0,:])
            phi2y = (1./J1)*(-A3[1]*dR1B[0,:])
            phi2z = (1./J1)*(-A3[2]*dR1B[0,:])

            phi1d1x = (1./J1)*(A3[0]*dR2B[2,:]+A3d1[0]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[0]*dR1B[1,:]))
            phi1d1y = (1./J1)*(A3[1]*dR2B[2,:]+A3d1[1]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[1]*dR1B[1,:]))
            phi1d1z = (1./J1)*(A3[2]*dR2B[2,:]+A3d1[2]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[2]*dR1B[1,:]))

            phi2d1x = (1./J1)*(-A3[0]*dR2B[0,:]-A3d1[0]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[0]*dR1B[0,:]))
            phi2d1y = (1./J1)*(-A3[1]*dR2B[0,:]-A3d1[1]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[1]*dR1B[0,:]))
            phi2d1z = (1./J1)*(-A3[2]*dR2B[0,:]-A3d1[2]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d1)+np.dot(np.cross(A3,A1),A2d1)))*(A3[2]*dR1B[0,:]))

            phi1d2x = (1./J1)*(A3[0]*dR2B[1,:]+A3d2[0]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[0]*dR1B[1,:]))
            phi1d2y = (1./J1)*(A3[1]*dR2B[1,:]+A3d2[1]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[1]*dR1B[1,:]))
            phi1d2z = (1./J1)*(A3[2]*dR2B[1,:]+A3d2[2]*dR1B[1,:] - ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[2]*dR1B[1,:]))

            phi2d2x = (1./J1)*(-A3[0]*dR2B[2,:]-A3d2[0]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[0]*dR1B[0,:]))
            phi2d2y = (1./J1)*(-A3[1]*dR2B[2,:]-A3d2[1]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[1]*dR1B[0,:]))
            phi2d2z = (1./J1)*(-A3[2]*dR2B[2,:]-A3d2[2]*dR1B[0,:] + ((1./J1)*(np.dot(np.cross(A2,A3),A1d2)+np.dot(np.cross(A3,A1),A2d2)))*(A3[2]*dR1B[0,:]))

            B2[0,0::3] = A3d1[0]*dR1B[0,:] + np.dot(np.cross(A1d1,A3),A1)*phi1x + np.dot(np.cross(A2,A3),A1)*phi2d1x + np.dot(np.cross(A2d1,A3),A1)*phi2x
            B2[0,1::3] = A3d1[1]*dR1B[0,:] + np.dot(np.cross(A1d1,A3),A1)*phi1y + np.dot(np.cross(A2,A3),A1)*phi2d1y + np.dot(np.cross(A2d1,A3),A1)*phi2y
            B2[0,2::3] = A3d1[2]*dR1B[0,:] + np.dot(np.cross(A1d1,A3),A1)*phi1z + np.dot(np.cross(A2,A3),A1)*phi2d1z + np.dot(np.cross(A2d1,A3),A1)*phi2z

            B2[1,0::3] = A3d2[0]*dR1B[1,:] + np.dot(np.cross(A1,A3),A2)*phi1d2x + np.dot(np.cross(A1d2,A3),A2)*phi1x + np.dot(np.cross(A2d2,A3),A2)*phi2x
            B2[1,1::3] = A3d2[1]*dR1B[1,:] + np.dot(np.cross(A1,A3),A2)*phi1d2y + np.dot(np.cross(A1d2,A3),A2)*phi1y + np.dot(np.cross(A2d2,A3),A2)*phi2y
            B2[1,2::3] = A3d2[2]*dR1B[1,:] + np.dot(np.cross(A1,A3),A2)*phi1d2z + np.dot(np.cross(A1d2,A3),A2)*phi1z + np.dot(np.cross(A2d2,A3),A2)*phi2z

            B2[2,0::3] = A3d2[0]*dR1B[0,:] + np.dot(np.cross(A1,A3),A2)*phi1d1x + np.dot(np.cross(A1d1,A3),A2)*phi1x + np.dot(np.cross(A2d1,A3),A2)*phi2x + A3d1[0]*dR1B[1,:] + np.dot(np.cross(A1d2,A3),A1)*phi1x + np.dot(np.cross(A2,A3),A1)*phi2d2x + np.dot(np.cross(A2d2,A3),A1)*phi2x
            B2[2,1::3] = A3d2[1]*dR1B[0,:] + np.dot(np.cross(A1,A3),A2)*phi1d1y + np.dot(np.cross(A1d1,A3),A2)*phi1y + np.dot(np.cross(A2d1,A3),A2)*phi2y + A3d1[1]*dR1B[1,:] + np.dot(np.cross(A1d2,A3),A1)*phi1y + np.dot(np.cross(A2,A3),A1)*phi2d2y + np.dot(np.cross(A2d2,A3),A1)*phi2y
            B2[2,2::3] = A3d2[2]*dR1B[0,:] + np.dot(np.cross(A1,A3),A2)*phi1d1z + np.dot(np.cross(A1d1,A3),A2)*phi1z + np.dot(np.cross(A2d1,A3),A2)*phi2z + A3d1[2]*dR1B[1,:] + np.dot(np.cross(A1d2,A3),A1)*phi1z + np.dot(np.cross(A2,A3),A1)*phi2d2z + np.dot(np.cross(A2d2,A3),A1)*phi2z

            K_loc += (w_gp[i]*J)*(np.matmul(B1.transpose(),np.matmul(((E*t)/(1-nu**2))*C,B1)) +  np.matmul(B2.transpose(),np.matmul((((E*t**3)/(12*(1-nu**2)))*C),B2)))

        return K_loc

    def local_mass(self,iu,iv,el_CP,el_w,x_gp,w_gp):
        M_loc = np.zeros((3*self.nen[0],3*self.nen[0]))
        rho = 7800
        t = 0.00484
        for i in range(0,len(w_gp)):
            cu1 = self.U[0][0][iu] + self.U[0][0][iu+1]
            cu2 = self.U[0][0][iu+1] - self.U[0][0][iu]
            u = 0.5*(cu1+cu2*x_gp[i,0])
            cv1 = self.U[0][1][iv] + self.U[0][1][iv+1]
            cv2 = self.U[0][1][iv+1] - self.U[0][1][iv]
            v = 0.5*(cv1+cv2*x_gp[i,1])
            R = self.ratBasisFuns2D(0,u,v,el_w,1,1)
            S = self.surfaceDerivs(R,el_CP,1,1)
            A1 = S[:,1,0]
            A2 = S[:,0,1]
            J1 = np.linalg.norm(np.cross(A1,A2))
            A3 = np.cross(A1,A2)/J1
            J2 = np.array([[0.5*(self.U[0][0][iu+1]-self.U[0][0][iu]),0], [0, 0.5*(self.U[0][1][iv+1]-self.U[0][1][iv])]])
            J = J1*np.linalg.det(J2)
            # Second Derivatives
            N = np.zeros((3,3*self.nen[0]))
            N[0,0::3] = np.reshape(R[:,:,0,0],self.nen[0],order='F')
            N[1,1::3] = np.reshape(R[:,:,0,0],self.nen[0],order='F')
            N[2,2::3] = np.reshape(R[:,:,0,0],self.nen[0],order='F')
            # First Derivatives
            C = np.array([[rho,0,0],[0,rho,0],[0,0,rho]])

            M_loc += (w_gp[i]*J)*(N.transpose().dot((t*C).dot(N)))

        return M_loc
