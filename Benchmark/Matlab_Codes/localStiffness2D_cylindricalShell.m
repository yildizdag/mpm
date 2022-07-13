function K_loc = localStiffness2D_cylindricalShell(iu,iv,Nurbs2D,kk)
local_dof = Nurbs2D.local_dof_dry;
K_loc = zeros(local_dof*Nurbs2D.nen{kk},local_dof*Nurbs2D.nen{kk});
nu = Nurbs2D.nu;
E = Nurbs2D.E;
t = Nurbs2D.t;
CP = Nurbs2D.cPoints{kk}(:, iu-Nurbs2D.order{kk}(1)+1:iu, iv-Nurbs2D.order{kk}(2)+1:iv);
     
for i = 1:Nurbs2D.ngp{kk} 
    cu1 = (Nurbs2D.knots.U{kk}(iu) + Nurbs2D.knots.U{kk}(iu+1));    
    cu2 = (Nurbs2D.knots.U{kk}(iu+1) - Nurbs2D.knots.U{kk}(iu));
    u = 0.5*(cu1+cu2*Nurbs2D.xgp{kk}(i,1));
    cv1 = (Nurbs2D.knots.V{kk}(iv) + Nurbs2D.knots.V{kk}(iv+1));
    cv2 = (Nurbs2D.knots.V{kk}(iv+1) - Nurbs2D.knots.V{kk}(iv));
    v = 0.5*(cv1+cv2*Nurbs2D.xgp{kk}(i,2));
    dNu = dersbasisfuns(iu,u,Nurbs2D.order{kk}(1)-1,2,Nurbs2D.knots.U{kk});
    dNv = dersbasisfuns(iv,v,Nurbs2D.order{kk}(2)-1,2,Nurbs2D.knots.V{kk});
    [dR,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{kk}(1),Nurbs2D.order{kk}(2),CP,2,2);
    % Covariant Basis Vectors:
    A1 = dS(:,2,1);
    A2 = dS(:,1,2);
    J1 = norm(cross(A1,A2));
    A3 = cross(A1,A2)/J1;
    % Derivatives of Basis Vectors:
    A1d1 = dS(:,3,1);
    A2d2 = dS(:,1,3);
    A1d2 = dS(:,2,2);
    A2d1 = A1d2;
	A3d1 = (1/J1).*(cross(A1d1,A2)+cross(A1,A2d1)-(dot(A1d1,cross(A2,A3))+dot(A2d1,cross(A3,A1))).*A3);
	A3d2 = (1/J1).*(cross(A1d2,A2)+cross(A1,A2d2)-(dot(A1d2,cross(A2,A3))+dot(A2d2,cross(A3,A1))).*A3);
    % First Fundemental Form:
    A11 = dot(A1,A1); A12 = dot(A1,A2); A22 = dot(A2,A2);
    % Contravariant Components:
    Ac = [A11 A12; A12 A22]\[1 0; 0 1];
    C = [Ac(1,1)^2, nu*Ac(1,1)*Ac(2,2)+(1-nu)*Ac(1,2)^2, Ac(1,1)*Ac(1,2);
            nu*Ac(1,1)*Ac(2,2)+(1-nu)*Ac(1,2)^2, Ac(2,2)^2, Ac(2,2)*Ac(1,2);
            Ac(1,1)*Ac(1,2), Ac(2,2)*Ac(1,2), 0.5*((1-nu)*Ac(1,1)*Ac(2,2)+(1+nu)*Ac(1,2)^2)];
    % From Parametric to Parent:
    J2 = diag([Nurbs2D.knots.U{kk}(iu+1) - Nurbs2D.knots.U{kk}(iu),...
                        Nurbs2D.knots.V{kk}(iv+1) - Nurbs2D.knots.V{kk}(iv)])/2;
    % Jacobian:
    J = J1*det(J2);
    % Second Derivatives:
    dR2B = zeros(3,Nurbs2D.nen{kk});
    dR2B(1,:) = reshape(dR(:,:,3,1),1,Nurbs2D.nen{kk});
    dR2B(2,:) = reshape(dR(:,:,1,3),1,Nurbs2D.nen{kk});
    dR2B(3,:) = reshape(dR(:,:,2,2),1,Nurbs2D.nen{kk});
    % First Derivatives:
    dR1B = zeros(2,Nurbs2D.nen{kk});
    dR1B(1,:) = reshape(dR(:,:,2,1),1,Nurbs2D.nen{kk});
    dR1B(2,:) = reshape(dR(:,:,1,2),1,Nurbs2D.nen{kk});

    B1 = zeros(3,local_dof*Nurbs2D.nen{kk});        
    B1(1,1:3:end) = dR1B(1,:).*A1(1); B1(1,2:3:end) = dR1B(1,:).*A1(2); B1(1,3:3:end) =  dR1B(1,:).*A1(3);
    B1(2,1:3:end) = dR1B(2,:).*A2(1); B1(2,2:3:end) = dR1B(2,:).*A2(2); B1(2,3:3:end) =  dR1B(2,:).*A2(3);
    B1(3,1:3:end) = dR1B(1,:).*A2(1)+dR1B(2,:).*A1(1);
    B1(3,2:3:end) = dR1B(1,:).*A2(2)+dR1B(2,:).*A1(2);
    B1(3,3:3:end) = dR1B(1,:).*A2(3)+dR1B(2,:).*A1(3);
        
    %Build Bending Matrix:
    B2 = zeros(3,local_dof*Nurbs2D.nen{kk}); 
	
    phi1x = (1/J1).*(A3(1).*dR1B(2,:));
	phi1y = (1/J1).*(A3(2).*dR1B(2,:));
	phi1z = (1/J1).*(A3(3).*dR1B(2,:));
	
    phi2x = (1/J1).*(-A3(1).*dR1B(1,:));
	phi2y = (1/J1).*(-A3(2).*dR1B(1,:));
	phi2z = (1/J1).*(-A3(3).*dR1B(1,:));
        
	phi1d1x = (1/J1).*(A3(1).*dR2B(3,:)+A3d1(1).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(1).*dR1B(2,:)));
	phi1d1y = (1/J1).*(A3(2).*dR2B(3,:)+A3d1(2).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(2).*dR1B(2,:)));
	phi1d1z = (1/J1).*(A3(3).*dR2B(3,:)+A3d1(3).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(3).*dR1B(2,:)));
        
	phi2d1x = (1/J1).*(-A3(1).*dR2B(1,:)-A3d1(1).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(1).*dR1B(1,:)));
	phi2d1y = (1/J1).*(-A3(2).*dR2B(1,:)-A3d1(2).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(2).*dR1B(1,:)));
	phi2d1z = (1/J1).*(-A3(3).*dR2B(1,:)-A3d1(3).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d1)+dot(cross(A3,A1),A2d1))).*(A3(3).*dR1B(1,:)));
        
	phi1d2x = (1/J1).*(A3(1).*dR2B(2,:)+A3d2(1).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(1).*dR1B(2,:)));
	phi1d2y = (1/J1).*(A3(2).*dR2B(2,:)+A3d2(2).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(2).*dR1B(2,:)));
	phi1d2z = (1/J1).*(A3(3).*dR2B(2,:)+A3d2(3).*dR1B(2,:) - ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(3).*dR1B(2,:)));
        
	phi2d2x = (1/J1).*(-A3(1).*dR2B(3,:)-A3d2(1).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(1).*dR1B(1,:)));
	phi2d2y = (1/J1).*(-A3(2).*dR2B(3,:)-A3d2(2).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(2).*dR1B(1,:)));
	phi2d2z = (1/J1).*(-A3(3).*dR2B(3,:)-A3d2(3).*dR1B(1,:) + ((1/J1)*(dot(cross(A2,A3),A1d2)+dot(cross(A3,A1),A2d2))).*(A3(3).*dR1B(1,:)));
        
    %SIGMA_11
	B2(1,1:3:end) = A3d1(1).*dR1B(1,:) + dot(cross(A1d1,A3),A1).*phi1x + dot(cross(A2,A3),A1).*phi2d1x + dot(cross(A2d1,A3),A1).*phi2x;
	B2(1,2:3:end) = A3d1(2).*dR1B(1,:) + dot(cross(A1d1,A3),A1).*phi1y + dot(cross(A2,A3),A1).*phi2d1y + dot(cross(A2d1,A3),A1).*phi2y; 
	B2(1,3:3:end) = A3d1(3).*dR1B(1,:) + dot(cross(A1d1,A3),A1).*phi1z + dot(cross(A2,A3),A1).*phi2d1z + dot(cross(A2d1,A3),A1).*phi2z;
    
    %SIGMA_22
	B2(2,1:3:end) = A3d2(1).*dR1B(2,:) + dot(cross(A1,A3),A2).*phi1d2x + dot(cross(A1d2,A3),A2).*phi1x + dot(cross(A2d2,A3),A2).*phi2x;
	B2(2,2:3:end) = A3d2(2).*dR1B(2,:) + dot(cross(A1,A3),A2).*phi1d2y + dot(cross(A1d2,A3),A2).*phi1y + dot(cross(A2d2,A3),A2).*phi2y; 
	B2(2,3:3:end) = A3d2(3).*dR1B(2,:) + dot(cross(A1,A3),A2).*phi1d2z + dot(cross(A1d2,A3),A2).*phi1z + dot(cross(A2d2,A3),A2).*phi2z;

    %SIGMA_12
    B2(3,1:3:end) = A3d2(1).*dR1B(1,:) + dot(cross(A1,A3),A2).*phi1d1x + dot(cross(A1d1,A3),A2).*phi1x + dot(cross(A2d1,A3),A2).*phi2x +...
                                    A3d1(1).*dR1B(2,:) + dot(cross(A1d2,A3),A1).*phi1x + dot(cross(A2,A3),A1).*phi2d2x + dot(cross(A2d2,A3),A1).*phi2x;
    B2(3,2:3:end) = A3d2(2).*dR1B(1,:) + dot(cross(A1,A3),A2).*phi1d1y + dot(cross(A1d1,A3),A2).*phi1y + dot(cross(A2d1,A3),A2).*phi2y +...
                                    A3d1(2).*dR1B(2,:) + dot(cross(A1d2,A3),A1).*phi1y + dot(cross(A2,A3),A1).*phi2d2y + dot(cross(A2d2,A3),A1).*phi2y; 
    B2(3,3:3:end) = A3d2(3).*dR1B(1,:) + dot(cross(A1,A3),A2).*phi1d1z + dot(cross(A1d1,A3),A2).*phi1z + dot(cross(A2d1,A3),A2).*phi2z +...
                                    A3d1(3).*dR1B(2,:) + dot(cross(A1d2,A3),A1).*phi1z + dot(cross(A2,A3),A1).*phi2d2z + dot(cross(A2d2,A3),A1).*phi2z; 
       
	K_loc = K_loc + (Nurbs2D.wgp{kk}(i)*J) .* ((transpose(B1)*(((E*t)/(1-nu^2)).*C)*B1) + (transpose(B2)*(((E*t^3)/(12*(1-nu^2))).*C)*B2));
        
end