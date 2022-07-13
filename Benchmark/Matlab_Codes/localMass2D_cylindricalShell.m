function M_loc = localMass2D_cylindricalShell(iu,iv,Nurbs2D,kk)
local_dof = Nurbs2D.local_dof_dry;
M_loc = zeros(local_dof*Nurbs2D.nen{kk},local_dof*Nurbs2D.nen{kk});
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
    %Covariant Basis Vectors:
    A1 = dS(:,2,1);
    A2 = dS(:,1,2);
    J1 = norm(cross(A1,A2));
    A3 = cross(A1,A2)/J1;
    % From Parametric to Parent:
    J2 = diag([Nurbs2D.knots.U{kk}(iu+1) - Nurbs2D.knots.U{kk}(iu),...
                            Nurbs2D.knots.V{kk}(iv+1) - Nurbs2D.knots.V{kk}(iv)])/2;
    % Jacobian:
    J = J1*det(J2);

    R1 = reshape(dR(:,:,1,1),1,Nurbs2D.nen{kk});
    N = zeros(3,local_dof*Nurbs2D.nen{kk});
    N(1,1:3:end) = R1; N(2,2:3:end) = R1; N(3,3:3:end) = R1;
    
    % First Derivatives:
%     dR1B = zeros(2,Nurbs2D.nen{kk});
%     dR1B(1,:) = reshape(dR(:,:,2,1),1,Nurbs2D.nen{kk});
%     dR1B(2,:) = reshape(dR(:,:,1,2),1,Nurbs2D.nen{kk});
    
%     B1 = zeros(1,local_dof*Nurbs2D.nen{kk});
%     B2 = zeros(1,local_dof*Nurbs2D.nen{kk});  
%     B1(1,1:3:end) = dR1B(1,:).*A3(1);
%     B1(1,2:3:end) = dR1B(1,:).*A3(2);
%     B1(1,3:3:end) =  dR1B(1,:).*A3(3);
%     B2(1,1:3:end) = dR1B(2,:).*A3(1);
%     B2(1,2:3:end) = dR1B(2,:).*A3(2);
%     B2(2,3:3:end) =  dR1B(2,:).*A3(3);

    C = [Nurbs2D.rho, 0, 0; 0, Nurbs2D.rho, 0; 0, 0, Nurbs2D.rho];
%     C1 = dot(cross(A2,A3),cross(A2,A3))*Nurbs2D.rho/J1;
%     C2 = dot(cross(A1,A3),cross(A1,A3))*Nurbs2D.rho/J1;

    %M_loc = M_loc + (Nurbs2D.wgp{kk}(i)*J) .* ( (transpose(N)*(t .* C)*N) + (transpose(N1)*((t^3/12) .* Nurbs2D.rho)*N1) + (transpose(N2)*((t^3/12) .* Nurbs2D.rho)*N2) );
    M_loc = M_loc + (Nurbs2D.wgp{kk}(i)*J) .* ( (transpose(N)*(t .* C)*N) );
    %M_loc = M_loc + (Nurbs2D.wgp{kk}(i)*J) .* ( transpose(N)*((t.*C)*N) + transpose(B1)*(((t^3/12)*C1)*B1) + transpose(B2)*(((t^3/12)*C2)*B2));

end