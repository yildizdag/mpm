function M = globalMass2D_cylindricalShell(Nurbs2D)
local_dof = Nurbs2D.local_dof_dry;
storeSparse = ((Nurbs2D.nen{1})^2)*Nurbs2D.nel{1};
% STORING:
I = zeros(storeSparse,1);
J = zeros(storeSparse,1);
M = zeros(storeSparse,1);
ntriplets = 0;
for kk = 1:1
    for el = 1:Nurbs2D.nel{kk}
        iu = Nurbs2D.INC{kk}(Nurbs2D.IEN{kk}(1,el),1);   
        iv = Nurbs2D.INC{kk}(Nurbs2D.IEN{kk}(1,el),2);
        lm_loc  = reshape(fliplr(Nurbs2D.LM{kk}(:,:,el)),local_dof*Nurbs2D.nen{kk},1);
        [i,~,s] = find(lm_loc);      
        M_loc = localMass2D_cylindricalShell(iu,iv,Nurbs2D,kk);
        for k = 1:numel(i)
            for l = 1:numel(i)
                ntriplets = ntriplets + 1;
                I(ntriplets) = s(k);
                J(ntriplets) = s(l);
                M(ntriplets) = M_loc(i(k),i(l));
            end
        end
    end
end
M = sparse(I(1:ntriplets),J(1:ntriplets),M(1:ntriplets),Nurbs2D.dryDofs,Nurbs2D.dryDofs);
M = (M+M')./2;