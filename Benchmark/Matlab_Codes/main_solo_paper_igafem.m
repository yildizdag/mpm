%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN CODE for SHELL FREE VIBRATION
% KIRCHHOFF-LOVE MODEL
% HALF CYLINDER w/ ANTISYMMETRIC BC
%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the Geometry imported from Rhino:
FileName = 'lindholm_5000_';
numDryPatch = 1;
%Degrees of Freedom per each node:
local_dof_dry = 3;
%Number of Modes to extract:
modeNum = 24;

%Material Properties:
E = 2.1E11;
nu = 0.3;
rho = 7800;
t = 0.00484;

%CREATE 2D IGA MESH (reads FileName):
Nurbs2D = iga2DmeshDry(FileName,numDryPatch,local_dof_dry);
%%%%%%%%%%%%%%%%%%%%%
%iga2DmeshPlotNURBS(Nurbs2D);

%Save Inputs:
Nurbs2D.E = E;
Nurbs2D.nu = nu;
Nurbs2D.rho = rho;
Nurbs2D.t = t;
Nurbs2D.modeNum = modeNum;
%%%%%%%%%%%%%%%%%%%%%

%Save and Load .mat File:
save('Nurbs2D.mat','Nurbs2D');
clear;clc;
load('Nurbs2D.mat');
%%%%%%%%%%%%%%%%%%%%%
%VACUUM ANALYSIS:
K = globalStiffness2D_cylindricalShell(Nurbs2D);
M = globalMass2D_cylindricalShell(Nurbs2D);

%APPLY BOUNDARY CONDITIONS
BCPoints = Nurbs2D.nnp{1}-2*Nurbs2D.number{1}(1)+1:Nurbs2D.nnp{1};
BounNodes = unique([3.*BCPoints-2, 3.*BCPoints-1, 3.*BCPoints]);

%Solution
modeNum = Nurbs2D.modeNum;
%Apply BCs:
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
[V,freq] = eigs(K,M,modeNum,'sm');
freq2 = (sqrt(freq))/(2*pi);
freq2 = (diag(freq2));