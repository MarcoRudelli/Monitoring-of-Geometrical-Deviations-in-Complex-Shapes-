function [A, B, indList] = LinearFEM(mesh, boundary)

%%
% Compute linear FEM for triangulations based on 
%   "Laplace Spectra for Shape Recognition" by Martin Reuter
%
% Inputs:
% mesh: mesh structure as defined in the "matlabmesh" package. 
% boundary: boundary condition to apply. 'D' for Dirichlet, 'N' or default
% for Neumann.
% 
% Outputs:
% A, B: generalized eigenvalue problem is AU=lambda*BU. 
% indList: indices of the original vertices corresponding to the leading rows and
% columns in A and B. For Neumann condition, indList=1:npoints.
%

%% Step 0: Read from mesh
v = mesh.v;
f = mesh.f;
nfaces = size(f,1);

%% Step 1: Local A and B entries on the unity triangle 
localAuu = [1/2 -1/2 0; -1/2 1/2 0; 0 0 0];
localAvv = [1/2 0 -1/2; 0 0 0; -1/2 0 1/2];
localAuv = [1/2 0 -1/2; -1/2 0 1/2; 0 0 0];
localAvu = [1/2 -1/2 0; 0 0 0; -1/2 1/2 0];
%localAuvvu = localAuv+localAvu;
localB = (ones(3)+eye(3))/24;

%% Step 2: Compute tau1=P2-P1 and tau2=P3-P1 for each triangle
P1 = v(f(:,1),:); 
P2 = v(f(:,2),:); 
P3 = v(f(:,3),:);

tau1 = P2-P1;
tau2 = P3-P1;

tau11 = sum(tau1.^2,2); % dot product of tau1 and tau1
tau22 = sum(tau2.^2,2); % dot product of tau2 and tau2
tau12 = sum(tau1.*tau2,2); % dot product of tau1 and tau2
W = sqrt(tau11.*tau22-tau12.^2); % sqrt of det(g), g is the metric tensor
% sum(W==0) % W=0 means tau1 and tau2 parallel. Shouldn't happen in triangles

%% Step 3: A and B entries in all local triangles
% each triangle-wise value is repeated three times because local A and B are 3*3
t11 = repmat(tau11,1,3);
t11 = reshape(t11',3*nfaces,1);

t22 = repmat(tau22,1,3);
t22 = reshape(t22',3*nfaces,1);

t12 = repmat(tau12,1,3);
t12 = reshape(t12',3*nfaces,1);

Wrep = repmat(W,1,3);
Wrep = reshape(Wrep',3*nfaces,1);

% each local matrix is repeated nfaces times
Auu = repmat(localAuu,nfaces,1).*t22;
Avv = repmat(localAvv,nfaces,1).*t11;
Auv = repmat(localAuv,nfaces,1).*t12;
Avu = repmat(localAvu,nfaces,1).*t12;
Aentries = (Auu+Avv-Auv-Avu)./Wrep;
Bentries = repmat(localB,nfaces,1).*Wrep;

%% Step 4: Find global indices in Aentries and Bentries
I_ind = repmat(reshape(f',3*nfaces,1),3,1);
temp = repmat(f,1,1,3);
temp = permute(temp,[3,1,2]);
J_ind = reshape(temp,9*nfaces,1);

%% Step 5: Construct sparse matrices
A = sparse(I_ind, J_ind, reshape(Aentries,9*nfaces,1));
B = sparse(I_ind, J_ind, reshape(Bentries,9*nfaces,1));
indList = 1:size(v,1);

% apply the Dirichlet boundary condition if needed
if exist('boundary','var') && boundary=='D'
    [row, col] = find(mesh.e==1);  
    BoundaryPts = unique([row;col]); % all points on boundary edges
    A(BoundaryPts,:) = [];
    A(:,BoundaryPts) = [];
    B(BoundaryPts,:) = [];
    B(:,BoundaryPts) = [];
    indList(BoundaryPts) = [];
end

