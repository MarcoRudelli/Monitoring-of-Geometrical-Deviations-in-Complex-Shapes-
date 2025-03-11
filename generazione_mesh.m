% This code contains examples to apply the modified DFEWMA chart to the LB
% spectra

%%{
%% Case I: simulated parts in folder "Sequence of parts for SPC"
% parameters

numEigs = 16; % # of eigenvalues
numParts = 150;
Data = zeros(numParts*100,numEigs);

noise = 0.04;
Data2 = zeros(numParts*100,1);

%sigma12 = 0.02^2; sigma22 = noise^2-sigma12;
%rx = 5; ry = 5; rz = 5;

mesh1 = readMesh('CAD-model-Cut.obj');
model = mesh1.v;

%model1 = mesh1.v;
mesh2 = readMesh('CAD-model-Cut_missing_strut.obj');
model2 = mesh2.v;
%[row2, col2] = find(model ~= model2);
%diff2 = model(row2,:) - model2(row2,:);
%model2(row2,:) = model2(row2,:)+0.95*diff2;
%mesh2.v = model2;

for j=1:100
    % calculate LB spectra
    for i=1:numParts
        if(i <= 100)
            mesh = mesh1; %readMesh(['PartNo', num2str(i),'.obj']); % matlabmesh package needed
        else
            mesh = mesh2;
        end
        P = mesh.v; P = P+noise*randn(size(P)); mesh.v = P; % add noise
    
        [A1, B1] = LinearFEM(mesh); % B is N*1 here
        % L_tilde=B^(-1)L shares the same eigenvalue with B^(-1/2)LB^(-1/2)
        % L_symm = sqrt(B).*L_tilde./sqrt(B'); 
        Data(i+(j-1)*150,:) = real(eigs(A1, B1, numEigs, -1e-8, 'IsSymmetricDefinite', false));

        mesh_b = mesh;
        part = mesh_b.v;
        
        % apply isometric transformations and permutation of point indices to
        % simulate real scanned objects; this step is not necessary to the LB
        % spectrum method as it is isometric invariant:
        per = randperm(size(part,1));
        per_inv = NaN(1, size(part,1));      
        per_inv(per) = 1:size(part,1); 
        part = part(per,:); % random permutation of point indices
        v1 = 0.01*(2*rand-1); v2 = 0.01*(2*rand-1); v3 = 0.01*(2*rand-1); % rotation angles
        R1 = [1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
        R2 = [cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
        R3 = [cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];
        part = part*R3'*R2'*R1'; % apply rotation transformation
        part = part + randn(1,3); % apply translation transformation
        
        % apply the ICP algorithm (ICP package needed)
        [~,~,icp_part,res,~]=icp(model,part); % res is the "smallest" deviation from target
        % res needs to be manually added as the second last output of function "icp"
    
        V_mat = vertcat([size(icp_part(:,per_inv), 2), 0, 0], transpose(icp_part(:,per_inv)));
        T_mat = vertcat([size(mesh_b.f, 1), 0, 0], mesh_b.f - 1);
        writematrix(V_mat, ['/Volumes/EXTERNAL_USB/parts_simulated_SIM2.01/SCAN_verts_', num2str(i+(j-1)*150),'.csv']);
        writematrix(T_mat, ['/Volumes/EXTERNAL_USB/parts_simulated_SIM2.01/SCAN_triangles_', num2str(i+(j-1)*150),'.csv']);
        Data2(i+(j-1)*150) = res;
    
        i+(j-1)*150
    end

    writematrix(Data, ['/Volumes/EXTERNAL_USB/parts_simulated_SIM2.01/parts_simulated_res_LBS_SIM2.01_', num2str(i),'.csv']);
    writematrix(Data2, ['/Volumes/EXTERNAL_USB/parts_simulated_SIM2.01/parts_simulated_res_ICP_SIM2.01_', num2str(i),'.csv']);
end


%{
% parameters for the DFEWMA chart
minwin = 1; maxwin = 10; dim = numEigs; m0 = 100; kmax = numParts-m0;
tau = 100; delta = 0; alpha = 0.005; nbound = 20000; lambda = 0.01;
[runlength, climit, Tval] = RobustDFEWMA_sd(minwin, maxwin, dim, kmax, ...
    m0, tau, delta, alpha, nbound, lambda, Data);
runlength

% plot the chart
figure
plot(Tval, 'k*-'); hold on
plot(climit, 'r:','LineWidth',2); hold on
line([runlength runlength],[0 Tval(runlength)],'Color','blue',...
    'LineWidth',2);
text(runlength, -2, num2str(runlength),'Color','blue');
legend({'Observed test statistic','Control limit'})
%}
%%}

%{
%% Case II: cylindrical parts with uncorrelated isotropic noise
% parameters
num_sim = 100;
r = 10; h = 50; noise = 0.05;
numEigs = 16; % # of eigenvalues
numParts = 150*num_sim;
Data = zeros(numParts,numEigs);
eigsopts = struct('disp',0);

sigma12 = 0.02^2; sigma22 = noise^2-sigma12;
rx = 2.6; ry = 2.6; rz = 16.7;

Data2 = zeros(numParts,1);
% a noise-free big mesh as the model to register the part
CyMesh1 = CylinderSimulation(r, h, 0, noise, 20000); 
model = CyMesh1.v; 
% writematrix(CyMesh1.v, 'cylinder_cad_verts.csv');
% writematrix(CyMesh1.f, 'cylinder_cad_triangles.csv');
% writeMesh(CyMesh1, 'cylinder_cad.obj')

% simulate cylindrical parts and calculate their LB spectra
for i=1:numParts
    % delta = chi2rnd(3) / 10000;
    if i<=100*num_sim % simulate in-control parts
        delta = 0;
        r = 10;
    else
        delta = 0; % + 0.0005;
        r = 10.01;
    end
    CyMesh = CylinderSimulation(r, h, delta, noise, 1994+randi(11));
    
    %%% % add noise
    %%% P = CyMesh.v; P = P+noise*randn(size(P)); CyMesh.v = P; 

     % add spatially correlated, non-isotropic noise
    P = CyMesh.v; numPoints = size(P,1);
    Sigma_x = sigma12*exp(-abs(P(:,1)-P(:,1)')/rx)+sigma22*eye(numPoints);
    Sigma_y = sigma12*exp(-abs(P(:,2)-P(:,2)')/ry)+sigma22*eye(numPoints);
    Sigma_z = sigma12*exp(-abs(P(:,3)-P(:,3)')/rz)+sigma22*eye(numPoints);
    P(:,1) = P(:,1)+mvnrnd(zeros(1,numPoints),Sigma_x,1)';
    P(:,2) = P(:,2)+mvnrnd(zeros(1,numPoints),Sigma_y,1)';
    P(:,3) = P(:,3)+mvnrnd(zeros(1,numPoints),Sigma_z,1)';
    CyMesh.v=P;


    [A1, B1] = LinearFEM(CyMesh); % B is N*1 here
    Data(i,:) = eigs(A1, B1, numEigs, -1e-8, 'IsSymmetricDefinite', 1);
    
   


    CyMesh2 = CyMesh;
    part = CyMesh2.v;
    
    % apply isometric transformations and permutation of point indices to
    % simulate real scanned objects; this step is not necessary to the LB
    % spectrum method as it is isometric invariant:
    per = randperm(size(part,1));
    per_inv = NaN(1, size(part,1));      
    per_inv(per) = 1:size(part,1); 
    part = part(per,:); % random permutation of point indices
    v1 = 0.6*(2*rand-1); v2 = 0.6*(2*rand-1); v3 = 0.6*(2*rand-1); % rotation angles
    R1 = [1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
    R2 = [cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
    R3 = [cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];
    part = part*R3'*R2'*R1'; % apply rotation transformation
    part = part + randn(1,3); % apply translation transformation
    
    % apply the ICP algorithm (ICP package needed)
    [~,~,icp_part,res,~]=icp(model,part); % res is the "smallest" deviation from target
    % res needs to be manually added as the second last output of function "icp"

    V_mat = vertcat([size(icp_part(:,per_inv), 2), 0, 0], transpose(icp_part(:,per_inv)));
    T_mat = vertcat([size(CyMesh2.f, 1), 0, 0], CyMesh2.f - 1);
    writematrix(V_mat, ['/Volumes/EXTERNAL_USB/sim_2.11/cylinder_simulated_SIM2.11/SCAN_verts_', num2str(i),'.csv']);
    writematrix(T_mat, ['/Volumes/EXTERNAL_USB/sim_2.11/cylinder_simulated_SIM2.11/SCAN_triangles_', num2str(i),'.csv']);
    Data2(i) = res;

    i
end
writematrix(Data, ['/Volumes/EXTERNAL_USB/sim_2.11/cylinder_simulated_SIM2.11/cylinder_simulated_res_LBS_SIM2.11_', num2str(i),'.csv']);
writematrix(Data2, ['/Volumes/EXTERNAL_USB/sim_2.11/cylinder_simulated_SIM2.11/cylinder_simulated_res_ICP_SIM2.11_', num2str(i),'.csv']);

%{
% parameters for the DFEWMA chart
minwin = 1; maxwin = 10; dim = numEigs; m0 = 100; kmax = numParts-m0;
tau = 100; delta = 0; alpha = 0.005; nbound = 20000; lambda = 0.01;
[runlength, climit, Tval] = RobustDFEWMA_sd(minwin, maxwin, dim, kmax, ...
    m0, tau, delta, alpha, nbound, lambda, Data);
runlength

% plot the chart
figure
plot(Tval, 'k*-'); hold on
plot(climit, 'r:','LineWidth',2); hold on
line([runlength runlength],[0 Tval(runlength)],'Color','blue',...
    'LineWidth',2);
text(runlength, -2, num2str(runlength),'Color','blue');
legend({'Observed test statistic','Control limit'})


%ICP
dim = 1;
[runlength, climit, Tval] = RobustDFEWMA_sd(minwin, maxwin, dim, kmax, ...
    m0, tau, delta, alpha, nbound, lambda, Data2);
runlength

% plot the chart
figure
plot(Tval, 'k*-'); hold on
plot(climit, 'r:','LineWidth',2); hold on
line([runlength runlength],[0 Tval(runlength)],'Color','blue',...
    'LineWidth',2);
text(runlength, -2, num2str(runlength),'Color','blue');
legend({'Observed test statistic','Control limit'})
%}

%}

%{
%% Case III: cylindrical parts with spatially correlated, non-isotropic noise
% parameters
r = 10; h = 50; noise = 0.05; 
sigma12 = 0.02^2; sigma22 = noise^2-sigma12;
rx = 2.6; ry = 2.6; rz = 16.7; 
numEigs = 16; % # of eigenvalues
numParts = 150*100;
Data = zeros(numParts,numEigs);
eigsopts = struct('disp',0);


Data2 = zeros(numParts,1);
% a noise-free big mesh as the model to register the part
CyMesh1 = CylinderSimulation(r, h, 0, noise, 20000); 
model = CyMesh1.v; 

% simulate cylindrical parts and calculate their LB spectra
for i=1:numParts
    if i<=100*100 % simulate in-control parts
        delta = 0;
    else
        delta = 0.0005;
    end
    CyMesh = CylinderSimulation(r, h, delta, noise, 1994+randi(11));
    % add spatially correlated, non-isotropic noise
    P = CyMesh.v; numPoints = size(P,1);
    Sigma_x = sigma12*exp(-abs(P(:,1)-P(:,1)')/rx)+sigma22*eye(numPoints);
    Sigma_y = sigma12*exp(-abs(P(:,2)-P(:,2)')/ry)+sigma22*eye(numPoints);
    Sigma_z = sigma12*exp(-abs(P(:,3)-P(:,3)')/rz)+sigma22*eye(numPoints);
    P(:,1) = P(:,1)+mvnrnd(zeros(1,numPoints),Sigma_x,1)';
    P(:,2) = P(:,2)+mvnrnd(zeros(1,numPoints),Sigma_y,1)';
    P(:,3) = P(:,3)+mvnrnd(zeros(1,numPoints),Sigma_z,1)';
    CyMesh.v=P;
    [A1, B1] = LinearFEM(CyMesh); % B is N*1 here
    % L_tilde=B^(-1)L shares the same eigenvalue with B^(-1/2)LB^(-1/2)
    % L_symm = sqrt(B).*L_tilde./sqrt(B'); 
    Data(i,:) = eigs(A1, B1, numEigs, -1e-8, 'IsSymmetricDefinite', 1);


    CyMesh2 = CyMesh;
    part = CyMesh2.v;

    % apply isometric transformations and permutation of point indices to
    % simulate real scanned objects; this step is not necessary to the LB
    % spectrum method as it is isometric invariant:
    per = randperm(size(part,1));
    per_inv = NaN(1, size(part,1));      
    per_inv(per) = 1:size(part,1); 
    part = part(per,:); % random permutation of point indices
    v1 = 0.6*(2*rand-1); v2 = 0.6*(2*rand-1); v3 = 0.6*(2*rand-1); % rotation angles
    R1 = [1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
    R2 = [cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
    R3 = [cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];
    part = part*R3'*R2'*R1'; % apply rotation transformation
    part = part + randn(1,3); % apply translation transformation
    
    % apply the ICP algorithm (ICP package needed)
    [~,~,icp_part,res,~]=icp(model,part); % res is the "smallest" deviation from target
    % res needs to be manually added as the second last output of function "icp"
    
    V_mat = vertcat([size(icp_part(:,per_inv), 2), 0, 0], transpose(icp_part(:,per_inv)));
    T_mat = vertcat([size(CyMesh2.f, 1), 0, 0], CyMesh2.f - 1);
    writematrix(V_mat, ['cylinder_simulated_SIM2.2/SCAN_verts_', num2str(i),'.csv']);
    writematrix(T_mat, ['cylinder_simulated_SIM2.2/SCAN_triangles_', num2str(i),'.csv']);
    Data2(i) = res;

    i
end
writematrix(Data, ['cylinder_simulated_SIM2.2/cylinder2_simulated_res_LBS_SIM2.2_', num2str(i),'.csv']);
writematrix(Data2, ['cylinder_simulated_SIM2.2/cylinder2_simulated_res_ICP_SIM2.2_', num2str(i),'.csv']);

%{
% parameters for the DFEWMA chart
minwin = 1; maxwin = 10; dim = numEigs; m0 = 100; kmax = numParts-m0;
tau = 100; delta = 0; alpha = 0.005; nbound = 20000; lambda = 0.01;
[runlength, climit, Tval] = RobustDFEWMA_sd(minwin, maxwin, dim, kmax, ...
    m0, tau, delta, alpha, nbound, lambda, Data);
runlength

% plot the chart
figure
plot(Tval, 'k*-'); hold on
plot(climit, 'r:','LineWidth',2); hold on
line([runlength runlength],[0 Tval(runlength)],'Color','blue',...
    'LineWidth',2);
text(runlength, -2, num2str(runlength),'Color','blue');
legend({'Observed test statistic','Control limit'})


%ICP
dim = 1;
[runlength, climit, Tval] = RobustDFEWMA_sd(minwin, maxwin, dim, kmax, ...
    m0, tau, delta, alpha, nbound, lambda, Data);
runlength

% plot the chart
figure
plot(Tval, 'k*-'); hold on
plot(climit, 'r:','LineWidth',2); hold on
line([runlength runlength],[0 Tval(runlength)],'Color','blue',...
    'LineWidth',2);
text(runlength, -2, num2str(runlength),'Color','blue');
legend({'Observed test statistic','Control limit'})
%}

%}
