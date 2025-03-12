% Generazione di 150*100 mesh 

numEigs = 15; % numero autovalori oltre al primo (il primo è sempre nullo)
numSim = 10; % numero di replicazioni simulate
numParts = 150; % numero totale di unità per ogni simulazione
m0 = 100; % numero di unità IC

Data = zeros(numParts, numEigs+1);

noise = 0.1; % sd errore casuale spazialmente incorrelato

% parametri errore spazialmente correlato
% sigma12 = 0.02^2; sigma22 = noise^2-sigma12;
% rx = 5; ry = 5; rz = 5;

q = 0.95; % parametro per regolare l'entità del difetto OC

mesh1 = readMesh('acceptable.obj');
model = mesh1.v;

mesh2 = readMesh('chipped1.obj');
model2 = mesh2.v;
[row2, col2] = find(model ~= model2);
diff2 = model(row2,:) - model2(row2,:);
model2(row2,:) = model2(row2,:) + q * diff2;
mesh2.v = model2;

for j=1:numSim
    if j==1
        mkdir parts_simulated
        mkdir results_LBS
    end

    for i=1:numParts
        if(i <= m0)
            mesh = mesh1; % IC
        else
            mesh = mesh2; % OC
        end

        P = mesh.v; P = P+noise*randn(size(P)); mesh.v = P; % aggiunta rumore
    
        % calcolo autovalori LBS
        [A1, B1] = LinearFEM(mesh);
        Data(i,:) = real(eigs(A1, B1, numEigs+1, -1e-8, 'IsSymmetricDefinite', false));

        mesh_b = mesh;
        part = mesh_b.v;
        
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
        % https://www.mathworks.com/matlabcentral/fileexchange/12627-iterative-closest-point-method?s_tid=prof_contriblnk
        % Need to add “res” and “vi” as the two additional outputs from function “icp”
        [~,~,icp_part,res,~]=icp(model,part); % res is the "smallest" deviation from target
        % res needs to be manually added as the second last output of function "icp"
    
        % creo matrici vertici e triangoli (la prima riga contiene il numero dei vertici/numero dei triangoli)
        % ripristinando la permutazione iniziale dei vertici
        V_mat = vertcat([size(icp_part(:,per_inv), 2), 0, 0], transpose(icp_part(:,per_inv)));
        T_mat = vertcat([size(mesh_b.f, 1), 0, 0], mesh_b.f - 1);

        writematrix(V_mat, ['parts_simulated/SCAN_verts_', num2str(i+(j-1)*150),'.csv']);
        writematrix(T_mat, ['parts_simulated/SCAN_triangles_', num2str(i+(j-1)*150),'.csv']);
    
        i+(j-1)*150
    end

    writematrix(Data(:,2:(numEigs+1)), ['results_LBS/parts_simulated_res_LBS_', num2str(j),'.csv']);
end



