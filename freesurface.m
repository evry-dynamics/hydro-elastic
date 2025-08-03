clear
clc

WaveType = 'symmetric';
ScaleFacor = 0.01;

NodeFile = 'NLIST_10_10_0.36.lis';
ElementFile = 'ELIST_10_10_0.36.lis';
[node_s,element_s] = ReadFile_s(NodeFile,ElementFile);
node_s(:,2:4) = node_s(:,2:4)*ScaleFacor;

g = 9.80665; % gravitational acceleration (m/s^2)
hh = 0.2; % fluid depth(m)


NodeFile = 'NLIST_f_10_10_5.lis';
ElementFile = 'ELIST_f_10_10_5.lis';
% NodeFile = 'NLIST_f_20_20_0.1m.lis';
% ElementFile = 'ELIST_f_20_20_0.1m.lis';
[node_f,element_f] = ReadFile(NodeFile,ElementFile);
% 通过z方向的平移，将流体从固体下方移到固体上方
% node_f(:,4) = node_f(:,4) + 2; 
node_f(:,2:4) = node_f(:,2:4)*ScaleFacor;
% 调整流体深度
node_f(:,4) = node_f(:,4)*2;

%% Material
% plate (PVC-聚氯乙烯Polyvinyl chloride)
Es1 = 6e6; % young's modulus
vs1 = 0.47; % poisson's ratio
rho_s1 = 1.29*10^(3); % density

% % stub (轧制纯铜)
% Es2 = 108e9; % young's modulus
% vs2 = 0.33; % poisson's ratio
% rho_s2 = 8.9*10^(3); % density

% % plate (epoxy)
% Es1 = 3.3e6; % young's modulus
% vs1 = 0.33; % poisson's ratio
% rho_s1 = 1*10^(3); % density
% 
% 
% stub (wu)
Es2 = 411e9; % young's modulus
vs2 = 0.28; % poisson's ratio
rho_s2 = 19.3*10^(3); % density

% lame constant for two material
lambda_s1 = Es1*vs1/((1+vs1)*(1-2*vs1)); 
nu_s1 = Es1/(2*(1+vs1)); 
lambda_s2 = Es2*vs2/((1+vs2)*(1-2*vs2)); 
nu_s2 = Es2/(2*(1+vs2)); 
% % the material of solid plate (Aluminum)
% Es1 = 70*10^9; % young's modulus
% vs1 = 0.33; % poisson's ratio
% rho_s1 = 2.7*10^(3); % density

% the material of fluid domain (water)
rho_f = 1000; % Fluid density
% rho_f = 2280; % Fluid density
c_f   = 1460; % Fluid acoustic velocity
% g = 9.80665; % gravitational acceleration (m/s^2)
%% Setting
% mode order
n = 8;
% plate thickness and stub height
% t = 1.5*ScaleFacor;
t = 2*ScaleFacor;
h = 0;
% h = 10*ScaleFacor;
% imaginary unit
I = sqrt(-1);
% -------------------------------------------------------------------------
%                  Node classification and rearrangement
% -------------------------------------------------------------------------
% Solid domain
[nodeS_FSI,DOF_INTER,DOF_EDGE_SORT,DOF_INDENPENT,DOF_A_S,...
    edge_dof_num,index_L_row,index_B_row,index_BL_row,index_R_row,index_T_row,index_INTER_row,...
    index_BR_row,index_TL_row,index_TR_row,index_L_col,index_B_col,index_BL_col,index_INTER_col]...
    = SortNode_S(node_s);
% Fluid domain
% 找出自由液面的节点，并且按照逆时针进行排列，组成自由液面4边形单元
[nodeF_FSI,node_Freesurface,element_Freesurface,node_FG,DOF_INTER_F,DOF_EDGE_SORT_F,DOF_INDENPENT_F,DOF_A_F,...
    edge_dof_num_F,index_L_row_F,index_B_row_F,index_BL_row_F,index_R_row_F,index_T_row_F,...
    index_BR_row_F,index_TL_row_F,index_TR_row_F,index_INTER_row_F,...
    index_L_col_F,index_B_col_F,index_BL_col_F,index_INTER_col_F] ...
    = SortNode_F(node_f);
% Free surface
[DOF_INTER_G,DOF_EDGE_SORT_G,DOF_INDENPENT_G,DOF_A_G,...
    edge_dof_num_G,index_L_row_G,index_B_row_G,index_BL_row_G,index_R_row_G,index_T_row_G,index_INTER_row_G,...
    index_BR_row_G,index_TL_row_G,index_TR_row_G,index_L_col_G,index_B_col_G,index_BL_col_G,index_INTER_col_G]...
    = SortNode_G(node_Freesurface);
% FSI Matrix
R = FSIMatrix(node_s,nodeS_FSI,node_f,nodeF_FSI);

%% Degree of freedom
% Solid 
dof_s = 5*length(node_s);
dof_f = length(node_f);
dof_g = length(node_Freesurface);

%% Mesh
% PlotPhononicsquare(node_s,element_s)
% PlotPhononicsquare_G(node_Freesurface,element_Freesurface)
% % PlotPhononicTriangular(node,element)
% % ansys mesh: mat 2-blue, 1-yellow.
% % comsol mesh: plate-1-b, stub-2-y



%% Calculate fluid stiffness and mass matrix
[Kf,Mf] = FEMSolver_Fluid(node_f,element_f,c_f);

%% Calculate free surface stiffness and mass matrix
[Kg,Mg] = FEMSolver_Freesurface(dof_f,dof_g,node_FG,node_Freesurface,element_Freesurface,rho_f,g);

%% Calculate solid stiffness and mass matrix
Ks = sparse(dof_s,dof_s); % inital global stiffness matrix
Ms = sparse(dof_s,dof_s); % inital global mass matrix
for i = 1:size(element_s,1)
    x1 = node_s(element_s(i,2),2);
    x2 = node_s(element_s(i,3),2);
    x3 = node_s(element_s(i,4),2);
    x4 = node_s(element_s(i,5),2);

    y1 = node_s(element_s(i,2),3);
    y2 = node_s(element_s(i,3),3);
    y3 = node_s(element_s(i,4),3);
    y4 = node_s(element_s(i,5),3);
    % material identifier
%     mattype = element(i,5);
    thick = element_s(i,6);
    % degree of freedom of element vector(index of assemble global matrix)
    edof = [element_s(i,2)*5-4,element_s(i,2)*5-3,element_s(i,2)*5-2,element_s(i,2)*5-1,element_s(i,2)*5,...
    element_s(i,3)*5-4,element_s(i,3)*5-3,element_s(i,3)*5-2,element_s(i,3)*5-1,element_s(i,3)*5,...
    element_s(i,4)*5-4,element_s(i,4)*5-3,element_s(i,4)*5-2,element_s(i,4)*5-1,element_s(i,4)*5,...
    element_s(i,5)*5-4,element_s(i,5)*5-3,element_s(i,5)*5-2,element_s(i,5)*5-1,element_s(i,5)*5];


    switch thick
        case 1 % stub
            nu = nu_s2;
            rho = rho_s2;
            v = vs2;
            E = Es2;
            Db = GeneralizedElasticMatrix(E,v,nu,t,h,'stub','Bend');
            Ds = GeneralizedElasticMatrix(E,v,nu,t,h,'stub','Shear'); 
            kb = ElementStiffnessMatrix(x1,x2,x3,x4,y1,y2,y3,y4,Db,'Bend');
            ks = ElementStiffnessMatrix(x1,x2,x3,x4,y1,y2,y3,y4,Ds,'Shear');
            k = kb + ks;
            m = ElementMassMatrix(x1,x2,x3,x4,y1,y2,y3,y4,rho,t,h,'stub');

        case 2 % plate
            nu = nu_s1;
            rho = rho_s1;
            v = vs1;
            E = Es1;
            Db = GeneralizedElasticMatrix(E,v,nu,t,h,'plate','Bend');
            Ds = GeneralizedElasticMatrix(E,v,nu,t,h,'plate','Shear'); 
            kb = ElementStiffnessMatrix(x1,x2,x3,x4,y1,y2,y3,y4,Db,'Bend');
            ks = ElementStiffnessMatrix(x1,x2,x3,x4,y1,y2,y3,y4,Ds,'Shear');
            k = kb + ks;
            m = ElementMassMatrix(x1,x2,x3,x4,y1,y2,y3,y4,rho,t,h,'plate');
            
    end
    Ks(edof,edof) = Ks(edof,edof) + k;
    Ms(edof,edof) = Ms(edof,edof) + m;
end

%% Calculate Bloch Vector
x_max = max(node_s(:,2)); 
y_max = max(node_s(:,3));
% kx_vector = [0,pi/(x_max*40):pi/(x_max*40):pi/(x_max)];  % bloch vector in x direction
kx_vector = [0:0.001:0.1];
ky_vector = [(zeros(1,length(kx_vector)))];  % bloch vector in y direction

% [kx_vector,ky_vector] = CalculateBlochVector(node_s,WaveType);
% kx_vector = kx_vector(2:2:60); 
% ky_vector = ky_vector(2:2:60);
%% 施加周期性边界
dof_s = edge_dof_num;
dof_g = edge_dof_num_G;
for i = 1 : length(kx_vector)
    kx = kx_vector(i); % x coordinate in bloch vector
    ky = ky_vector(i); % y coordinate in bolch vector

    lambda1 = exp(I*kx*max(node_s(:,2))); % periodic boundary factor
    lambda2 = exp(I*ky*max(node_s(:,3))); % periodic boundary factor

    T_F = sparse(length(DOF_A_F),edge_dof_num_F); % inital transfer matrix which product right side 
    T_F = TransferMatrix_F(T_F,...
           index_L_row_F,index_R_row_F,index_B_row_F,index_T_row_F,index_BL_row_F,...
           index_BR_row_F,index_TL_row_F,index_TR_row_F,index_INTER_row_F,...
           index_L_col_F,index_B_col_F,index_BL_col_F,index_INTER_col_F,...
           lambda1,lambda2);

    T_S = sparse(length(DOF_A_S),length(DOF_INDENPENT)); 
    T_S = TransferMatrix_S(T_S,...
           index_L_row,index_R_row,index_B_row,index_T_row,...
           index_BL_row,index_BR_row,index_TR_row,index_TL_row,index_INTER_row,index_L_col,...
           index_B_col,index_BL_col,index_INTER_col,lambda1,lambda2);

    T_G = sparse(length(DOF_A_G),edge_dof_num_G); % inital transfer matrix which product right side 
    T_G = TransferMatrix_G(T_G,...
           index_L_row_G,index_R_row_G,index_B_row_G,index_T_row_G,index_BL_row_G,...
           index_BR_row_G,index_TL_row_G,index_TR_row_G,index_INTER_row_G,...
           index_L_col_G,index_B_col_G,index_BL_col_G,index_INTER_col_G,...
           lambda1,lambda2);

    % modify element stiffness matrix
    % modify mass matrix
%     MA = rho_f*T_S'*R(DOF_A_S,DOF_A_F)*T_F/(T_F'*Kf(DOF_A_F,DOF_A_F)*T_F)*T_F'*R(DOF_A_S,DOF_A_F)'*T_S;
%     MA = sparse(length(DOF_INDENPENT),length(DOF_INDENPENT)); %如果不想加流体，只计算纯固体，拿这一行替代上一行
    
    Ks_re = T_S'*Ks(DOF_A_S,DOF_A_S)*T_S;
    Ms_re = T_S'*Ms(DOF_A_S,DOF_A_S)*T_S;
%     Ms_re = Ms_re+MA;
    R_re = T_S'*R(DOF_A_S,DOF_A_F)*T_F;
    Kf_re = T_F'*Kf(DOF_A_F,DOF_A_F)*T_F;
    Mg_re = T_G'*Mg(DOF_A_G,DOF_A_F)*T_F;
    Kg_re = T_G'*Kg(DOF_A_G,DOF_A_G)*T_G;

    
    % Assembly the global matrix
    K = sparse(dof_s+dof_g,dof_s+dof_g);
    M = sparse(dof_s+dof_g,dof_s+dof_g);
    M11=Ms_re+rho_f*R_re/Kf_re*R_re';
    M12=R_re/Kf_re*Mg_re';
    M22=(1/rho_f)*Mg_re/Kf_re*Mg_re';
    M(1:dof_s,1:dof_s) = M11;
    M(1:dof_s,dof_s+1:dof_s+dof_g) = M12;
    M(dof_s+1:dof_s+dof_g,1:dof_s) = M12';
    M(dof_s+1:dof_s+dof_g,dof_s+1:dof_s+dof_g) = M22;
    
    K(1:dof_s,1:dof_s) = Ks_re;
    K(dof_s+1:dof_s+dof_g,dof_s+1:dof_s+dof_g) = Kg_re;
    
    % calculate frequency of gravity method
    [V_g,lamb_g]=eigs(K,M,dof_g+n,'sm');
    lamb_g = diag(lamb_g);
    freq_g = sqrt(lamb_g)/(2*pi);
    W_gf(:,i) = real(freq_g(1:dof_g));
    W_g(:,i) = real(freq_g(dof_g+1:dof_g+n));

%     % calculate frequency of added mass method
%     [V_a,lamb_a]=eigs(Ks_re,M11,n,'sm');
%     lamb_a = diag(lamb_a);
%     freq_a = sqrt(lamb_a)/(2*pi);
%     W_a(:,i) = real(freq_a);
    
%     % calculate frequency of pure solid
%     [V_s,lamb_s]=eigs(Ks_re,Ms_re,n,'sm');
%     lamb_s = diag(lamb_s);
%     freq_s = sqrt(lamb_s)/(2*pi);
%     W_s(:,i) = real(freq_s);
    disp(['当前验算点为' num2str(i)])
end

%% Gravity wave

k = kx_vector;
freq = sqrt(g.*k.*tanh(k.*hh))/(2*pi);

%% Plot figures
