clear
clc
WaveType = 'symmetric';
% -------------------------------------------------------------------------
%                       Define the problem
% -------------------------------------------------------------------------
%% Mesh
ScaleFacor = 0.01; % (20mm)
% Input node and element
% Solid plate
NodeFile = 'NLIST_10_10_0.36.lis';
ElementFile = 'ELIST_10_10_0.36.lis';
[node_s,element_s] = ReadFile_s(NodeFile,ElementFile);
node_s(:,2:4) = node_s(:,2:4)*ScaleFacor;
% Fluid
NodeFile = 'NLIST_f_10_10_3.lis';
ElementFile = 'ELIST_f_10_10_3.lis';
[node_f,element_f] = ReadFile(NodeFile,ElementFile);
node_f(:,2:4) = node_f(:,2:4)*ScaleFacor;
% Fluid depth
node_f(:,4) = node_f(:,4)/6*10;

%% Material
% Basic plate (PVC-Polyvinyl chloride)
Es1 = 6e6; % young's modulus
vs1 = 0.47; % poisson's ratio
rho_s1 = 1.29*10^(3); % density
% Inclusion (wu)
Es2 = 411e9; % young's modulus
vs2 = 0.28; % poisson's ratio
rho_s2 = 19.3*10^(3); % density
% Lame constant for two material
lambda_s1 = Es1*vs1/((1+vs1)*(1-2*vs1)); 
nu_s1 = Es1/(2*(1+vs1)); 
lambda_s2 = Es2*vs2/((1+vs2)*(1-2*vs2)); 
nu_s2 = Es2/(2*(1+vs2)); 
% Fluid domain (water)
rho_f = 1000; % Fluid density
c_f   = 1460; % Fluid acoustic velocity
g = 9.8; % gravitational acceleration (m/s^2)
%% Setting
% Mode order
n = 8;
% Plate thickness and stub height
t = 2*ScaleFacor;
h = 0;
% h = 10*ScaleFacor;
% Imaginary unit
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
% FSI coupling matrix
R = FSIMatrix(node_s,nodeS_FSI,node_f,nodeF_FSI);

%% Degree of freedom
% Solid 
dof_s = 5*length(node_s);
% Fluid
dof_f = length(node_f);
% Free surface
dof_g = length(node_Freesurface);

%% Calculate fluid stiffness and mass matrix
[Kf,Mf] = FEMSolver_Fluid(node_f,element_f,c_f);

%% Calculate free surface stiffness and mass matrix
[Kg,Mg] = FEMSolver_Freesurface(dof_f,dof_g,node_FG,node_Freesurface,element_Freesurface,rho_f,g);

%% Calculate solid stiffness and mass matrix
[Ks,Ms] = FEMSolver_Solid(dof_s,node_s,element_s,nu_s1,rho_s1,vs1,Es1,nu_s2,rho_s2,vs2,Es2,t,h);

%% Calculate Bloch Vector
[kx_vector,ky_vector] = CalculateBlochVector(node_s,WaveType);

%% Bloch periodic conditions
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
    
    Ks_re = T_S'*Ks(DOF_A_S,DOF_A_S)*T_S;
    Ms_re = T_S'*Ms(DOF_A_S,DOF_A_S)*T_S;

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
    
    % Calculate frequency
    [V_g,lamb_g]=eigs(K,M,dof_g+n,'sm');
    lamb_g = diag(lamb_g);
    freq_g = sqrt(lamb_g)/(2*pi);
    W_gf(:,i) = real(freq_g(1:dof_g));
    W_g(:,i) = real(freq_g(dof_g+1:dof_g+n));

    disp(['当前验算点为' num2str(i)])
end

%% Plot figures
ScaleFacor = 1;
xmin = 0; xmax=60;
ymax = 240;
yt = 20;
xt = xmax-xmin;
freq = 0:1:250;

%% BG - Hydro-elastic model
figure
for i = 1:n
    pg=plot([xmin:xmax],W_g(i,:)*ScaleFacor,'LineWidth',1.5); %'k-','*'
    hold on
end
title('Hydro-elastic model','FontName','Times New Roman','FontSize',15)
xlim([xmin xmax])
ylim([0 ymax])
yticks(0:yt:ymax)
xticks([xmin xt/3 xt/3*2 xmax])
xticklabels({'\Gamma ','X','M','\Gamma'})
ylabel('Frequency (Hz)','FontName','Times New Roman','FontSize',15)
gap1_b=max(W_g(6,:));
gap1_t=min(W_g(7,:));
gap2_b=max(W_g(7,:));
gap2_t=min(W_g(8,:));
BG_g = [gap1_b, gap1_t
    gap2_b, gap2_t];
% plot gap
f1 = [1 2 3 4];
% v1 =[xmin gap1_b;xmax gap1_b;xmax gap1_t;xmin gap1_t];
% patch( 'Faces',f1,'vertices',v1,'Facecolor',"#D95319",'FaceAlpha',0.2,'edgecolor','none' );
v2 =[xmin gap2_b;xmax gap2_b;xmax gap2_t;xmin gap2_t];
patch( 'Faces',f1,'vertices',v2,'Facecolor',"#4DBEEE",'FaceAlpha',0.2,'edgecolor','none' );
axis on
box on
grid on