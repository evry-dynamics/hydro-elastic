function [Kf,Mf] = FEMSolver_Fluid(node_f,element_f,c_f)
% ======================== State ==========================================
% this function will use finite element method to calculate global
% stiffness matrix and mass matrix of fluid part.
%========================= Program ======================================== 
    % total degree of freedom
    fdof = 1 * size(node_f,1);       %each node only has one dof in fluid
    % total element number
    te = size(element_f,1);          %total row number of element
    % node number of one element
%     NodeNum = size(element_f,2) - 1; %total column number of element
    % inital K,M
    Kf = sparse(fdof,fdof);
    Mf = sparse(fdof,fdof);
    
    % begin element loop
    for i = 1 : te
        % current loop element
        e = element_f(i,:);
        % eight nodes in current element
        n1 = e(2);
        n2 = e(3);
        n3 = e(4);
        n4 = e(5);
        n5 = e(6);
        n6 = e(7);
        n7 = e(8);
        n8 = e(9);
        % coordinate for current element
        x1 = node_f(n1,2);y1 = node_f(n1,3);z1 = node_f(n1,4);
        x2 = node_f(n2,2);y2 = node_f(n2,3);z2 = node_f(n2,4);
        x3 = node_f(n3,2);y3 = node_f(n3,3);z3 = node_f(n3,4);
        x4 = node_f(n4,2);y4 = node_f(n4,3);z4 = node_f(n4,4);
        x5 = node_f(n5,2);y5 = node_f(n5,3);z5 = node_f(n5,4);
        x6 = node_f(n6,2);y6 = node_f(n6,3);z6 = node_f(n6,4);
        x7 = node_f(n7,2);y7 = node_f(n7,3);z7 = node_f(n7,4);
        x8 = node_f(n8,2);y8 = node_f(n8,3);z8 = node_f(n8,4);

        % element degree of freedom
        edof(1) = n1 ;
        edof(2) = n2 ;
        edof(3) = n3 ;
        edof(4) = n4 ;
        edof(5) = n5 ;
        edof(6) = n6 ;
        edof(7) = n7 ;
        edof(8) = n8 ;
        edof = edof(:);%这一步将矩阵转化为一列

        [k,m] = ElementStiffness_F(x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8, ...
                                    z1,z2,z3,z4,z5,z6,z7,z8,c_f);
        Kf(edof,edof) = Kf(edof,edof) + k;
        Mf(edof,edof) = Mf(edof,edof) + m;
                    
        clear edof
    end
end