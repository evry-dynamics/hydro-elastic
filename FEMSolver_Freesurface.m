function [K,M] = FEMSolver_Freesurface(dof_f,dof_g,node_FG,node_Freesurface,element_Freesurface,rho_f,g)

K = sparse(dof_g,dof_g); % inital global stiffness matrix
M = sparse(dof_g,dof_f); % inital global mass matrix
for i = 1:size(element_Freesurface,1)
    x1 = node_Freesurface(element_Freesurface(i,2),2);
    x2 = node_Freesurface(element_Freesurface(i,3),2);
    x3 = node_Freesurface(element_Freesurface(i,4),2);
    x4 = node_Freesurface(element_Freesurface(i,5),2);

    y1 = node_Freesurface(element_Freesurface(i,2),3);
    y2 = node_Freesurface(element_Freesurface(i,3),3);
    y3 = node_Freesurface(element_Freesurface(i,4),3);
    y4 = node_Freesurface(element_Freesurface(i,5),3);

    % degree of freedom of element vector(index of assemble global matrix)
    % 节点对应自由液面的自由度编号
    edof_g = [element_Freesurface(i,2:5)];
    % 节点对应整体流域的自由度编号
    edof_f = [node_FG(element_Freesurface(i,2:5))]';

    p = 1/sqrt(3);
    GaussPoint = [-p,-p;-p,p;p,-p;p,p];
    k = zeros(4,4);
    m = zeros(4,4);
    for ig = 1:size(GaussPoint,1)
        s = GaussPoint(ig,1);
        t = GaussPoint(ig,2);
        J = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
        J_det = det(J);
        [N1,N2,N3,N4] = ShapeFunction(s,t);
        N = [N1,N2,N3,N4];
        k = k + (rho_f * g) * (N' * N) * J_det;
        m = m + rho_f * (N' * N) * J_det;
    end

    K(edof_g,edof_g) = K(edof_g,edof_g) + k;
    M(edof_g,edof_f) = M(edof_g,edof_f) + m;
end

