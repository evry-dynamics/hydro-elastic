function [Ks,Ms] = FEMSolver_Solid(dof_s,node_s,element_s,nu_s1,rho_s1,vs1,Es1,nu_s2,rho_s2,vs2,Es2,t,h)
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