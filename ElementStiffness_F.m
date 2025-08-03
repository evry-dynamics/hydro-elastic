function [k,m] = ElementStiffness_F(x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8, ...
                                    z1,z2,z3,z4,z5,z6,z7,z8,c_f)
% ================================ State ==================================
% input:
% x1,x2,x3,y1,y2,y3: node coor
% C,A,ea,es,Ea,Es: generilized contitutive 
% rho: density matrix
% ga,gd: feedback control 
% output:
% k: element stiffness matrix
% ================================ Program ================================
        k = zeros(8,8);
        m = zeros(8,8);
        p = 1/sqrt(3);
        GassPoint = [-p,-p,-p;
                     -p,-p, p;
                     -p, p,-p;
                     -p, p, p;
                      p,-p,-p;
                      p,-p, p;
                      p, p,-p;
                      p, p, p];
        for i = 1 : size(GassPoint,1)
            s = GassPoint(i,1);
            t = GassPoint(i,2);
            w = GassPoint(i,3);
            [J,J_det] = JacobianMatrix_3D(s,t,w,x1,x2,x3,x4,x5,x6,x7,x8, ...
                                       y1,y2,y3,y4,y5,y6,y7,y8, ...
                                       z1,z2,z3,z4,z5,z6,z7,z8);
            %J_det = det(J);
%           Bb = StrainMatrix(s,t,x1,x2,x3,x4,y1,y2,y3,y4,J,'classical');
            G_F = StrainMatrix_F(s,t,w,J);
            N = ShapeFunction_3D(s,t,w);

            k = k + G_F' * G_F* J_det ;
            m = m + 1/c_f^2*(N'*N)*J_det;
        end
        
        
end