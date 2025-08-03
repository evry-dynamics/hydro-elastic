function B = StrainMatrix(s,t,x1,x2,x3,x4,y1,y2,y3,y4,type)
[J] = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
[N1,N2,N3,N4] = ShapeFunction(s,t);
[dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] = FirstGlobal(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
    switch type
        case 'Bend'            
            B = [dN1dx,     0,     0,     0,     0, dN2dx,     0,     0,     0,     0, dN3dx,     0,     0,     0,     0, dN4dx,     0,     0,     0,     0
                     0,     0,     0, dN1dx,     0,     0,     0,     0, dN2dx,     0,     0,     0,     0, dN3dx,     0,     0,     0,     0, dN4dx,     0
                     0, dN1dy,     0,     0,     0,     0, dN2dy,     0,     0,     0,     0, dN3dy,     0,     0,     0,     0, dN4dy,     0,     0,     0
                     0,     0,     0,     0, dN1dy,     0,     0,     0,     0, dN2dy,     0,     0,     0,     0, dN3dy,     0,     0,     0,     0, dN4dy
                 dN1dy, dN1dx,     0,     0,     0, dN2dy, dN2dx,     0,     0,     0, dN3dy, dN3dx,     0,     0,     0, dN4dy, dN4dx,     0,     0,     0
                     0,     0,     0, dN1dy, dN1dx,     0,     0,     0, dN2dy, dN2dx,     0,     0,     0, dN3dy, dN3dx,     0,     0,     0, dN4dy, dN4dx];
        case 'Shear'
            [J1] = Jacobian(0,-1,x1,x2,x3,x4,y1,y2,y3,y4);
            [J2] = Jacobian(1,0,x1,x2,x3,x4,y1,y2,y3,y4);
            [J3] = Jacobian(0,1,x1,x2,x3,x4,y1,y2,y3,y4);
            [J4] = Jacobian(-1,0,x1,x2,x3,x4,y1,y2,y3,y4);
            C = zeros(8,8);
            C(1:2,1:2) = J1;
            C(3:4,3:4) = J2;
            C(5:6,5:6) = J3;
            C(7:8,7:8) = J4;
            [N1,N2,N3,N4] = ShapeFunction(0,-1);
            [dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] =...
                    FirstGlobal(0,-1,x1,x2,x3,x4,y1,y2,y3,y4);
            Bs1 = [    0,     0, dN1dx,   -N1,     0,     0,     0, dN2dx,   -N2,     0,     0,     0, dN3dx,   -N3,     0,     0,     0, dN4dx,   -N4,     0
                     0,     0, dN1dy,     0,   -N1,     0,     0, dN2dy,     0,   -N2,     0,     0, dN3dy,     0,   -N3,     0,     0, dN4dy,     0,   -N4];

            [N1,N2,N3,N4] = ShapeFunction(1,0);
            [dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] =...
                FirstGlobal(1,0,x1,x2,x3,x4,y1,y2,y3,y4);
            Bs2 = [    0,     0, dN1dx,   -N1,     0,     0,     0, dN2dx,   -N2,     0,     0,     0, dN3dx,   -N3,     0,     0,     0, dN4dx,   -N4,     0
                     0,     0, dN1dy,     0,   -N1,     0,     0, dN2dy,     0,   -N2,     0,     0, dN3dy,     0,   -N3,     0,     0, dN4dy,     0,   -N4];


            [N1,N2,N3,N4] = ShapeFunction(0,1);
            [dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] =...
                FirstGlobal(0,1,x1,x2,x3,x4,y1,y2,y3,y4);
            Bs3 = [    0,     0, dN1dx,   -N1,     0,     0,     0, dN2dx,   -N2,     0,     0,     0, dN3dx,   -N3,     0,     0,     0, dN4dx,   -N4,     0
                     0,     0, dN1dy,     0,   -N1,     0,     0, dN2dy,     0,   -N2,     0,     0, dN3dy,     0,   -N3,     0,     0, dN4dy,     0,   -N4];

            [N1,N2,N3,N4] = ShapeFunction(-1,0);
            [dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] =...
                FirstGlobal(-1,0,x1,x2,x3,x4,y1,y2,y3,y4);
            Bs4 = [    0,     0, dN1dx,   -N1,     0,     0,     0, dN2dx,   -N2,     0,     0,     0, dN3dx,   -N3,     0,     0,     0, dN4dx,   -N4,     0
                     0,     0, dN1dy,     0,   -N1,     0,     0, dN2dy,     0,   -N2,     0,     0, dN3dy,     0,   -N3,     0,     0, dN4dy,     0,   -N4];

            Bs = [Bs1;Bs2;Bs3;Bs4];
            APT = 0.5 * [(1-t),0,0,0,(1+t),0,0,0;
                    0,0,0,(1+s),0,0,0,(1-s)];
            B = inv(J) * APT * C * Bs;
%             without shear correction
%             B = [    0,     0, dN1dx,   -N1,     0,     0,     0, dN2dx,   -N2,     0,     0,     0, dN3dx,   -N3,     0,     0,     0, dN4dx,   -N4,     0
%                      0,     0, dN1dy,     0,   -N1,     0,     0, dN2dy,     0,   -N2,     0,     0, dN3dy,     0,   -N3,     0,     0, dN4dy,     0,   -N4];
    end
end