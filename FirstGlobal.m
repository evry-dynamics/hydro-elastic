function [dN1dx,dN1dy,dN2dx,dN2dy,dN3dx,dN3dy,dN4dx,dN4dy] = FirstGlobal(s,t,x1,x2,x3,x4,y1,y2,y3,y4)
J = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
J_inv = inv(J);

dN1ds = t/4 - 1/4;
dN1dt = s/4 - 1/4;

dN2ds = 1/4 - t/4;
dN2dt = - s/4 - 1/4;

dN3ds = t/4 + 1/4;
dN3dt = s/4 + 1/4;

dN4ds = - t/4 - 1/4;
dN4dt = 1/4 - s/4;

dN1 = J_inv * [dN1ds;dN1dt];
dN2 = J_inv * [dN2ds;dN2dt];
dN3 = J_inv * [dN3ds;dN3dt];
dN4 = J_inv * [dN4ds;dN4dt];

dN1dx = dN1(1);
dN1dy = dN1(2);
dN2dx = dN2(1);
dN2dy = dN2(2);
dN3dx = dN3(1);
dN3dy = dN3(2);
dN4dx = dN4(1);
dN4dy = dN4(2);
end
