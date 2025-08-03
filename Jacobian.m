function [J] = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4)
dN1ds = t/4 - 1/4;
dN1dt = s/4 - 1/4;

dN2ds = 1/4 - t/4;
dN2dt = - s/4 - 1/4;

dN3ds = t/4 + 1/4;
dN3dt = s/4 + 1/4;

dN4ds = - t/4 - 1/4;
dN4dt = 1/4 - s/4;

dxds = dN1ds * x1 + dN2ds * x2 + dN3ds * x3 + dN4ds * x4;
dxdt = dN1dt * x1 + dN2dt * x2 + dN3dt * x3 + dN4dt * x4;

dyds = dN1ds * y1 + dN2ds * y2 + dN3ds * y3 + dN4ds * y4;
dydt = dN1dt * y1 + dN2dt * y2 + dN3dt * y3 + dN4dt * y4;

J = [dxds,dyds;dxdt,dydt];

end