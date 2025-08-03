function [dN1ds,dN1dt,dN1dw, ...
          dN2ds,dN2dt,dN2dw, ...
          dN3ds,dN3dt,dN3dw, ...
          dN4ds,dN4dt,dN4dw, ...
          dN5ds,dN5dt,dN5dw, ...
          dN6ds,dN6dt,dN6dw, ...
          dN7ds,dN7dt,dN7dw, ...
          dN8ds,dN8dt,dN8dw] = First_3D(s,t,w)

dN1ds = -1/8*(1-t)*(1-w);
dN2ds =  1/8*(1-t)*(1-w);
dN3ds =  1/8*(1+t)*(1-w);
dN4ds = -1/8*(1+t)*(1-w);
dN5ds = -1/8*(1-t)*(1+w);
dN6ds =  1/8*(1-t)*(1+w);
dN7ds =  1/8*(1+t)*(1+w);
dN8ds = -1/8*(1+t)*(1+w);

dN1dt = -1/8*(1-s)*(1-w);
dN2dt = -1/8*(1+s)*(1-w);
dN3dt =  1/8*(1+s)*(1-w);
dN4dt =  1/8*(1-s)*(1-w);
dN5dt = -1/8*(1-s)*(1+w);
dN6dt = -1/8*(1+s)*(1+w);
dN7dt =  1/8*(1+s)*(1+w);
dN8dt =  1/8*(1-s)*(1+w);

dN1dw = -1/8*(1-s)*(1-t);
dN2dw = -1/8*(1+s)*(1-t);
dN3dw = -1/8*(1+s)*(1+t);
dN4dw = -1/8*(1-s)*(1+t);
dN5dw =  1/8*(1-s)*(1-t);
dN6dw =  1/8*(1+s)*(1-t);
dN7dw =  1/8*(1+s)*(1+t);
dN8dw =  1/8*(1-s)*(1+t);



end