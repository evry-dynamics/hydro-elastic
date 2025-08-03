function G_F = StrainMatrix_F(s,t,w,J)

[dN1ds,dN1dt,dN1dw,dN2ds,dN2dt,dN2dw,dN3ds,dN3dt,dN3dw,dN4ds,dN4dt,dN4dw, ...
 dN5ds,dN5dt,dN5dw,dN6ds,dN6dt,dN6dw,dN7ds,dN7dt,dN7dw,dN8ds,dN8dt,dN8dw] = First_3D(s,t,w);

dNdstw=[dN1ds,dN2ds,dN3ds,dN4ds,dN5ds,dN6ds,dN7ds,dN8ds;
        dN1dt,dN2dt,dN3dt,dN4dt,dN5dt,dN6dt,dN7dt,dN8dt;
        dN1dw,dN2dw,dN3dw,dN4dw,dN5dw,dN6dw,dN7dw,dN8dw;];

G_F=inv(J)*dNdstw;

end