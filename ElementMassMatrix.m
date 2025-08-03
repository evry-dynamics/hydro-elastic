function  m = ElementMassMatrix(x1,x2,x3,x4,y1,y2,y3,y4,rho,t,h,thick)
switch thick
    case 'plate'
        rho_maxtrix = [rho*t,     0,     0,            0,            0
                           0, rho*t,     0,            0,            0
                           0,     0, rho*t,            0,            0
                           0,     0,     0, (rho*t^3)/12,            0
                           0,     0,     0,            0, (rho*t^3)/12];
    case 'stub'
%       不对称结构-stub位于plate上方
        rho_maxtrix = [       rho*(h + t),                  0,           0,                 -(h*rho*(h + t))/2,                                  0
                                        0,        rho*(h + t),           0,                                  0,                 -(h*rho*(h + t))/2
                                        0,                  0, rho*(h + t),                                  0,                                  0
                       -(h*rho*(h + t))/2,                  0,           0, (rho*t^3)/24 + (rho*(h + t/2)^3)/3,                                  0
                                        0, -(h*rho*(h + t))/2,           0,                                  0, (rho*t^3)/24 + (rho*(h + t/2)^3)/3];
end

p = 1/sqrt(3);
GaussPoint = [-p,-p;-p,p;p,-p;p,p];
m=zeros(20,20);
for i = 1:size(GaussPoint,1)
    s = GaussPoint(i,1);
    t = GaussPoint(i,2);
    
    J = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
    J_det = det(J);
    
    [N1,N2,N3,N4] = ShapeFunction(s,t);
    
    N = [N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0,0,0;
        0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0,0;
        0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0;
        0,0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0;
        0,0,0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4];
    
    m = m + N' * rho_maxtrix * N * J_det;
end
% m_sum = sum(m,'all');
% m_diag = diag(m,0);
% m_diag_sum = sum(m_diag);
% aplha = m_sum/m_diag_sum;
% m = diag(aplha * m_diag);

end