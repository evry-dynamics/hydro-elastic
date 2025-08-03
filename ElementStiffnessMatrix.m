function k = ElementStiffnessMatrix(x1,x2,x3,x4,y1,y2,y3,y4,D,Type)

% full integral
p = 1/sqrt(3);
GaussPoint = [-p,-p;-p,p;p,-p;p,p];
k = zeros(20,20);
for i = 1:size(GaussPoint,1)
    J = Jacobian(GaussPoint(i,1),GaussPoint(i,2),x1,x2,x3,x4,y1,y2,y3,y4);
    J_det = det(J);
    
    B = StrainMatrix(GaussPoint(i,1),GaussPoint(i,2),x1,x2,x3,x4,y1,y2,y3,y4,Type);
    k = k + B' * D * B * J_det;
end
end