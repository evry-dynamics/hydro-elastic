function [N1,N2,N3,N4] = ShapeFunction(s,t)
N1 = 0.25 * (1 - s) * (1 - t);
N2 = 0.25 * (1 + s) * (1 - t);
N3 = 0.25 * (1 + s) * (1 + t);
N4 = 0.25 * (1 - s) * (1 + t);
end