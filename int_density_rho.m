clear all

syms z h t rho
R = [1 0 0 -z 0;
    0 1 0 0 -z;
    0 0 1 0 0];


rho_gen = R'*rho*R;
% plate
rho_p = int(rho_gen, z, [-t/2, t/2])
% plate + stub
rho_h = int(rho_gen, z, [-t/2, t/2+h])


R = [0 0 1 0 0;
    1 0 0 -z 0;
    0 1 0 0 -z];


rho_gen = R'*rho*R;
% plate
rho_p = int(rho_gen, z, [-t/2, t/2])
% plate + stub
rho_h = int(rho_gen, z, [-t/2, t/2+h])