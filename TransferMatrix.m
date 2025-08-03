function T = TransferMatrix(T,DOF_L,DOF_R,DOF_B,DOF_T,DOF_BL,DOF_BR,DOF_TR,DOF_TL,index_L,index_B,index_BL,lambda1,lambda2)
% ============================================= %
% This function calculate transfer matrix to modify stiffness matrix and
% mass matrix 
% refer:Krattiger D , Hussein M I . Generalized Bloch mode 
% synthesis for accelerated calculation of elastic band structures[J]. Journal of Computational Physics, 2018, 357:183-205.
% 
% param T: sparse matrix which filled all zeros
% param DOF_L: degree of freedom of left boundary 
% param DOF_R: degree of freedom of right boundary
% param DOF_B: degree of freedom of bottom boundary
% param DOF_T: degree of freedom of top boundary
% param DOF_BL: degree of freedom of lower-left boundary 
% param DOF_BR: degree of freedom of lower-right boundary
% param DOF_TR: degree of freedom of top-right boundary
% param DOF_TL: degree of freedom of top-left boundary
% param index_L: column index of left boundary in matrix
% param index_B: column index of right boundary in matrix
% param index_BL: column index of lower-right boundary in matrix
% param lambda1: periodic boundary factor according to x direction
% param lambda2: periodic boundary factor according to y direction
% 
% return T: Transfer matrix
% ============================================= %
    T(DOF_L,index_L) =  eye(length(DOF_L),length(DOF_L));
    T(DOF_R,index_L) =  lambda1*eye(length(DOF_R),length(DOF_L));
    T(DOF_B,index_B) =  eye(length(DOF_B),length(DOF_B));
    T(DOF_T,index_B) =  lambda2 *eye(length(DOF_T),length(DOF_B));
    T(DOF_BL,index_BL) =  eye(length(DOF_BL),length(DOF_BL));
    T(DOF_BR,index_BL) =  lambda1*eye(length(DOF_BR),length(DOF_BL));
    T(DOF_TR,index_BL) =  lambda1*lambda2 *eye(length(DOF_TR),length(DOF_BL));
    T(DOF_TL,index_BL) =  lambda2 *eye(length(DOF_TL),length(DOF_BL));
end