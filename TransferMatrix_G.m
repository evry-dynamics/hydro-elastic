function T_G = TransferMatrix_G(T_G,...
             index_L_row_G,index_R_row_G,index_B_row_G,index_T_row_G,index_BL_row_G,...
             index_BR_row_G,index_TL_row_G,index_TR_row_G,index_INTER_row_G,...
             index_L_col_G,index_B_col_G,index_BL_col_G,index_INTER_col_G,...
             lambda1,lambda2)

T_G(index_L_row_G,index_L_col_G) =  eye(length(index_L_row_G),length(index_L_row_G));
T_G(index_R_row_G,index_L_col_G) =  lambda1*eye(length(index_R_row_G),length(index_L_row_G));
T_G(index_B_row_G,index_B_col_G) =  eye(length(index_B_row_G),length(index_B_row_G));
T_G(index_T_row_G,index_B_col_G) =  lambda2*eye(length(index_T_row_G),length(index_B_row_G));    
T_G(index_BL_row_G,index_BL_col_G) =eye(length(index_BL_row_G),length(index_BL_row_G));
T_G(index_BR_row_G,index_BL_col_G) =lambda1*eye(length(index_BR_row_G),length(index_BL_row_G));
T_G(index_TR_row_G,index_BL_col_G) =lambda1*lambda2*eye(length(index_TR_row_G),length(index_BL_row_G));
T_G(index_TL_row_G,index_BL_col_G) =lambda2*eye(length(index_TL_row_G),length(index_BL_row_G));
T_G(index_INTER_row_G,index_INTER_col_G) =  eye(length(index_INTER_row_G),length(index_INTER_row_G));



end