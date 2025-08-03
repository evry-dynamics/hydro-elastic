function T_F = TransferMatrix_F(T_F,...
             index_L_row_F,index_R_row_F,index_B_row_F,index_T_row_F,index_BL_row_F,...
             index_BR_row_F,index_TL_row_F,index_TR_row_F,index_INTER_row_F,...
             index_L_col_F,index_B_col_F,index_BL_col_F,index_INTER_col_F,...
             lambda1,lambda2)

T_F(index_L_row_F,index_L_col_F) =  eye(length(index_L_row_F),length(index_L_row_F));
T_F(index_R_row_F,index_L_col_F) =  lambda1*eye(length(index_R_row_F),length(index_L_row_F));
T_F(index_B_row_F,index_B_col_F) =  eye(length(index_B_row_F),length(index_B_row_F));
T_F(index_T_row_F,index_B_col_F) =  lambda2*eye(length(index_T_row_F),length(index_B_row_F));    
T_F(index_BL_row_F,index_BL_col_F) =eye(length(index_BL_row_F),length(index_BL_row_F));
T_F(index_BR_row_F,index_BL_col_F) =lambda1*eye(length(index_BR_row_F),length(index_BL_row_F));
T_F(index_TR_row_F,index_BL_col_F) =lambda1*lambda2*eye(length(index_TR_row_F),length(index_BL_row_F));
T_F(index_TL_row_F,index_BL_col_F) =lambda2*eye(length(index_TL_row_F),length(index_BL_row_F));
T_F(index_INTER_row_F,index_INTER_col_F) =  eye(length(index_INTER_row_F),length(index_INTER_row_F));



end