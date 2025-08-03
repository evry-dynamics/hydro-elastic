function [a,b,node,element,DOF_INTER,DOF_EDGE_SORT,DOF_INDENPENT,DOF_A,edge_dof_num,index_L_row,index_B_row,index_BL_row,index_R_row,index_T_row,...
    index_BR_row,index_TL_row,index_TR_row,index_L_col,index_B_col,index_BL_col] = SortNode(node,element)
% ========================================== %
% This function classifes node according to the lower left, 
% lower, upper, upper left, left, right, upper right and lower right boundaries
% I know this code looks confusing,but the function is only used classify 
% node information,which is convenient for me to impose periodic boundary 
% condition.The specific construction of the code is not important. 
%
% param node: node matrix 
% param element: element matrix 
% 
% return node: node matrix after classify 
% return element: element matrix after classify 
% return DOF_INDENPENT: independent degree of freedom vector
% return DOF_A: degree of freedom in edge
% return edge_dof_num: the number of degree of freedom in edge
% return index_L_row: row index of left boundary in matrix
% return index_B_row: row index of bottom boudary in matrix
% return index_BL_row: row index of lower left boundary in matrix
% return index_R_row: row index of right boundary in matrix
% return index_T_row: row index of top boundary in matix 
% return index_BR_row: row index of lower right boundary in matrix
% return index_TL_row: row index of top left boundary in matrix 
% return index_TR_row: row index of top right boundary in matrix 
% return index_L_col: column index of left boundary in matix
% return index_B_col: column index of bottom boundary in matrix
% return index_BL_col: column index of lower left boundary in matrix
% ========================================== %

co = node(:,2:3);
xmax = max(co(:,1));  ymin = min(co(:,2));
xmin = min(co(:,1));  ymax = max(co(:,2)); 

a = abs(xmax - xmin);
b = abs(ymax - ymin);

id_left = find( co(:,1)==xmin );
[out,idx] = sort( co(id_left,2) );
id_left =id_left(idx);
node_left = id_left(2:end-1);
DOF_L = reshape([5*node_left-4 5*node_left-3 5*node_left-2 5*node_left-1 5*node_left]',5*length(node_left),1);

id_BL = id_left(1);
DOF_BL = [id_BL*5-4;id_BL*5-3;id_BL*5-2;id_BL*5-1;id_BL*5];
id_TL = id_left(end);
DOF_TL = [id_TL*5-4;id_TL*5-3;id_TL*5-2;id_TL*5-1;id_TL*5];


id_right = find( co(:,1)==xmax );
[out,idx] = sort( co(id_right,2) );
id_right =id_right(idx);
node_right = id_right(2:end-1);
DOF_R = reshape([5*node_right-4 5*node_right-3 5*node_right-2 5*node_right-1 5*node_right]',5*length(node_right),1);

id_BR = id_right(1);
DOF_BR = [id_BR*5-4;id_BR*5-3;id_BR*5-2;id_BR*5-1;id_BR*5];
id_TR = id_right(end);
DOF_TR = [id_TR*5-4;id_TR*5-3;id_TR*5-2;id_TR*5-1;id_TR*5];


id_upper = find( co(:,2)==ymax );
[out,idx] = sort( co(id_upper,1) );
id_upper =id_upper(idx);
node_upper = id_upper(2:end-1);
DOF_T = reshape([5*node_upper-4 5*node_upper-3 5*node_upper-2 5*node_upper-1 5*node_upper]',5*length(node_upper),1);

id_lower = find( co(:,2)==ymin );
[out,idx] = sort( co(id_lower,1) );
id_lower =id_lower(idx);
node_lower = id_lower(2:end-1);
DOF_B = reshape([5*node_lower-4 5*node_lower-3 5*node_lower-2 5*node_lower-1 5*node_lower]',5*length(node_lower),1);

% 所有自由度向量
DOF_ALL = reshape([5*node(:,1)-4 5*node(:,1)-3 5*node(:,1)-2 5*node(:,1)-1 5*node(:,1)]',5*length(node(:,1)),1);
% 内部节点自由度向量
DOF_INTER = setdiff(DOF_ALL,[DOF_BL;DOF_TR;DOF_TL;DOF_BR;DOF_B;DOF_L;DOF_T;DOF_R]);
%% 形成P矩阵
% 独立节点自由度向量
DOF_EDGE = [DOF_L;DOF_B;DOF_BL];
DOF_EDGE_SORT = DOF_EDGE ;
edge_dof_num = length(DOF_EDGE); % 边节点个数
% 独立节点自由度
DOF_INDENPENT = [DOF_L;DOF_B;DOF_BL;DOF_INTER];

% 总边界节点自由度
DOF_A = [DOF_L;DOF_R;DOF_B;DOF_T;DOF_BL;DOF_BR;DOF_TR;DOF_TL];
% DOF_A = sort(DOF_A); % 自由度重排
% 寻找总边界节点自由度在转化矩阵中的行索引位置
for i = 1:length(DOF_L)
    index_L_row(i) = find(DOF_A==DOF_L(i));
end
for i = 1:length(DOF_B)
    index_B_row(i) = find(DOF_A==DOF_B(i));
end
for i = 1:length(DOF_BL)
    index_BL_row(i) = find(DOF_A==DOF_BL(i));
end
for i = 1:length(DOF_R)
    index_R_row(i) = find(DOF_A==DOF_R(i));
end
for i = 1:length(DOF_T)
    index_T_row(i) = find(DOF_A==DOF_T(i));
end
for i = 1:length(DOF_BR)
    index_BR_row(i) = find(DOF_A==DOF_BR(i));
end
for i = 1:length(DOF_TL)
    index_TL_row(i) = find(DOF_A==DOF_TL(i));
end
for i = 1:length(DOF_TR)
    index_TR_row(i) = find(DOF_A==DOF_TR(i));
end
% 寻找边节点自由度在转化矩阵中的列索引位置
for i = 1:length(DOF_L)
    index_L_col(i) = find(DOF_EDGE_SORT==DOF_L(i));
end
for i = 1:length(DOF_B)
    index_B_col(i) = find(DOF_EDGE_SORT==DOF_B(i));
end
for i = 1:length(DOF_BL)
    index_BL_col(i) = find(DOF_EDGE_SORT==DOF_BL(i));
end
end