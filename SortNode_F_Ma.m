function [nodeF_FSI,node_Freesurface,element_Freesurface,node_FG,DOF_INTER,DOF_EDGE_SORT,DOF_INDENPENT,DOF_A,...
edge_dof_num,index_L_row,index_B_row,index_BL_row,index_R_row,index_T_row,...
index_BR_row,index_TL_row,index_TR_row,index_INTER_row,...
index_L_col,index_B_col,index_BL_col,index_INTER_col] = SortNode_F_Ma(node_f)



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

co = node_f(:,2:4);
xmax = max(co(:,1));  ymin = min(co(:,2));  zmin = min(co(:,3));
xmin = min(co(:,1));  ymax = max(co(:,2));  zmax = max(co(:,3));


a = abs(xmax - xmin);
b = abs(ymax - ymin);
c = abs(zmax - zmin);

%% 找出Free surface的点及自由度并按x&y坐标从小到大排序
id_Freesurface = find( co(:,3)==zmax );   %找到Freesurface边界上的节点编号
[~,idx] = sort( co(id_Freesurface,1) );   %对Freesurface边界上的节点按x进行从小到大排序
id_Freesurface = id_Freesurface(idx);     %排序后编号
[~,idx] = sort( co(id_Freesurface,2) );   %对Freesurface边界上的节点按y进行从小到大排序
id_Freesurface = id_Freesurface(idx);     %排序后编号
node_Freesurface = [(1:length(id_Freesurface))',co(id_Freesurface,1),co(id_Freesurface,2)];        %Freesurface边界的节点编号
node_FG = id_Freesurface; % 自由液面的节点对应于整体流域的节点编号


%% 找出BottomLeft的点及自由度并按z坐标从小到大排序
id_BL = find( co(:,1)==xmin & co(:,2)==ymin  );%找到BottomLeft边界上的节点编号
[~,idx] = sort( co(id_BL,3) );%对BottomLeft边界上的节点进行从小到大排序
id_BL =id_BL(idx);           %排序后编号
node_BL = id_BL;           
DOF_BL = setdiff(node_BL,node_FG);           %BottomLeft边界点的自由度

NZ=length(DOF_BL);

id_NX = find( co(:,2)==xmin & co(:,3)==zmin);%找到BottomLeft边界上的节点编号
NX=length(id_NX);
% NX=b/c*10*(NZ-1)+1; % 流体高度放大10倍
% NX=21; % 流体高度100mm，缩小10倍

%% 找出BottomRight的点及自由度并按z坐标从小到大排序
id_BR = find( co(:,1)==xmax & co(:,2)==ymin  );%找到BottomrRight边界上的节点编号
[~,idx] = sort( co(id_BR,3) );%对BottomrRight边界上的节点进行从小到大排序
id_BR =id_BR(idx);           %排序后编号
node_BR = id_BR;   
DOF_BR = setdiff(node_BR,node_FG);           %BottomrRight边界点的自由度

%% 找出TopRight的点及自由度并按z坐标从小到大排序
id_TR = find( co(:,1)==xmax & co(:,2)==ymax  );%找到TopRight边界上的节点编号
[~,idx] = sort( co(id_TR,3) );%对TopRight边界上的节点进行从小到大排序
id_TR =id_TR(idx);           %排序后编号
node_TR = id_TR;   
DOF_TR = setdiff(node_TR,node_FG);           %TopRight边界点的自由度

%% 找出TopLeft的点及自由度并按z坐标从小到大排序
id_TL = find( co(:,1)==xmin & co(:,2)==ymax  );%找到TopLeft边界上的节点编号
[~,idx] = sort( co(id_TL,3) );%对TopLeft边界上的节点进行从小到大排序
id_TL =id_TL(idx);           %排序后编号
node_TL = id_TL;   
DOF_TL = setdiff(node_TL,node_FG);           %TopLeft边界点的自由度

%% 找出Left的点及自由度并按z&y坐标从小到大排序
id_L = find( co(:,1)==xmin );%找到Left边界上的节点编号
[~,idx] = sort( co(id_L,3) );%对Left边界上的节点按z进行从小到大排序
id_L = id_L(idx);             %排序后编号
[~,idx] = sort( co(id_L,2) );%对Left边界上的节点按y进行从小到大排序
id_L = id_L(idx);             %排序后编号
node_L = id_L;   
node_L = setdiff(node_L,[node_L(1:1:NZ);node_L(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_L,3) );%对Left边界上的节点按z进行从小到大排序
node_L = node_L(idx);             %排序后编号
[~,idx] = sort( co(node_L,2) );%对Left边界上的节点按y进行从小到大排序
node_L = node_L(idx);             %排序后编号
DOF_L = setdiff(node_L,node_FG);                   %Left边界点的自由度

%% 找出Right的点及自由度并按z&y坐标从小到大排序
id_R = find( co(:,1)==xmax );%找到Right边界上的节点编号
[~,idx] = sort( co(id_R,3) );%对Right边界上的节点按z进行从小到大排序
id_R = id_R(idx);             %排序后编号
[~,idx] = sort( co(id_R,2) );%对Right边界上的节点按y进行从小到大排序
id_R = id_R(idx);             %排序后编号
node_R = id_R;   
node_R = setdiff(node_R,[node_R(1:1:NZ);node_R(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_R,3) );%对Right边界上的节点按z进行从小到大排序
node_R = node_R(idx);             %排序后编号
[~,idx] = sort( co(node_R,2) );%对Right边界上的节点按y进行从小到大排序
node_R = node_R(idx);             %排序后编号
DOF_R = setdiff(node_R,node_FG);                   %Right边界点的自由度

%% 找出Bottom的点及自由度并按z&x坐标从小到大排序
id_B = find( co(:,2)==ymin );%找到Bottom边界上的节点编号
[~,idx] = sort( co(id_B,3) );%对Bottom边界上的节点按z进行从小到大排序
id_B = id_B(idx);             %排序后编号
[~,idx] = sort( co(id_B,1) );%对Bottom边界上的节点按x进行从小到大排序
id_B = id_B(idx);             %排序后编号
node_B = id_B;   
node_B = setdiff(node_B,[node_B(1:1:NZ);node_B(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_B,3) );%对Bottom边界上的节点按z进行从小到大排序
node_B = node_B(idx);             %排序后编号
[~,idx] = sort( co(node_B,1) );%对Bottom边界上的节点按x进行从小到大排序
node_B = node_B(idx);             %排序后编号
DOF_B = setdiff(node_B,node_FG);                   %Bottom边界点的自由度

%% 找出Top的点及自由度并按z&x坐标从小到大排序
id_T = find( co(:,2)==ymax );%找到Top边界上的节点编号
[~,idx] = sort( co(id_T,3) );%对Top边界上的节点按z进行从小到大排序
id_T = id_T(idx);             %排序后编号
[~,idx] = sort( co(id_T,1) );%对Top边界上的节点按x进行从小到大排序
id_T = id_T(idx);             %排序后编号
node_T = id_T;   
node_T = setdiff(node_T,[node_T(1:1:NZ);node_T(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_T,3) );%对Top边界上的节点按z进行从小到大排序
node_T = node_T(idx);             %排序后编号
[~,idx] = sort( co(node_T,1) );%对Top边界上的节点按x进行从小到大排序
node_T = node_T(idx);             %排序后编号
DOF_T = setdiff(node_T,node_FG);                   %Top边界点的自由度                                 

%% 找出FSI的点及自由度并按x&y坐标从小到大排序
% id_FSI = find( co(:,3)==zmax );   %找到FSI边界上的节点编号 %流体位于固体下方
id_FSI = find( co(:,3)==zmin );   %找到FSI边界上的节点编号 %流体位于固体上方
[~,idx] = sort( co(id_FSI,1) );   %对FSI边界上的节点按x进行从小到大排序
id_FSI = id_FSI(idx);             %排序后编号
[~,idx] = sort( co(id_FSI,2) );   %对FSI边界上的节点按y进行从小到大排序
id_FSI = id_FSI(idx);             %排序后编号
nodeF_FSI = id_FSI;               %FSI边界的节点编号
              

%% 对Free surface边界上的节点按照逆时针顺序进行排列
nx=size(node_Freesurface,1)^0.5;ny=size(node_Freesurface,1)^0.5; %%只适用于x和y方向单元个数相等的情况 Free surface边界上x和y方向上的节点数。nx=ny=21
serialG = reshape(node_Freesurface(:,1),nx,ny); % 行为按y方向从小到大排列的点，列为按x方向从小到大排列的点。
% serialG = serialG'; %转置后，行为按x方向从小到大排列的点，列为按y方向从小到大排列的点。

n=size(node_Freesurface,1)^0.5-1;  % Free surface边界上的每行每列的单元数
element_Freesurface=zeros(n^2,5); % 初始化单元矩阵
a=0; % 重构Freesurface边界上的单元时的计数参数，此时初始化。
for i=1:n
    for j=1:n
        a=a+1;
        %每个Freesurface边界单元的节点编号为逆时针排列
        element_Freesurface(a,2) = serialG(i,j);
        element_Freesurface(a,3) = serialG(i+1,j);
        element_Freesurface(a,4) = serialG(i+1,j+1);
        element_Freesurface(a,5) = serialG(i,j+1);
    end
end
element_Freesurface(:,1)=(1:a)';


% 所有自由度向量
DOF_ALL   = setdiff(node_f(:,1),node_FG); %除掉自由液面自由度
% 内部节点自由度向量
DOF_INTER = setdiff(DOF_ALL,[DOF_BL;DOF_TR;DOF_TL;DOF_BR;DOF_B;DOF_L;DOF_T;DOF_R]);

%% 形成P矩阵
% 独立节点自由度向量
DOF_EDGE       = [DOF_L;DOF_B;DOF_BL];
DOF_EDGE_SORT  = DOF_EDGE ;
%edge_dof_num   = length(DOF_EDGE); % 边节点自由度个数
% 独立节点自由度
DOF_INDENPENT  = [DOF_L;DOF_B;DOF_BL;DOF_INTER];
edge_dof_num   = length(DOF_INDENPENT);

% 总边界节点自由度
DOF_A  = [DOF_L;DOF_R;DOF_B;DOF_T;DOF_BL;DOF_BR;DOF_TR;DOF_TL;DOF_INTER];
% DOF_A = sort(DOF_A); % 自由度重排
% 寻找总边界节点自由度在转化矩阵中的行索引位置
for i = 1:length(DOF_L)
    index_L_row(i) = find(DOF_A==DOF_L(i));
end
for i = 1:length(DOF_R)
    index_R_row(i) = find(DOF_A==DOF_R(i));
end
for i = 1:length(DOF_B)
    index_B_row(i) = find(DOF_A==DOF_B(i));
end
for i = 1:length(DOF_T)
    index_T_row(i) = find(DOF_A==DOF_T(i));
end
for i = 1:length(DOF_BL)
    index_BL_row(i) = find(DOF_A==DOF_BL(i));
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
for i = 1:length(DOF_INTER)
    index_INTER_row(i) = find(DOF_A==DOF_INTER(i));
end
% 寻找边节点自由度在转化矩阵中的列索引位置
for i = 1:length(DOF_L)
    index_L_col(i) = find(DOF_INDENPENT==DOF_L(i));
end
for i = 1:length(DOF_B)
    index_B_col(i) = find(DOF_INDENPENT==DOF_B(i));
end
for i = 1:length(DOF_BL)
    index_BL_col(i) = find(DOF_INDENPENT==DOF_BL(i));
end
for i = 1:length(DOF_INTER)
    index_INTER_col(i) = find(DOF_INDENPENT==DOF_INTER(i));
end
end