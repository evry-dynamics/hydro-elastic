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

%% �ҳ�Free surface�ĵ㼰���ɶȲ���x&y�����С��������
id_Freesurface = find( co(:,3)==zmax );   %�ҵ�Freesurface�߽��ϵĽڵ���
[~,idx] = sort( co(id_Freesurface,1) );   %��Freesurface�߽��ϵĽڵ㰴x���д�С��������
id_Freesurface = id_Freesurface(idx);     %�������
[~,idx] = sort( co(id_Freesurface,2) );   %��Freesurface�߽��ϵĽڵ㰴y���д�С��������
id_Freesurface = id_Freesurface(idx);     %�������
node_Freesurface = [(1:length(id_Freesurface))',co(id_Freesurface,1),co(id_Freesurface,2)];        %Freesurface�߽�Ľڵ���
node_FG = id_Freesurface; % ����Һ��Ľڵ��Ӧ����������Ľڵ���


%% �ҳ�BottomLeft�ĵ㼰���ɶȲ���z�����С��������
id_BL = find( co(:,1)==xmin & co(:,2)==ymin  );%�ҵ�BottomLeft�߽��ϵĽڵ���
[~,idx] = sort( co(id_BL,3) );%��BottomLeft�߽��ϵĽڵ���д�С��������
id_BL =id_BL(idx);           %�������
node_BL = id_BL;           
DOF_BL = setdiff(node_BL,node_FG);           %BottomLeft�߽������ɶ�

NZ=length(DOF_BL);

id_NX = find( co(:,2)==xmin & co(:,3)==zmin);%�ҵ�BottomLeft�߽��ϵĽڵ���
NX=length(id_NX);
% NX=b/c*10*(NZ-1)+1; % ����߶ȷŴ�10��
% NX=21; % ����߶�100mm����С10��

%% �ҳ�BottomRight�ĵ㼰���ɶȲ���z�����С��������
id_BR = find( co(:,1)==xmax & co(:,2)==ymin  );%�ҵ�BottomrRight�߽��ϵĽڵ���
[~,idx] = sort( co(id_BR,3) );%��BottomrRight�߽��ϵĽڵ���д�С��������
id_BR =id_BR(idx);           %�������
node_BR = id_BR;   
DOF_BR = setdiff(node_BR,node_FG);           %BottomrRight�߽������ɶ�

%% �ҳ�TopRight�ĵ㼰���ɶȲ���z�����С��������
id_TR = find( co(:,1)==xmax & co(:,2)==ymax  );%�ҵ�TopRight�߽��ϵĽڵ���
[~,idx] = sort( co(id_TR,3) );%��TopRight�߽��ϵĽڵ���д�С��������
id_TR =id_TR(idx);           %�������
node_TR = id_TR;   
DOF_TR = setdiff(node_TR,node_FG);           %TopRight�߽������ɶ�

%% �ҳ�TopLeft�ĵ㼰���ɶȲ���z�����С��������
id_TL = find( co(:,1)==xmin & co(:,2)==ymax  );%�ҵ�TopLeft�߽��ϵĽڵ���
[~,idx] = sort( co(id_TL,3) );%��TopLeft�߽��ϵĽڵ���д�С��������
id_TL =id_TL(idx);           %�������
node_TL = id_TL;   
DOF_TL = setdiff(node_TL,node_FG);           %TopLeft�߽������ɶ�

%% �ҳ�Left�ĵ㼰���ɶȲ���z&y�����С��������
id_L = find( co(:,1)==xmin );%�ҵ�Left�߽��ϵĽڵ���
[~,idx] = sort( co(id_L,3) );%��Left�߽��ϵĽڵ㰴z���д�С��������
id_L = id_L(idx);             %�������
[~,idx] = sort( co(id_L,2) );%��Left�߽��ϵĽڵ㰴y���д�С��������
id_L = id_L(idx);             %�������
node_L = id_L;   
node_L = setdiff(node_L,[node_L(1:1:NZ);node_L(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_L,3) );%��Left�߽��ϵĽڵ㰴z���д�С��������
node_L = node_L(idx);             %�������
[~,idx] = sort( co(node_L,2) );%��Left�߽��ϵĽڵ㰴y���д�С��������
node_L = node_L(idx);             %�������
DOF_L = setdiff(node_L,node_FG);                   %Left�߽������ɶ�

%% �ҳ�Right�ĵ㼰���ɶȲ���z&y�����С��������
id_R = find( co(:,1)==xmax );%�ҵ�Right�߽��ϵĽڵ���
[~,idx] = sort( co(id_R,3) );%��Right�߽��ϵĽڵ㰴z���д�С��������
id_R = id_R(idx);             %�������
[~,idx] = sort( co(id_R,2) );%��Right�߽��ϵĽڵ㰴y���д�С��������
id_R = id_R(idx);             %�������
node_R = id_R;   
node_R = setdiff(node_R,[node_R(1:1:NZ);node_R(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_R,3) );%��Right�߽��ϵĽڵ㰴z���д�С��������
node_R = node_R(idx);             %�������
[~,idx] = sort( co(node_R,2) );%��Right�߽��ϵĽڵ㰴y���д�С��������
node_R = node_R(idx);             %�������
DOF_R = setdiff(node_R,node_FG);                   %Right�߽������ɶ�

%% �ҳ�Bottom�ĵ㼰���ɶȲ���z&x�����С��������
id_B = find( co(:,2)==ymin );%�ҵ�Bottom�߽��ϵĽڵ���
[~,idx] = sort( co(id_B,3) );%��Bottom�߽��ϵĽڵ㰴z���д�С��������
id_B = id_B(idx);             %�������
[~,idx] = sort( co(id_B,1) );%��Bottom�߽��ϵĽڵ㰴x���д�С��������
id_B = id_B(idx);             %�������
node_B = id_B;   
node_B = setdiff(node_B,[node_B(1:1:NZ);node_B(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_B,3) );%��Bottom�߽��ϵĽڵ㰴z���д�С��������
node_B = node_B(idx);             %�������
[~,idx] = sort( co(node_B,1) );%��Bottom�߽��ϵĽڵ㰴x���д�С��������
node_B = node_B(idx);             %�������
DOF_B = setdiff(node_B,node_FG);                   %Bottom�߽������ɶ�

%% �ҳ�Top�ĵ㼰���ɶȲ���z&x�����С��������
id_T = find( co(:,2)==ymax );%�ҵ�Top�߽��ϵĽڵ���
[~,idx] = sort( co(id_T,3) );%��Top�߽��ϵĽڵ㰴z���д�С��������
id_T = id_T(idx);             %�������
[~,idx] = sort( co(id_T,1) );%��Top�߽��ϵĽڵ㰴x���д�С��������
id_T = id_T(idx);             %�������
node_T = id_T;   
node_T = setdiff(node_T,[node_T(1:1:NZ);node_T(NZ*(NX-1)+1:1:NZ*NX)]);
[~,idx] = sort( co(node_T,3) );%��Top�߽��ϵĽڵ㰴z���д�С��������
node_T = node_T(idx);             %�������
[~,idx] = sort( co(node_T,1) );%��Top�߽��ϵĽڵ㰴x���д�С��������
node_T = node_T(idx);             %�������
DOF_T = setdiff(node_T,node_FG);                   %Top�߽������ɶ�                                 

%% �ҳ�FSI�ĵ㼰���ɶȲ���x&y�����С��������
% id_FSI = find( co(:,3)==zmax );   %�ҵ�FSI�߽��ϵĽڵ��� %����λ�ڹ����·�
id_FSI = find( co(:,3)==zmin );   %�ҵ�FSI�߽��ϵĽڵ��� %����λ�ڹ����Ϸ�
[~,idx] = sort( co(id_FSI,1) );   %��FSI�߽��ϵĽڵ㰴x���д�С��������
id_FSI = id_FSI(idx);             %�������
[~,idx] = sort( co(id_FSI,2) );   %��FSI�߽��ϵĽڵ㰴y���д�С��������
id_FSI = id_FSI(idx);             %�������
nodeF_FSI = id_FSI;               %FSI�߽�Ľڵ���
              

%% ��Free surface�߽��ϵĽڵ㰴����ʱ��˳���������
nx=size(node_Freesurface,1)^0.5;ny=size(node_Freesurface,1)^0.5; %%ֻ������x��y����Ԫ������ȵ���� Free surface�߽���x��y�����ϵĽڵ�����nx=ny=21
serialG = reshape(node_Freesurface(:,1),nx,ny); % ��Ϊ��y�����С�������еĵ㣬��Ϊ��x�����С�������еĵ㡣
% serialG = serialG'; %ת�ú���Ϊ��x�����С�������еĵ㣬��Ϊ��y�����С�������еĵ㡣

n=size(node_Freesurface,1)^0.5-1;  % Free surface�߽��ϵ�ÿ��ÿ�еĵ�Ԫ��
element_Freesurface=zeros(n^2,5); % ��ʼ����Ԫ����
a=0; % �ع�Freesurface�߽��ϵĵ�Ԫʱ�ļ�����������ʱ��ʼ����
for i=1:n
    for j=1:n
        a=a+1;
        %ÿ��Freesurface�߽絥Ԫ�Ľڵ���Ϊ��ʱ������
        element_Freesurface(a,2) = serialG(i,j);
        element_Freesurface(a,3) = serialG(i+1,j);
        element_Freesurface(a,4) = serialG(i+1,j+1);
        element_Freesurface(a,5) = serialG(i,j+1);
    end
end
element_Freesurface(:,1)=(1:a)';


% �������ɶ�����
DOF_ALL   = setdiff(node_f(:,1),node_FG); %��������Һ�����ɶ�
% �ڲ��ڵ����ɶ�����
DOF_INTER = setdiff(DOF_ALL,[DOF_BL;DOF_TR;DOF_TL;DOF_BR;DOF_B;DOF_L;DOF_T;DOF_R]);

%% �γ�P����
% �����ڵ����ɶ�����
DOF_EDGE       = [DOF_L;DOF_B;DOF_BL];
DOF_EDGE_SORT  = DOF_EDGE ;
%edge_dof_num   = length(DOF_EDGE); % �߽ڵ����ɶȸ���
% �����ڵ����ɶ�
DOF_INDENPENT  = [DOF_L;DOF_B;DOF_BL;DOF_INTER];
edge_dof_num   = length(DOF_INDENPENT);

% �ܱ߽�ڵ����ɶ�
DOF_A  = [DOF_L;DOF_R;DOF_B;DOF_T;DOF_BL;DOF_BR;DOF_TR;DOF_TL;DOF_INTER];
% DOF_A = sort(DOF_A); % ���ɶ�����
% Ѱ���ܱ߽�ڵ����ɶ���ת�������е�������λ��
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
% Ѱ�ұ߽ڵ����ɶ���ת�������е�������λ��
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