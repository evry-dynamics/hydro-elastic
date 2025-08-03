function R = FSIMatrix(node_s,nodeS_FSI,node_f,nodeF_FSI)

%% 构建相同的流体界面单元 elementF_fsi 和 固体界面单元 elementS_fsi 矩阵
nx=size(node_s,1)^0.5;ny=size(node_s,1)^0.5;%只适用于x和y方向单元个数相等的情况 %流固耦合界面上x和y方向上的节点数。
serialS = reshape(nodeS_FSI,nx,ny); % 行为按y方向从小到大排列的点，列为按x方向从小到大排列的点。
% serialS = serialS'; %转置后，行为按x方向从小到大排列的点，列为按y方向从小到大排列的点。
serialF = reshape(nodeF_FSI,nx,ny);
% serialF = serialF';

n=size(node_s,1)^0.5-1;  %流固耦合界面上的每行每列的单元数
elementS_fsi=zeros(n^2,4);%初始化固体面上的fsi单元矩阵
elementF_fsi=zeros(n^2,4);%初始化流体面上的fsi单元矩阵
a=0;%重构流固耦合界面上的单元时的计数参数，此时初始化。
for i=1:n
    for j=1:n
        a=a+1;
        %每个流固耦合单元的节点编号为逆时针排列
        elementS_fsi(a,1) = serialS(i,j);
        elementS_fsi(a,2) = serialS(i+1,j);
        elementS_fsi(a,3) = serialS(i+1,j+1);
        elementS_fsi(a,4) = serialS(i,j+1);

        elementF_fsi(a,1) = serialF(i,j);
        elementF_fsi(a,2) = serialF(i+1,j);
        elementF_fsi(a,3) = serialF(i+1,j+1);
        elementF_fsi(a,4) = serialF(i,j+1);
    end
end


%% 面法线向量矩阵 normal
normal=zeros(length(elementS_fsi),5);
for i=1:length(elementS_fsi)
    normal(i,3) = -1;
end

%% gauss点计算流固耦合矩阵 R
R  = zeros(5*length(node_s),length(node_f));%初始化矩阵R
te = size(elementS_fsi,1);%总的流固耦合单元的个数
for i = 1:te
    %分别获取第i个流固耦合单元的流固节点编号
    eS = elementS_fsi(i,:);
    eF = elementF_fsi(i,:);
    
    nS1 = eS(1);    nS2 = eS(2);    nS3 = eS(3);    nS4 = eS(4);
    nF1 = eF(1);    nF2 = eF(2);    nF3 = eF(3);    nF4 = eF(4);
    %提取一个单元里，四个节点的坐标
    x1 = node_s(nS1,2);x2 = node_s(nS2,2);x3 = node_s(nS3,2);x4 = node_s(nS4,2);
    y1 = node_s(nS1,3);y2 = node_s(nS2,3);y3 = node_s(nS3,3);y4 = node_s(nS4,3);

    %根据固体节点编号，计算其自由度编号 edofS
    edofS = [nS1*5-4,nS1*5-3,nS1*5-2,nS1*5-1,nS1*5,...
            nS2*5-4,nS2*5-3,nS2*5-2,nS2*5-1,nS2*5,...
            nS3*5-4,nS3*5-3,nS3*5-2,nS3*5-1,nS3*5,...
            nS4*5-4,nS4*5-3,nS4*5-2,nS4*5-1,nS4*5];
    edofS = edofS';
%     edofS=zeros(5,4);
%     for j = 1 : 5
%         edofS(j,1) = nS1 * 5 - 5 + j;
%         edofS(j,2) = nS2 * 5 - 5 + j;
%         edofS(j,3) = nS3 * 5 - 5 + j;
%         edofS(j,4) = nS4 * 5 - 5 + j;
%     end
%     edofS = edofS(:);%这一步将矩阵转化为一列
    %根据流体节点编号，计算其自由度编号 edofF
    edofF(1)=nF1;edofF(2)=nF2;edofF(3)=nF3;edofF(4)=nF4;
    edofF = edofF(:);%这一步将矩阵转化为一列

    normalE = normal(i,:);   %normal of this interface element 

    % 计算单元 r 矩阵
    p = 1/sqrt(3);
    GassPoint = [-p,-p;
                 -p,  p;
                  p, -p;
                  p,  p];%四边形等参元，此处采用4个高斯积分点参数，权重1
    r=zeros(length(edofS),length(edofF));%初始化
    for k = 1 : size(GassPoint,1)
        s = GassPoint(k,1);%第k个高斯积分点的横坐标
        t = GassPoint(k,2);%第k个高斯积分点的纵坐标
        J = Jacobian(s,t,x1,x2,x3,x4,y1,y2,y3,y4);
        J_det = det (J);
        [N1,N2,N3,N4] = ShapeFunction(s,t);
        Nf=[N1,N2,N3,N4];%流体形函数矩阵
        Ns=[N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0,0,0;
            0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0,0;
            0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0,0;
            0,0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4,0;
            0,0,0,0,N1,0,0,0,0,N2,0,0,0,0,N3,0,0,0,0,N4];%固体形函数矩阵
        NN = Ns'*normalE'*Nf; % (20*5) * (5*1) * (1*4)
        
        %gaussian integration
        r=r+NN*J_det;
    end
    R(edofS,edofF)=R(edofS,edofF)+r;%根据单元自由度，将单元r并入到整体R
end
end