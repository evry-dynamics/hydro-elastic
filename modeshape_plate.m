function [f]=modeshape_plate(u,node,element,order,method)

u = normalize(u,'range');

node_x = zeros(4,size(element,1));
node_y = zeros(4,size(element,1));
U = zeros(4,size(element,1));

for i=1: size(element,1)
    nu = element(i,2:5); % 单元的节点编号 dof = [n1, n2, n3, n4]
    dof_z = 5*nu-2; % 节点的z方向自由度编号 nz = [dof_nz1, dof_nz2, dof_nz3, dof_nz4]
    node_x(:,i) = node(nu,2); % 节点x坐标 node_x = [n1_x, n2_x, n3_x, n4_x]
    node_y(:,i) = node(nu,3); % 节点y坐标 node_y = [n1_y, n2_y, n3_y, n4_y]
    U(:,i) = u(dof_z); % 节点的z方向位移 U = [n1_u, n2_u, n3_u, n4_u] 
end

set(0,'defaultfigurecolor','w') %设置背景色为白色
%% 绘制固体振型
f=figure;
patch(node_x,node_y,U); %填充
shading interp; %色彩平滑
colormap(jet);
c=colorbar;
c.Ruler.TickLabelFormat = '%.2f';
% clim([0 1])  % 颜色范围
title(method)
% title([method,' Solid uz', ' order=', num2str(order)],'Interpreter', 'none')
% ylim([0,h+0.1])
% xlim([-0.2,0.2])
yticks([])
xticks([])
box
set(gca,'DataAspectRatio',[1 1 1])