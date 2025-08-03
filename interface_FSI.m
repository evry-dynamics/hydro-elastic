function [node_s_FSI,node_f_FSI,nv] = interface_FSI(node_s,node_f)
% interface nodes on the solid side
id_FSI_s = node_s;
[~,index] = sortrows(id_FSI_s(:,2));
node_s_FSI = node_s(index,:);
% interface nodes on the fluid side
id_FSI_f = find( node_f(:,4)==max(node_f(:,4)) ); 
[~,index] = sortrows(node_f(id_FSI_f,2));
node_f_FSI = node_f(index,:);