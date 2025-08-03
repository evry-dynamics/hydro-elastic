function [node,element] = ReadFile_s(NodeFile,ElementFile)
% ========================================== %
% This matlab file is the part of MindlinPlatePeriodic Folder,
% the function reads the node and element information of Fem
% :param  NodeFile: Node information in finite element analysis 
% :param ElementFile: Element information in finite element analysis
% :return node: node matrix
% :return element: element matrix
% ========================================== %
    node_data = textread(NodeFile);
    element_data = textread(ElementFile);
    node = node_data(:,1:4);
    element = element_data(:,[1,7:10,2]);
end