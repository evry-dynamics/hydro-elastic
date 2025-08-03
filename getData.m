function data = getData(model,dataName)

%---------------------------------------------------------------- 
% function getData : 
% Extract data with the name 'dataName' from the model, for ex. 
% node coordinates 'coor'. If several blocks are defined for the 
% same data, for example when several materials are defined, 
% then the exported data will contain all the blocks in a cell 
% structure.
%----------------------------------------------------------------
% extract data containing "dataName" in the title


% (?i) case insensative
matchLines = regexp(model(:,1),['(?i)^.*',dataName,'.*$'],'match','once');

index = find(~cellfun(@isempty, matchLines));
lindex = length(index);

data = {};
if lindex == 0
    disp(['No "',dataName,'" matches in the model.'])
    return;
end

data = model(index, 2);
% disp(['Extracted ',num2str(length(index)),' "',...
%               dataName,'" object(s) from the model.'])
end














% life in a boat