function model = readInput(fileName)

%---------------------------------------------------------------- 
% function readMesh : 
% Read the model text file and save all the data sets in a cell 
% structure called 'model'.
%
% Author : CONG Yu - Universite Paris Saclay  
%          v1 12/2019
%----------------------------------------------------------------

% Pre determining the number of lines in the file  
% clearvars;   fileName = 'cas_pont.inp'

% disp(['Reading the model "',fileName,'"...'])

fid=fopen(fileName);
count = 0;
while true
    if ~ischar( fgetl(fid) ); break; end
    count = count + 1;
end
frewind(fid)

% Read file content to tlines cells
tlines=cell(count,1);
tline = fgetl(fid);
ii = 1;
while ischar(tline)
    tlines{ii,1} = tline;
    tline = fgetl(fid);
    ii = ii + 1;
end
fclose(fid);
% ---------------------------------------


% Remove empty lines
tlines = tlines(~cellfun(@isempty,tlines));


% Remove lines with only whitespaces
comLines = regexp(tlines,'^\s*$','match','once');
tlines = tlines(cellfun(@isempty, comLines));


% Remove comment lines starting with **
comLines = regexp(tlines,'^(\*){2,}','match','once');
tlines = tlines(cellfun(@isempty, comLines));

% ---------------------------------------


% Find block title lines starting with *
titLines = regexp(tlines,'^\*.*$','match','once');

% Mask for lines of title
titLineMask = ~cellfun(@isempty, titLines);

% Extract the start and end of each block
blockseed = [ find(titLineMask); length(tlines)+1 ]';
nblock = length(blockseed)-1;

% Allocate result cells
model = cell(nblock,2);

for i = 1:nblock
    
    % Get data block sizes
    interval = blockseed(i)+1 : blockseed(i+1)-1;
   
    % allocate data block 
    dataCell   = cell(length(interval),1 );
    
    % Convert lines of data one by one
    for il = 1:length(interval)
        
        tline = strtrim(tlines{interval(il)});
        
        [res,done] = str2num(tline);
        
        if done
          dataCell{il,1} = res;
        else
          dataCell{il,1} = tline;
        end
          
    end
      
    % Store data to result
    model{i,1} = char(titLines(blockseed(i)));  % title block
    model{i,2} = dataCell;                      % data block
end


%--
end
