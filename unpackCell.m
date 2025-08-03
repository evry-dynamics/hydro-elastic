function  [CELL]= unpackCell(cell_0)

ncell_0 = size(cell_0,1);

% count number of subcells
ncell_1 = 0;
for ii = 1:ncell_0
    ncell_1 = ncell_1 + size(cell_0(ii,1),1);
end
CELL = cell(ncell_1,1);

n = 1;
for ii = 1:ncell_0
    cell_1 = cell_0{ii,1};
    
    if iscell(cell_1)
    
        for jj = 1 : size(cell_1,1)
            CELL{n,1} =  cell_1{jj,1};
            n = n+1;
        end
        
    else
        CELL{ii,1} = cell_1;
    end
        
end



% pack up all my load, so long, good bye 

end
