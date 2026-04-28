function save_results(filename, varargin)
% Check if the file already exists and delete it to avoid appending to an old file
if isfile(filename)
    delete(filename);
end

% Write the structures `S1` and `S2` to the HDF5 file under separate groups
for i = 1:numel(varargin)
    writeStructToHDF5(varargin{i}, inputname(i+1), filename);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to write a structure to an HDF5 file under a specified group
function writeStructToHDF5(structData, groupName, filename)

% Loop over structure's field names
fields = fieldnames(structData);
for i = 1:numel(fields)
    fieldName = fields{i};
    data = structData.(fieldName);
    fullPath = ['/' groupName '/' fieldName];

    % Operate according to data type
    if iscell(data)
        for j = 1:numel(data)
            cellData = data{j};
            cellPath = sprintf('%s/f%d', fullPath, j);
            if ischar(cellData)
                h5create(filename, cellPath, size(cellData), 'DataType', 'string');
                h5write(filename, cellPath, cellData);
            else
                h5create(filename, cellPath, size(cellData), 'DataType', class(cellData));
                h5write(filename, cellPath, cellData);
            end
        end
 
    elseif ischar(data)
        data = string(data);
        h5create(filename, fullPath, size(data), 'DataType', 'string');
        h5write(filename, fullPath, data);

    elseif isstruct(data)
        % call recursively
        writeStructToHDF5(data, fullPath, filename)

    elseif ~isempty(data)
        h5create(filename, fullPath, size(data), 'DataType', class(data));
        h5write(filename, fullPath, data);
    end
end
end
