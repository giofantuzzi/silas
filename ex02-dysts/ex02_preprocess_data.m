% Preprocess data for the dysts library
% Authors: Albert Alcalde and Giovanni Fantuzzi

% clean up
clear, clc

% Directories
dir_orig  = './python-scripts/data-raw';
dir_test  = './data-test-basic';
dir_train = './data-train-basic';

% Get available data files
data_files = dir([dir_orig, filesep, '*.hdf5']);
nFiles = numel(data_files);

% Loop over systems
for k = 1:nFiles
    
    % Get data from source file
    src_name = data_files(k).name;
    file_orig  = [dir_orig, filesep, src_name];
    dt = h5read(file_orig, '/dt');
    f = h5read(file_orig, '/f');
    t = h5read(file_orig, '/t');
    x = h5read(file_orig, '/x');

    % Number of datapoints and system dimension
    [nVars, nPoints] = size(x);

    % Destination files
    dst_name = sprintf('%02i_%s', nVars, src_name);
    file_train = [dir_train, filesep, dst_name];
    file_test  = [dir_test, filesep, dst_name];

    % save test data
    if ~exist(file_test, 'file')
        h5create(file_test, '/dt', [1 1]);
        h5create(file_test, '/f', [nVars, nPoints]);
        h5create(file_test, '/t', nPoints);
        h5create(file_test, '/x', [nVars, nPoints]);
        h5write(file_test, '/dt', dt);
        h5write(file_test, '/f', f(:,1:nPoints));
        h5write(file_test, '/t', t(1:nPoints));
        h5write(file_test, '/x', x(:,1:nPoints));

    end

    % save half of the test data as training data
    if ~exist(file_train, 'file')
        nPoints = ceil( 0.5 * nPoints );
        h5create(file_train, '/dt', [1 1]);
        h5create(file_train, '/f', [nVars, nPoints]);
        h5create(file_train, '/t', nPoints);
        h5create(file_train, '/x', [nVars, nPoints]);
        h5write(file_train, '/dt', dt);
        h5write(file_train, '/f', f(:,1:nPoints));
        h5write(file_train, '/t', t(1:nPoints));
        h5write(file_train, '/x', x(:,1:nPoints));

    end
end