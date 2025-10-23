
function data = loadVector2DToArrayInt(filename)
    % Open the file for reading
    fid = fopen(filename, 'rb');
    
    % Read dimensions
    rows = fread(fid, 1, 'int32');
    cols = fread(fid, 1, 'int32');
    
    % Read the data
    data = fread(fid, rows * cols, 'int32');
    
    % Reshape into 2D array
    data = reshape(data, [cols, rows]);
    data = data'; % Transpose to match original row-major order

    % Close the file
    fclose(fid);
end
