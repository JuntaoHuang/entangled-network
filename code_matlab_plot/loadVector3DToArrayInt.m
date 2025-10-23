
function data = loadVector3DToArrayInt(filename)

    % Open the file for reading
    fid = fopen(filename, 'rb');
    
    % Read dimensions
    x_dim = fread(fid, 1, 'int32');
    y_dim = fread(fid, 1, 'int32');
    z_dim = fread(fid, 1, 'int32');
    
    % Read the data
    data = fread(fid, x_dim * y_dim * z_dim, 'int32');
    
    % Reshape into 3D array
    data = reshape(data, [z_dim, y_dim, x_dim]);
    data = permute(data, [3, 2, 1]); % Rearrange dimensions to match original

    % Close the file
    fclose(fid);
        
end