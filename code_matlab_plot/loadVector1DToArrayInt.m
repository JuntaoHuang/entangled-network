
function vec = loadVector1DToArrayInt(filename)
    % Open the binary file
    fid = fopen(filename, 'rb');
    
    if fid == -1
        error('Error opening file: %s', filename);
    end
    
    % Read the size of the vector (stored as a size_t in C++)
    size = fread(fid, 1, 'uint64');
    
    % Read the vector data
    vec = fread(fid, size, 'int32');
    
    % Close the file
    fclose(fid);
end
