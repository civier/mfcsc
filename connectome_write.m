function connectome_write(filename,mat,property,value)
    global is_symmetrical_global;
    if is_symmetrical_global
        mat = copy_upper_to_lower(mat);
    else
        mat = zero_lower_half(mat);
    end
    % zero diagonal
    mat(logical(eye(size(mat)))) = 0;
    
    if exist('property')
        dlmwrite(filename,mat,property,value);
    else
        dlmwrite(filename,mat);
    end        
end
