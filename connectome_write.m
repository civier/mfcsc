function connectome_write(filename,mat)
    global is_symmetrical_global;
    if is_symmetrical_global
        mat = copy_upper_to_lower(mat);
    else
        mat = zero_lower_half(mat);
    end
    % zero diagonal
    mat(logical(eye(size(mat)))) = 0;
    
    dlmwrite(filename,mat);
end