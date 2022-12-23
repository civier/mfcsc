function x = mirror_mat(x)
    N = length(x);
    for i=1:N
        x(i,i) = 0;
        for j=1:N
            if isnan(x(i,j))
                x(i,j) = x(j,i);
            end
        end
    end
end
