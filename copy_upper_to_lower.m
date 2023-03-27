function x = copy_upper_to_lower(x)
    N = length(x);
    for row=2:N
        for col=1:(row-1)
            x(row,col) = x(col,row);
        end
    end
end
