function cc = remove_lower_half(cc)
    for y=1:size(cc,1)
        for x=1:y
            cc(y,x) = 0;
        end
    end
end
    