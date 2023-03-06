function [sc,labels] = reorder(in_sc,in_labels,FULL_IDX)

    FULL_IDX = FULL_IDX(1:length(in_sc));
    
    sc = in_sc(FULL_IDX,FULL_IDX);

    if exist('in_labels')
        labels = in_labels(FULL_IDX);
    end
    
    % nan lower part of functional and structual matrices
    % and transfer value to other half if necessary
    % This is important because the switch from MRtrix ordering to Rosenthal
    % ordering might move some values to the lower triangle
    %
    % NOTICE: I use 0s rather than NaNs for lower trianlge
    for i=1:size(sc,1)
        for j=i:size(sc,2)
            % in the converted natrix, values will be in eithr ij or ji (and the other will be NaN or 0); if
            % in ji, so move then to ij
            if sc(j,i) > 0 || sc(j,i) < 0 % this way it will catch any number or NaN (cannot use ~= because NaN~=0 is true)
                sc(i,j) = sc(j,i);
            end
            sc(j,i) = NaN;
        end
    end    
    
end