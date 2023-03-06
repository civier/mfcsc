function [sc,labels] = cmtk_to_mrtrix3_order(in_sc,in_labels)

FULL_IDX = get_index;
for i=1:length(FULL_IDX)
      REV_IDX(FULL_IDX(i)) = i;
end 

[sc, labels] = reorder(in_sc,in_labels, REV_IDX);