function [sc,labels] = mrtrix3_to_cmtk_order(in_sc,in_labels);

FULL_IDX = get_index;
    
[sc, labels] = reorder(in_sc,in_labels,FULL_IDX);
    
    