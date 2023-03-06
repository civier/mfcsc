Oren Civier
 
     XXX

       calulate_laterality
       ===================
 
       Note: difference between hemispheres will be only caclulated to
       connections which are within the mask.csv in both hemispheres.
 
       For HCP data in the respository, there are 512 connections that are
       1 in the left hemisphere, and 497 connections that are 1 in the
       right hemisphere. However, their overlap only gives 475 connections
       which are 1 in both hemispheres.
 
       B & E (flips) has the most significant p-values that are not in the
       mask. This is because such a big difference usually does not exist,
       and if it does, it is usually due to FC in the strong side being affected by
       strong indirect structural connections.
       D (L < R in both sides) has tons of significant p-values that are
       not in the mask. I'm not sure why? Maybe there are low FC in both
       because there is no effective connectivity in both, and the relatively stronger
       SC is not real; it is actually weaker than the indirect structural connectivity.
 
       The only difference from the article is in sig_neg_L_st_R (Fig. 3,
       panel D). If running the code, it is:
            ctx-lh-parsopercularis - Left-Putamen
            ctx-lh-parstriangularis - Left-Putamen
            ctx-lh-medialorbitofrontal - Left-Amygdala
       But in the paper, it is only the first two:
            ctx-lh-parsopercularis - Left-Putamen
            ctx-lh-parstriangularis - Left-Putamen           
       The difference is due a slight difference in the implementation
       between the code provided here and that used for the paper. Here the
       linear regression inside each individual includes all connections
       where direct SC is the shortest structural path. In contrast, in the
       paper, the linear regression only includes the connections where
       direct SC is the shortest structural path *in both hemispheres*.
%
