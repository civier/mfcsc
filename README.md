 
     DESCRIPTION

         mfcsc receive pairs of functional and structural connectivity
         matrices, one pair per participant. For each participant, mfsc
         calculates a connectivity matrix that gives the mFCSC metric
         value for every connection. 

         The mFCSC metric can be calculated either for ipsilateral or 
         contralateral connectionsis (but not both), and is ill-defined 
         for some connections. The 'mask.csv' binary matrix mask
         indicates which mFCSC metrics should be included in further analysis. 
         The cells that should not be consulted further are set to -99.

         The input connectivity matrices should be comma seperated values that describe symmetrical 
         connectomes (i.e., connections do not have directionality). 
         Only the upper right triangle of the connectivity matrices is
         consulated. This also excludes the main diagonal.

         The connectivity matrices must have an even number of regions N: regions numbers 
         1 to N/2 for one hemisphere, and regions number N/2+1 to N for the other.
         Note that the order of regions for each hemisphere may be different, i.e., 
         regions number i and N/2+i do not have to be homologous.

     USAGE

         1)  Make sure that your Matlab instllation includes the Curve Fitting toolbox

         2)  Organise input files and generate a FC-SC list (see ARGUMENTS below)

         3)  Verify that the each pair in the FC-SC list is made up of two matching connectomes

         4)  From Matlab IDE:
             Change to the folder with the code, and then enter the following command -
             mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,BCT_DIR)         

             From the command-line:
             Change to the folder with the code, and then enter the following command -
             matlab -batch "mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,BCT_DIR)"

     ARGUMENTS (mandatory)

         FC_SC_LIST

         Full path to a tab-separated file with two columns.
         The first column lists the files with the connectivity matrices of the FC connectomes
         The second column lists the connectivity matrix files for the matching SC connectomes

         FC_INPUT_DIR

         The directory containing the files with the FC connectivity matrices

         SC_INPUT_DIR

         The directory containing the files with the SC connectivity matrices

         OUPTUT_DIR

         The output directory where the mask and mFCSC files are to be
         saved

     FLAGS (optional)

         is_contra
         
         By default, mFCSC is calculated for ipsilateral connections
         (is_contra = false or omitted)
         set is_contra to true to calculate mFCSC for contralateral connections instead.

         is_keep_neg_fc

         By default, mFCSC removes cells that have negative mean FC before 
         fitting the model used to transform SC
         (is_remove_negative_fc = false or omitted)
         set is_keep_negative_fc to true to keep them

         is_symmetrical

         By default, only the upper right triangle of the output matrices is
         populated, with the lower triangle being zeroed out.
         (is_symmetrical = false or omitted)
         set is_symmetrical to true to save symmetrical matrices instead by
         mirroring the upper right triangle into the bottom left triangle
         In both cases, the main diagonal is zeroed out.

         is_figures

         By default, do not show figures.
         (is_figures = false or omitted)
         set is_figures to true to print verbose figures with information for QC and debg.
         Not tested. Use at your own risk!

         BCT_DIR

         If specified, Matlab will look for the Brain Connectome Toolbox (BCT)
         in this directory instead of the BCT version supplied with
         MFCSC (2017/01/15). BCT is required in order to calculate the mask.

     OUTPUT

         The main output of mfcsc consists of the file:
                 mask-final.csv -                final mask indicating the connections to which mFCSC is calculated
         
         As well as one connectome file for each participant:
                 mFCSC-file1-file2-masked.csv - connectome of mFCSC values for the
                                                participant whose FC and SC
                                                connectomes are stored in file1 and file2
                                                respectively (excluding
                                                the files extensions)

         There are also several misc files in the misc subdir:
             transformed_SC_avg -                  the average transfered SC connectome
             FC_avg -                              the average FC connectome
             mask-direct_SC_is_shortest_path.csv - mask of connections in which the path length of the direct connection (1/transformed_SC) is shorter
                                                   than any other indirect path between the two regions
             SC_avg -                              the average SC connectome

     NOTE
         
         Cells of individual SC connectivity matrices that have the value of 0 are
         not transfromed well using the model.
         In most cases these cells will not be included in the mask because the
         direct SC in the average SC connectome is not the shortest path;
         however, in case they are within the mask after all, a warning will
         be issued and they will be assigned the mFCSC value of -999.
         In this case, one approach is to manually exclude these cells from the
         mask by editing 'mask-final.csv' and using the ammended mask in
         further analyses.

	ACKNOWLEDGMENTS

	    When using MFCSC, the following articles should be cited:

		MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, and its application for 
		detecting hemispheric functional specialisations
		Civier O, Sourty M, Calamante F (in press) Scientific Reports https://doi.org/10.1038/s41598-022-17213-z

		Complex network measures of brain connectivity: Uses and interpretations.
		Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.



