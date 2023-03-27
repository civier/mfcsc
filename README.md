 

     mfcsc -- Mismatch between Functional Connectivity and Structural Connectivity  
     =============================================================================

     JOURNAL ARTICLE
     
          Civier O, Sourty M, Calamante F (2023) MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, 
          and its application for detecting hemispheric functional specialisations. Scientific Reports https://doi.org/10.1038/s41598-022-17213-z

     VERSION

          1.1

     USAGE

         From Neurodesk (no Matlab license required!):
             See https://osf.io/d7j9n/

         From Matlab IDE:
             mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dir)         

         From the command-line:
             Change to the folder with the code, and then enter the following command -
             matlab -batch "mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dor)"

         The 4 arguments in caps are mandatory. 

         IMPORTANT: before running mfcsc, it is recommended to open FC_SC_LIST 
         and visually check that the two columns have matching
         participant connectomes.

     PREREQUISITES:

         Matlab 9.8.0.1721703 (R2020a) Update 7 or newer (older versions of Matlab should work as well, but not tested)
         Matlab's Curve Fitting toolbox

             OR

         Neurodesktop 20230324 or newer 
         (older versions of Neurodesktop should work as well, but not tested; to use in older versions, load from the terminal using: 'ml mfcsc/1.1')

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

     INPUT

         The input connectivity matrices should be comma seperated values that describe symmetrical 
         connectomes (i.e., connections do not have directionality). 
         Only the upper right triangle of the connectivity matrices is
         consulated. This also excludes the main diagonal.

         The connectivity matrices must have an even number of regions N: regions numbers 
         1 to N/2 for one hemisphere, and regions number N/2+1 to N for the other.
         Note that the order of regions for each hemisphere may be different, i.e., 
         regions number i and N/2+i do not have to be homologous.

     OUTPUT

         IMPORTANT: before consulting the output, ensure that the file1 and file2 parts 
         of the filenames below represent matching connectivity
         matrices!
         
         The output of mfcsc consists of one connectome file for each participant:

                 mFCSC-file1-file2-masked.csv - connectome of mFCSC values for the
                                                participant whose FC and SC
                                                connectomes are stored in file1 and file2
                                                respectively (excluding
                                                the files extensions)

         The main output of mfcsc consists of the file:

                 mask-final.csv - final mask indicating the connections to which mFCSC is calculated



         There are also several misc files in the group_connectomes subdir:

             transformed_SC_avg - the average transfered SC connectome
             FC_avg - the average FC connectome
             mask-direct_SC_is_shortest_path.csv - mask of connections in which the path length of the direct connection (1/transformed_SC) is shorter than any other indirect path between the two regions
             SC_avg - the average SC connectome

     MANDATORY ARGUMENTS

         FC_SC_LIST (path)

         Path to a tab-separated file with two columns.
         The first column lists the files with the connectivity matrices of the FC connectomes
         The second column lists the connectivity matrix files for the matching SC connectomes

         An easy method to generate the FC_SC_LIST on Linux or MacOS is to:
         1) include a matching participant ID in the filenames of both FC and SC connectomes 
            (plus optional fixed prefixes and suffixes for each modality).
         2) put all FC connectomes in one folder and all SC connectomes in another,
            ensuring that there are no other files in these folders.
         3) change the working directory to the folder with FC connectomes, and run:
             ls | sort -n > /tmp/fc_list
         4) Change the working directory to the folder with the SC connectomes, and run:
             ls | sort -n > /tmp/sc_list
         5) Run:
             paste /tmp/fc_list /tmp/sc_list > path_to_filename
         6) Provide path_to_filename as the FC_SC_LIST argument

         FC_INPUT_DIR (path)

         The directory containing the files with the FC connectivity matrices

         SC_INPUT_DIR (path)

         The directory containing the files with the SC connectivity matrices

         OUPTUT_DIR (path)

         The output directory where the mask and mFCSC files are to be
         saved

     OPTIONAL ARGUMENTS

         not_in_mask_value (any number)

         Value that will be assigned to cells in the output matrices that are
         not in the mask. By deafult it is set to -99 to make sure
         people do not report the values in these cells.
         Can be set to 0 to prevent this value from affecting color scaling of plots.    

         is_contra (true or false)
         
         By default, mFCSC is calculated for ipsilateral connections
         (is_contra = false or omitted)
         set is_contra to true to calculate mFCSC for contralateral connections instead.

         is_keep_neg_fc (true or false)

         By default, mFCSC removes cells that have negative mean FC before 
         fitting the model used to transform SC
         (is_remove_negative_fc = false or omitted)
         set is_keep_negative_fc to true to keep them

         is_symmetrical (true or false)

         By default, only the upper right triangle of the output matrices is
         populated, with the lower triangle being zeroed out.
         (is_symmetrical = false or omitted)
         set is_symmetrical to true to save symmetrical matrices instead by
         mirroring the upper right triangle into the bottom left triangle
         In both cases, the main diagonal is zeroed out.

         is_figures (true or false)

         By default, do not show figures.
         (is_figures = false or omitted)
         set is_figures to true to print verbose figures with information for QC and debg.
         Not tested. Use at your own risk!

         bct_dir (path)

         If specified, Matlab will look for the Brain Connectome Toolbox (BCT)
         in this directory instead of the BCT version supplied with
         MFCSC (2017/01/15). BCT is required in order to calculate the mask.
         This argument is not available in the Neurodesk version, where only the supplied
         BCT can be used.

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

      EXAMPLE

         For the processing performed in the journal article, see:
         https://osf.io/d7j9n/ under "TESTING MFCSC INSTALLATION"

      DEVELOPER

         Oren Civier (orenciv@gmail.com)
         https://www.swinburne.edu.au/research/our-research/access-our-research/find-a-researcher-or-supervisor/researcher-profile/?id=ocivier

      CITATIONS

         When using mfcsc, authors should cite:

             Civier O, Sourty M, Calamante F (2023) MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, 
             and its application for detecting hemispheric functional specialisations. Scientific Reports https://doi.org/10.1038/s41598-022-17213-z

             Rubinov M, Sporns O (2010) Complex network measures of brain connectivity: Uses and interpretations. NeuroImage 52:1059-69.%

      ACKNOWLEDGMENTS

         National Health and Medical Research Council of Australia (grant numbers APP1091593 andAPP1117724)
         Australian Research Council (grant number DP170101815)
         National Imaging Facility (NIF), a National Collaborative Research Infrastructure Strategy (NCRIS) capability at Swinburne Neuroimaging, 
         Swinburne University of Technology.

