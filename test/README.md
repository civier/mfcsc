
      testlat -- test laterality
      ==========================

      USAGE

             testlat(MFCFC_INPUT_DIR,OUTPUT_DIR)

      PREREQUISITES

             Matlab's Statistics and Machine Learning Toolbox

      DESCRIPTION

             testlat compares the mismatch between FC and SC in the left hemisphere
             to the mismatch in the right hemisphere.

             testlat recevies a group of connectivity matrices (Desikan-Killiany ROIs) with MFCSC values
             and for each connection in the matrix, compare between left and right hemispheres.
             testlat returns the results of the statsitical tests, as well as which
             tests are significant after benferroni correction for multiple comparisons.

             testlat intends to analyse a homogenous group of participants. However,
             it also outputs the difference between hemispheres within each participant, to 
             assist with other statistical designs.

             testlat can also be used for connectivity matrices with other metrics other than mfcsc.
             
      INPUT

             The input files are connectivity matrices with a MFCSC value for each connection, 
             generated using the mfcsc tool.
             They should have 84 rows and 84 columns, in par with the Desikan-Killiany atlas.

             Rows and columns of the connetivity matrices need to be ordered according to the 
             order specified here: https://osf.io/q7v9t
             This is the order used by the MRtrix3 software: https://www.mrtrix.org/
             IMPORTANT: When using mfcsc to generate the MFCSC connectivity matrices, it is thus
             recommended to organise the rows and columns of the input functional and structural
             connectivity matrices already in this order.

      OUTPUT

             The output files will have half the rows and columns of the input matrices (42).
             They will be ordered identically to columns/rows 1-42 of the input matrices.
             Each cell represent a pair of bilateral (left and right) connections.
             Note that all output matrices only include values for cells which are within 
             the mask (mask_L_and_R.csv). Cells outside of the mask get the value of -99.

             The main output of testlat consists of the files:

                 sig_all.csv -    binary connectivity matrix. The cells where MFCSC
                                  is significantly different between the left and right connections
                                  of the pair are 1 (after correction for multiple comparisons over 
                                  all cells in the upper triangle of the matrix).

                 sig_all_labels.txt - ROI labels for all signficant pairs of connections.
                                  Labels used are those of the left connection in each pair.

                 pvalue_all.csv   - connectivity matrix of p-values. The p value is each cell
                                    is the result of the t-test comparing the MFCSC of the left connection
                                    to the MFCSC of the right connection in all input connectomes.
                                    p values are not corrected for multiple comparisons.

                 mask_L_and_R.csv - binary connectivity matrix. The cells where both left and right
                                    connections in the pair are in the whole-brain mask are 1.
                                    
                 LminusR-mean.csv - a connectivity matrix. Every cell is the average values of
                                    the right connection's MFCSC subtracted from the left connection. 

             Additional files are:

                 sig_*.csv - significant tests matrices for 6 different classes of relationships between
                         left and right average values. Results are corrected.
                 
                 pvalue_*.csv - p-values matrices for the 6 different classes of relationships.

                 sig_*_labels.txt - ROI labels for significant pairs of connections for the 6 different classes
                                of relationships.

                 individual_connectomes/LminusR* - a connectivity matrix for each input connectome.
                                                   Each cell is the difference between the MFCSC of left and right
                                                   connections in the pair for the individual participant.

      ARGUMENTS

             MFCSC_INPUT_DIR

             The directory containing the connectivity matrices with MFCSC values.

             OUTPUT_DIR

             The directory where the output files will be written to.

      DEVELOPER

             Oren Civier (orenciv@gmail.com)
             https://www.swinburne.edu.au/research/our-research/access-our-research/find-a-researcher-or-supervisor/researcher-profile/?id=ocivier

      CITATIONS

             When using testlat, authors should cite:
             Civier O, Sourty M, Calamante F (2023) MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, 
             and its application for detecting hemispheric functional specialisations. Scientific Reports https://doi.org/10.1038/s41598-022-17213-z

      ACKNOWLEDGMENTS

             National Health and Medical Research Council of Australia (grant numbers APP1091593 andAPP1117724)
             Australian Research Council (grant number DP170101815)
             National Imaging Facility (NIF), a National Collaborative Research Infrastructure Strategy (NCRIS) capability at Swinburne Neuroimaging, 
             Swinburne University of Technology.

