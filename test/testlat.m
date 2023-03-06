function testlat(input_folder,output_folder)
%
%      testlat -- test laterality
%      ==========================
%
%      USAGE
%
%             testlat(MFCFC_INPUT_DIR,OUTPUT_DIR)
%
%      DESCRIPTION
%
%             testlat compares the mismatch between FC and SC in the left hemisphere
%             to the mismatch in the right hemisphere.
%
%             testlat recevies a group of connectivity matrices (Desikan-Killiany ROIs) with MFCSC values
%             and for each connection in the matrix, compare between left and right hemispheres.
%             testlat returns the results of the statsitical tests, as well as which
%             tests are significant after benferroni correction for multiple comparisons.
%
%             testlat intends to analyse a homogenous group of participants. However,
%             it also outputs the difference between hemispheres within each participant, to 
%             assist with other statistical designs.
%
%             testlat can also be used for connectivity matrices with other metrics other than mfcsc.
%             
%      INPUT
%
%             The input files will have 84 rows and 84 columns, in par with the Desikan-Killiany atlas.
%             Rows and columns of the connetivity matrices need to be order according to the 
%             order specified here: https://osf.io/q7v9t
%             This is the order used by the MRtrix3 software: https://www.mrtrix.org/
%
%      OUTPUT
%
%             The output files will have half the rows and columns of the input matrices (42).
%             They will be ordered identically to columns/rows 1-42 of the input matrices.
%             Each cell represent a pair of bilateral (left and right) connections.
%             Note that all output matrices only include values for cells which are within 
%             the mask (mask_L_and_R.csv). Cells outside of the mask get the value of -99.
%
%             The main output of testlat consists of the files:
%
%                 sig_all.csv -    binary connectivity matrix. The cells where MFCSC
%                                  is significantly different between the left and right connections
%                                  of the pair are 1 (after correction for multiple comparisons over 
%                                  all cells in the upper triangle of the matrix).
%
%                 sig_all_labels.txt - ROI labels for all signficant pairs of connections.
%                                  Labels used are those of the left connection in each pair.
%
%                 pvalue_all.csv   - connectivity matrix of p-values. The p value is each cell
%                                    is the result of the t-test comparing the MFCSC of the left connection
%                                    to the MFCSC of the right connection in all input connectomes.
%                                    p values are not corrected for multiple comparisons.
%
%                 mask_L_and_R.csv - binary connectivity matrix. The cells where both left and right
%                                    connections in the pair are in the whole-brain mask are 1.
%                                    
%                 LminusR-mean.csv - a connectivity matrix. Every cell is the average values of
%                                    the right connection's MFCSC subtracted from the left connection. 
%
%             Additional files are:
%
%                 sig_* - significant tests matrices for 6 different classes of relationships between
%                         left and right average values. Results are corrected.
%                 
%                 pvalue_* - p-values matrices for the 6 different classes of relationships.
%
%                 sig_*_labels - ROI labels for significant pairs of connections for the 6 different classes
%                                of relationships.
%
%                 individual_connectomes/LminusR* - a connectivity matrix for each input connectome.
%                                                   Each cell is the difference between the MFCSC of left and right
%                                                   connections in the pair for the individual participant.
%
%      ARGUMENTS
%
%             MFCSC_INPUT_DIR
%
%             The directory containing the connectivity matrices with MFCSC values.
%
%             OUTPUT_DIR
%
%             The directory where the output files will be written to.
%
%
%      For HCP data in the respository, there are 512 connections that are
%      1 in the left hemisphere, and 497 connections that are 1 in the
%      right hemisphere. However, their overlap only gives 475 connections
%      which are 1 in both hemispheres.
%
%      B & E (flips) has the most significant p-values that are not in the
%      mask. This is because such a big difference usually does not exist,
%      and if it does, it is usually due to FC in the strong side being affected by
%      strong indirect structural connections.
%      D (L < R in both sides) has tons of significant p-values that are
%      not in the mask. I'm not sure why? Maybe there are low FC in both
%      because there is no effective connectivity in both, and the relatively stronger
%      SC is not real; it is actually weaker than the indirect structural connectivity.
%
%      The only difference from the article is in sig_neg_L_st_R (Fig. 3,
%      panel D). If running the code, it is:
%           ctx-lh-parsopercularis - Left-Putamen
%           ctx-lh-parstriangularis - Left-Putamen
%           ctx-lh-medialorbitofrontal - Left-Amygdala
%      But in the paper, it is only the first two:
%           ctx-lh-parsopercularis - Left-Putamen
%           ctx-lh-parstriangularis - Left-Putamen           
%      The difference is due a slight difference in the implementation
%      between the code provided here and that used for the paper. Here the
%      linear regression inside each individual includes all connections
%      where direct SC is the shortest structural path. In contrast, in the
%      paper, the linear regression only includes the connections where
%      direct SC is the shortest structural path *in both hemispheres*.
%


%dbstop if error;

%     DO_EXTRA = false;

% try to create output directory - so if fails, not continue anymore
if ~exist(output_folder, 'dir')
       mkdir(output_folder);
end

if ~exist([output_folder filesep 'individual_connectomes'], 'dir')
       mkdir([output_folder filesep 'individual_connectomes']);
end

%%%% Load labels from the same folder the testlat.m is in
[folder,~,~] = fileparts(mfilename('fullpath'));
fid = fopen([folder filesep 'MRtrix_labels_Desikan-Killiany.txt']);
in_labels = textscan(fid,'%s','delimiter','\n');
in_labels = in_labels{1};
fclose(fid);

% 1. need to load mismatches from files and convert to Rosenthal ordering
folder_path = input_folder; % Replace with the path to your folder
file_list = dir(fullfile(folder_path, '*.csv')); % List all .mat files in the folder


% first read mask
file_name = 'mask.csv';
file_path = fullfile(folder_path, file_name); % Construct the full file path

disp(['Loaded file ', file_name]);

mask= dlmread(file_path); % Load the data from the file
N = size(mask,1); % get number of regions
[mask_cmtk,cmtk_labels_mask_cmtk] = mrtrix3_to_cmtk_order(mask,in_labels);

up_left = mask_cmtk(1:N/2,1:N/2);
down_right = mask_cmtk(N/2+(1:N/2),N/2+(1:N/2));
up_left(isnan(up_left)) = 0;
down_right(isnan(down_right)) = 0;
        
mask_L_and_R_cmtk = up_left .* down_right; % make a new mask by multiplying both sides. Do not generate binary mask because then cannot be reordered
[mask_L_and_R,twice_converted_labels_mask] = cmtk_to_mrtrix3_order(mask_L_and_R_cmtk,cmtk_labels_mask_cmtk);        
mask_L_and_R(isnan(mask_L_and_R)) = 0; % turn NaN to 0 and then turn to logical matrix
mask_L_and_R = logical(mask_L_and_R);
        
        
for i = 1:length(file_list)
    file_name = file_list(i).name; % Get the name of the file
    file_path = fullfile(folder_path, file_name); % Construct the full file path
    disp(['Loaded file ', file_name]);
    
    if ~strcmp(file_name, 'mask.csv')
        mismatch = dlmread(file_path); % Load the data from the file

        N = size(mismatch,1);


        % 2. need to save L_gt_R for each subject V

        % change the order of regions from mrtrix3 to cmtk, because in the cmtk
        % ordering region i and 42+i are homologous, and it makes calculation easier
        [mismatch_cmtk,cmtk_labels_mismatch] = mrtrix3_to_cmtk_order(mismatch,in_labels);

        mismatches_cmtk(i,:,:) = mismatch_cmtk;

        % calculate (mismatch_left - mismatch_right) for the subject and save it to a .txt file    
        L_gt_R_cmtk = mismatch_cmtk(1:N/2,1:N/2) - mismatch_cmtk(N/2+(1:N/2),N/2+(1:N/2)); % automatically takes care of NaN, because operation between two NaNs is NaN    
        % BETTER NOT SAVE, BECAUSE ROSENTHAL ORDER - connectome_write([PREFIX '_' subj_name '.txt'],left_gt_right); %mirror_mat((-a) - (-b)));    

        % go back to original order, as to not confuse users
        [L_gt_R,twice_converted_labels_L_gt_R] = cmtk_to_mrtrix3_order(L_gt_R_cmtk,cmtk_labels_mismatch);

        % mark cells that were excluded
        L_gt_R(~mask_L_and_R) = -99;
        connectome_write([output_folder filesep 'individual_connectomes' filesep 'LminusR-' file_name],L_gt_R);

        L_gt_Rs(i,:,:) = L_gt_R;
    end
end

% save mean L_gt_R
mean_L_gt_R = squeeze(mean(L_gt_Rs));

connectome_write([output_folder filesep 'mask_L_and_R.csv'],mask_L_and_R);
mean_L_gt_R(~mask_L_and_R) = -99;
connectome_write([output_folder filesep 'LminusR-mean.csv'],mean_L_gt_R);

% save std L_gt_R
sd_L_gt_R = squeeze(std(L_gt_Rs));

% 3. Calculate sign files that I will need to use later on. Initially in
% cmtk, and then convert to mrtrix3
mean_mismatch_cmtk = squeeze(mean(mismatches_cmtk));
sign_L_cmtk = mean_mismatch_cmtk(1:N/2,1:N/2)./abs(mean_mismatch_cmtk(1:N/2,1:N/2));
sign_R_cmtk = mean_mismatch_cmtk(N/2+(1:N/2),N/2+(1:N/2))./abs(mean_mismatch_cmtk(N/2+(1:N/2),N/2+(1:N/2)));

[sign_L, twice_converted_labels_sign_L] = cmtk_to_mrtrix3_order(sign_L_cmtk,cmtk_labels_mismatch);
[sign_R, twice_converted_labels_sign_R] = cmtk_to_mrtrix3_order(sign_R_cmtk,cmtk_labels_mismatch);

% save ps
[~,ps] = ttest(squeeze(L_gt_Rs)); % uncorrected ps, as correction done in another function
ps = squeeze(ps); 

strs = {'pos_L_gt_R',
'neg_L_st_R',
'pos_L_st_R',
'neg_L_gt_R',
'flip_L_gt_R',
'flip_L_st_R',
'all'};

% 4. Separate ps files + Separate label files (to compare with paper)
pvalue{find(ismember(strs,'pos_L_gt_R'))} = (1 - ps) .* (mean_L_gt_R > 0) .* (sign_L == 1 & sign_R == 1); 
pvalue{find(ismember(strs,'neg_L_st_R'))} = (1 - ps) .* (mean_L_gt_R < 0) .* (sign_L == -1 & sign_R == -1); 
pvalue{find(ismember(strs,'pos_L_st_R'))} = (1 - ps) .* (mean_L_gt_R < 0) .* (sign_L == 1 & sign_R == 1);
pvalue{find(ismember(strs,'neg_L_gt_R'))} = (1 - ps) .* (mean_L_gt_R > 0) .* (sign_L == -1 & sign_R == -1);
pvalue{find(ismember(strs,'flip_L_gt_R'))} = (1 - ps) .* (sign_L == 1 & sign_R == -1); 
pvalue{find(ismember(strs,'flip_L_st_R'))} = (1 - ps) .* (sign_L == -1 & sign_R == 1); 
pvalue{find(ismember(strs,'all'))} = (1-ps);

% calculate benferroni correction on the upper triangle only
correction = ((N/2 * N/2) / 2 - (N/2/2));
alpha = 0.05 / correction;

all_sig = zeros(N/2, N/2);
overlap_sig = ones(N/2,N/2);

for str_i=1:7
    connectome_write([output_folder filesep 'pvalue_' strs{str_i} '.csv'],pvalue{str_i},'precision',10);
    sig{str_i} = pvalue{str_i} > 1-alpha;
    connectome_write([output_folder filesep 'sig_' strs{str_i} '.csv'],sig{str_i},'precision',10);

    
    fid = fopen([output_folder filesep 'sig_' strs{str_i} '_labels.txt'],'wt');
    [pair_is,pair_js] = find(sig{str_i});
    for pair_num = 1:length(pair_is)
        if mask_L_and_R(pair_is(pair_num),pair_js(pair_num)) % only print the ones in the mask
            fprintf(fid,'%s - %s\n',in_labels{pair_is(pair_num)},in_labels{pair_js(pair_num)});
        end
    end
    fclose(fid);

    if str_i~=7
        all_sig = all_sig | sig{str_i};
        overlap_sig = overlap_sig & sig{str_i};
    end
end

% check if significance maps sum. Only test items above diagonal
if ~any(triu(all_sig,1) ~= triu(sig{7},1)) 
    error('Missing values in significance maps. Connect developers');
end

% check no overlap between significance maps
if any(overlap_sig(mask_L_and_R))
    error('Overlap between significance maps. Connect developers');
end


