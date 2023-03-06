function calculate_laterality(input_folder,output_folder)
%
%      calulate_laterality
%      ===================
%
%      Note: difference between hemispheres will be only caclulated to
%      connections which are within the mask.csv in both hemispheres.
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

%%% This is the main script for the SCFC metric, previously it was called
%%% "show_transformations.m"





%%% HELENA - files names need to be SCXXX.csv for structural and
%%% FCXXX.csv for funcitonal. XXX is the number of the subject

%%% HELENA - notice that I clear all variables every time that I run
%%% script! It's safer
%clear all; 
dbstop if error;

%%% HELENA - do these steps:
%%% 1) run with 'SAVE_AVG = true' to save averages 
%%% 2) use curvefit to obtain a,b,c
%%% 3) update a,b,c below, and run again with 'SAVE_AVG = false' & 'SHOW_FIRST_SUBJS = true' in order to look on the transformations
%%% 4) run with 'SAVE_AVG = false' & 'SHOW_FIRST_SUBJS = false' & 'CALCULATE_SCFC = true' to calculate metric
%%%    It uses SUBJ_1_plot_connectome_90_mismatch_FS_reg.m
%%% 5) You need to run reorganize_all.m
%%%
%%%    Output related to mismatch will be: 
%%%    4_RUNS_REG_mismatch.csv - average mismatch of each conections
%%%    BIN_4_RUNS_REG_pvalue_ANY.csv - all connections where there is a
%%%    significant difference between left and right, Benferroni corrected
%%%
%%%    There is one file for each one of the 6 options shown in article:
%%%    BIN_MRtrix_4_RUNS_REG_pos_L_gt_R.csv
%%%    BIN_MRtrix_4_RUNS_REG_neg_L_st_R.csv
%%%    BIN_MRtrix_4_RUNS_REG_pos_L_st_R.csv
%%%    BIN_MRtrix_4_RUNS_REG_neg_L_gt_R.csv
%%%    BIN_MRtrix_4_RUNS_REG_flip_L_gt_R.csv
%%%    BIN_MRtrix_4_RUNS_REG_flip_L_st_R.csv
%%%
%%%    For each of these files, there is also a txt file with only half the
%%%    matrix (in Rosenthal ordering, if Oren's data):
%%%    4_RUNS_REG_pos_L_gt_R.txt - half size  % Rosenthal ordering
%%%
%%%    Output related to FC assymetry will be:
%%%    4_RUNS_REG_FC_L_gt_R.csv - 
%%%    4_RUNS_REG_FC_L_st_R.csv
%%%



%     DO_EXTRA = false;
%     %%% set these two if want to see specific connection
%     %roi1 = 
%     %roi2 =
%     
%     rs = [];
%     intersects = [];
%     slopes = [];
%     mismatches = [];
%     L_gt_Rs = [];
%     FC_L_gt_Rs = [];
%     order_i = 0;
% 
%     for i=1:length(subjs)
%         %for i= [4, 14, 19, 46]  %[*13, 14 , *15, 19] %35]   %19 - both have large mismatch %length(subjs)
%         
%         order_i = order_i + 1;
%         
%         % calculate SCFC for each subject
% 
%         [r,intersect,slope,mismatch,L_gt_R,FC_L_gt_R] = calculate_mismatch_for_subject(sc_all(i,:,:),fc_all(i,:,:),a,b,c,FC_PREFIX,subjs{i},order_i);
% 
%         rs = [rs r];
%         intersects = [intersects intersect];
%         slopes = [slopes slope];
%         mismatches(i,:,:,:) = mismatch;
%         L_gt_Rs(i,:,:,:) = L_gt_R;
%         FC_L_gt_Rs(i,:,:,:) = FC_L_gt_R;
%         clear r intersect slope mismatch L_gt_R;
%         order_i = order_i + 1;
%     end

% try to create output directory - so if fails, not continue anymore
if ~exist(output_folder, 'dir')
       mkdir(output_folder);
end

if ~exist([output_folder filesep 'individual_connectomes'], 'dir')
       mkdir([output_folder filesep 'individual_connectomes']);
end

%%%% JUST SO CAN VERIFY THAT CONVERSTION MOVE LABELS TO RIGHT PLACE %%% 
BASE = ['.' filesep];
fid = fopen([BASE 'MRtrix_labels_Desikan-Killiany.txt']);
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

% % different classes
% pos_L_gt_R = 1;
% neg_L_st_R = 2;
% pos_L_st_R = 3;
% neg_L_gt_R = 4;
% flip_L_gt_R = 5;
% flip_L_st_R = 6;
% all = 7;

% strings to use for filenames
% str{pos_L_gt_R} = 'pos_L_gt_R';
% str{neg_L_st_R} = 'neg_L_st_R'; 
% str{pos_L_st_R} = 'pos_L_st_R';
% str{neg_L_gt_R} = 'neg_L_gt_R';
% str{flip_L_gt_R} = 'flip_L_gt_R'; 
% str{flip_L_st_R} = 'flip_L_st_R'; 

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

% ====

% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_pos_L_gt_R.csv'],(1 - ps) .* (mean_L_gt_R > 0) .* (sign_L == 1 & sign_R == 1),'precision',10); 
% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_neg_L_st_R.csv'],(1 - ps) .* (mean_L_gt_R < 0) .* (sign_L == -1 & sign_R == -1),'precision',10); 
% 
% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_pos_L_st_R.csv'],(1 - ps) .* (mean_L_gt_R < 0) .* (sign_L == 1 & sign_R == 1),'precision',10); 
% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_neg_L_gt_R.csv'],(1 - ps) .* (mean_L_gt_R > 0) .* (sign_L == -1 & sign_R == -1),'precision',10); 
% 
% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_flip_L_gt_R.csv'],(1 - ps) .* (sign_L == 1 & sign_R == -1),'precision',10); 
% connectome_write([output_folder filesep '4_RUNS_REG_pvalue_flip_L_st_R.csv'],(1 - ps) .* (sign_L == -1 & sign_R == 1),'precision',10); 

ghjhjhg

% =======================================================================

% 5. Save all files to disc


    %%% TO SEE MISMATCHES FOR SPECIFIC CONNECTIONS USE
    % Notice that roi1 < roi2 in order, and that is according to Rosenthal
    % ordering
    if exist('roi1') & exist('roi12')
        figure; plot(mismatches(:,roi1,roi2) - mismatches(:,roi1+42,roi2+42),'.')
    end
    %%%%%

    figure;
    subplot(131);
    hist(rs);
    title('r');

    subplot(132);
    hist(intersects);
    title('intersect');

    subplot(133);
    hist(slopes);
    title('slope');

    % calcualate mean values and stats


    sd_L_gt_R = squeeze(std(L_gt_Rs));

    mean_FC_L_gt_R = squeeze(mean(FC_L_gt_Rs));
    [stub,FC_ps] = ttest(squeeze(FC_L_gt_Rs));
    FC_ps = squeeze(FC_ps); 
    effect_FC_L_gt_R = squeeze(mean(FC_L_gt_Rs) ./ std(FC_L_gt_Rs));

    mean_L_gt_R = squeeze(mean(L_gt_Rs));
    [stub,ps] = ttest(squeeze(L_gt_Rs)); % uncorrected ps, as correction done in another function
    ps = squeeze(ps); 
    effect_L_gt_R = squeeze(mean(L_gt_Rs) ./ std(L_gt_Rs));

    %ttest_L_gt_R = squeeze(ttest_L_gt_R);

    % corrected not saved to file, but used for visualization on diagonal later
    % in this function
    corrected_L_gt_R = ps < (0.05 / ((N/2*N/2)/2 - N/2/2)); % only if not nan, and p value < 0.05/number_of_comparisons



    
    
    %pos_ttest_L_gt_R = corrected_L_gt_R .* (effect_L_gt_R > 0);
    %neg_ttest_L_gt_R = corrected_L_gt_R .* (effect_L_gt_R < 0);
    %

    % save mismatches
    a = mean_mismatch;
        % nan upper part of functional and structual matrices
        % and transfer value to other half if necessary
        for i=1:N/2
            for j=i:N/2
                if a(j,i) > 0
                    a(i,j) = a(j,i);
                end
                a(j,i) = nan;

                if a(j,i) > 0
                    a(i,j) = a(j,i);
                end
                a(j,i) = nan;
            end
        end
    [MRtrix_mean_mismatch,stub] = prepare_sc_fc_MRtrix(a,zeros(N,N));    
    MRtrix_mean_mismatch(isnan(MRtrix_mean_mismatch)) = 0;
    connectome_write('4_RUNS_REG_mismatch.csv',MRtrix_mean_mismatch);
    clear a;
    clear MRtrix_mean_mismatch;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% From here on it is optional
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if DO_EXTRA
    % save difference (L gt R)
        b = [mean_L_gt_R zeros(N/2,N/2);
        zeros(N/2,N/2) mean_L_gt_R;];     
        % nan upper part of functional and structual matrices
        % and transfer value to other half if necessary
        for i=1:N/2
            for j=i:N/2
                if b(j,i) > 0
                    b(i,j) = b(j,i);
                end
                b(j,i) = nan;

                if b(j,i) > 0
                    b(i,j) = b(j,i);
                end
                b(j,i) = nan;
            end
        end
        [MRtrix_mean_L_gt_R,stub] = prepare_sc_fc_MRtrix(b,zeros(N,N));    
        MRtrix_mean_L_gt_R(isnan(MRtrix_mean_L_gt_R)) = 0;
        connectome_write('4_RUNS_REG_L_gt_R.csv',MRtrix_mean_L_gt_R);
        clear b;
        clear MRtrix_mean_L_gt_R;


        % save FC laterality
        c = [mean_FC_L_gt_R zeros(N/2,N/2);
            zeros(N/2,N/2) mean_FC_L_gt_R;];     
            % nan upper part of functional and structual matrices
            % and transfer value to other half if necessary
            for i=1:N/2
                for j=i:N/2
                    if c(j,i) > 0
                        c(i,j) = b(j,i);
                    end
                    c(j,i) = nan;

                    if c(j,i) > 0
                        c(i,j) = c(j,i);
                    end
                    c(j,i) = nan;
                end
            end
        [MRtrix_mean_FC_L_gt_R,stub] = prepare_sc_fc_MRtrix(c,zeros(N,N));    
        MRtrix_mean_FC_L_gt_R(isnan(MRtrix_mean_FC_L_gt_R)) = 0;
        connectome_write('4_RUNS_REG_FC_L_gt_R.csv',MRtrix_mean_FC_L_gt_R);

        % plot all data points
        figure; plot(mean_mismatch(N/2+(1:N/2),N/2+(1:N/2)),mean_mismatch(1:N/2,1:N/2),'*k'); % notice left hemi is Y axis
        hold on; plot([-0.5 0.5],[-0.5 0.5],'k'); 
        axis([-0.4 0.4 -0.4 0.4]);    

        % plot SD or SE and write effect size
        for i=1:N/2
            for j=1:N/2
                if corrected_L_gt_R(i,j) % only if not nan, and p value < 0.5/number_of_comparisons
                    color = 'r';
                    text(mean_mismatch(N/2+i,N/2+j),mean_mismatch(i,j),num2str(effect_L_gt_R(i,j)),'Color','r');
                else
                    color = 'k';
                end
                SE_FACTOR = 50;
                plot([mean_mismatch(N/2+i,N/2+j) mean_mismatch(N/2+i,N/2+j)], [mean_mismatch(i,j)-sd_L_gt_R(i,j)./sqrt(SE_FACTOR) mean_mismatch(i,j)+sd_L_gt_R(i,j)./sqrt(SE_FACTOR)],color);
            end
        end

        ylabel('Left mismatch');
        xlabel('Right mismatch');
    end        
end
