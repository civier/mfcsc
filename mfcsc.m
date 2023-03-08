function mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dir)
%%% 
%%%
%%%     mfcsc -- Mismatch between Functional Connectivity and Structural Connectivity  
%%%     =============================================================================
%%%
%%%     JOURNAL ARTICLE
%%%     
%%%          Civier O, Sourty M, Calamante F (2023) MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, 
%%%          and its application for detecting hemispheric functional specialisations. Scientific Reports https://doi.org/10.1038/s41598-022-17213-z
%%%
%%%     USAGE
%%%
%%%         From Neurodesk (no Matlab license required!):
%%%             See https://osf.io/d7j9n/
%%%
%%%         From Matlab IDE:
%%%             mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dir)         
%%%
%%%         From the command-line:
%%%             Change to the folder with the code, and then enter the following command -
%%%             matlab -batch "mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dor)"
%%%
%%%         The 4 arguments in caps are mandatory. 
%%%
%%%         IMPORTANT: before running mfcsc, it is recommended to open FC_SC_LIST 
%%%         and visually check that the two columns have matching
%%%         participant connectomes.
%%%
%%%     PREREQUISITES:
%%%
%%%         Matlab's Curve Fitting toolbox
%%%
%%%     DESCRIPTION
%%%
%%%         mfcsc receive pairs of functional and structural connectivity
%%%         matrices, one pair per participant. For each participant, mfsc
%%%         calculates a connectivity matrix that gives the mFCSC metric
%%%         value for every connection. 
%%%
%%%         The mFCSC metric can be calculated either for ipsilateral or 
%%%         contralateral connectionsis (but not both), and is ill-defined 
%%%         for some connections. The 'mask.csv' binary matrix mask
%%%         indicates which mFCSC metrics should be included in further analysis. 
%%%         The cells that should not be consulted further are set to -99.
%%%
%%%     INPUT
%%%
%%%         The input connectivity matrices should be comma seperated values that describe symmetrical 
%%%         connectomes (i.e., connections do not have directionality). 
%%%         Only the upper right triangle of the connectivity matrices is
%%%         consulated. This also excludes the main diagonal.
%%%
%%%         The connectivity matrices must have an even number of regions N: regions numbers 
%%%         1 to N/2 for one hemisphere, and regions number N/2+1 to N for the other.
%%%         Note that the order of regions for each hemisphere may be different, i.e., 
%%%         regions number i and N/2+i do not have to be homologous.
%%%
%%%     OUTPUT
%%%
%%%         IMPORTANT: before consulting the output, ensure that the file1 and file2 parts 
%%%         of the filenames below represent matching connectivity
%%%         matrices!
%%%         
%%%         The output of mfcsc consists of one connectome file for each participant:
%%%
%%%                 mFCSC-file1-file2-masked.csv - connectome of mFCSC values for the
%%%                                                participant whose FC and SC
%%%                                                connectomes are stored in file1 and file2
%%%                                                respectively (excluding
%%%                                                the files extensions)
%%%
%%%         The main output of mfcsc consists of the file:
%%%
%%%                 mask-final.csv - final mask indicating the connections to which mFCSC is calculated
%%%

%%%
%%%         There are also several misc files in the group_connectomes subdir:
%%%
%%%             transformed_SC_avg - the average transfered SC connectome
%%%             FC_avg - the average FC connectome
%%%             mask-direct_SC_is_shortest_path.csv - mask of connections in which the path length of the direct connection (1/transformed_SC) is shorter than any other indirect path between the two regions
%%%             SC_avg - the average SC connectome
%%%
%%%     MANDATORY ARGUMENTS
%%%
%%%         FC_SC_LIST
%%%
%%%         Path to a tab-separated file with two columns.
%%%         The first column lists the files with the connectivity matrices of the FC connectomes
%%%         The second column lists the connectivity matrix files for the matching SC connectomes
%%%
%%%         An easy method to generate the FC_SC_LIST on Linux or MacOS is to:
%%%         1) include a matching participant ID in the filenames of both FC and SC connectomes 
%%%            (plus optional fixed prefixes and suffixes for each modality).
%%%         2) put all FC connectomes in one folder and all SC connectomes in another,
%%%            ensuring that there are no other files in these folders.
%%%         3) change the working directory to the folder with FC connectomes, and run:
%%%             ls | sort -n > /tmp/fc_list
%%%         4) Change the working directory to the folder with the SC connectomes, and run:
%%%             ls | sort -n > /tmp/sc_list
%%%         5) Run:
%%%             paste /tmp/fc_list /tmp/sc_list > path_to_filename
%%%         6) Provide path_to_filename as the FC_SC_LIST argument
%%%
%%%         FC_INPUT_DIR
%%%
%%%         The directory containing the files with the FC connectivity matrices
%%%
%%%         SC_INPUT_DIR
%%%
%%%         The directory containing the files with the SC connectivity matrices
%%%
%%%         OUPTUT_DIR
%%%
%%%         The output directory where the mask and mFCSC files are to be
%%%         saved
%%%
%%%     OPTIONAL ARGUMENTS
%%%
%%%         not_in_mask_value
%%%
%%%         Value that will be assigned to cells in the output matrices that are
%%%         not in the mask. By deafult it is set to -99 to make sure
%%%         people do not report the values in these cells.
%%%         Can be set to 0 to prevent this value from affecting color scaling of plots.    
%%%
%%%         is_contra
%%%         
%%%         By default, mFCSC is calculated for ipsilateral connections
%%%         (is_contra = false or omitted)
%%%         set is_contra to true to calculate mFCSC for contralateral connections instead.
%%%
%%%         is_keep_neg_fc
%%%
%%%         By default, mFCSC removes cells that have negative mean FC before 
%%%         fitting the model used to transform SC
%%%         (is_remove_negative_fc = false or omitted)
%%%         set is_keep_negative_fc to true to keep them
%%%
%%%         is_symmetrical
%%%
%%%         By default, only the upper right triangle of the output matrices is
%%%         populated, with the lower triangle being zeroed out.
%%%         (is_symmetrical = false or omitted)
%%%         set is_symmetrical to true to save symmetrical matrices instead by
%%%         mirroring the upper right triangle into the bottom left triangle
%%%         In both cases, the main diagonal is zeroed out.
%%%
%%%         is_figures
%%%
%%%         By default, do not show figures.
%%%         (is_figures = false or omitted)
%%%         set is_figures to true to print verbose figures with information for QC and debg.
%%%         Not tested. Use at your own risk!
%%%
%%%         bct_dir
%%%
%%%         If specified, Matlab will look for the Brain Connectome Toolbox (BCT)
%%%         in this directory instead of the BCT version supplied with
%%%         MFCSC (2017/01/15). BCT is required in order to calculate the mask.
%%%         This argument is not available in the Neurodesk version, where only the supplied
%%%         BCT can be used.
%%%
%%%     NOTE
%%%         
%%%         Cells of individual SC connectivity matrices that have the value of 0 are
%%%         not transfromed well using the model.
%%%         In most cases these cells will not be included in the mask because the
%%%         direct SC in the average SC connectome is not the shortest path;
%%%         however, in case they are within the mask after all, a warning will
%%%         be issued and they will be assigned the mFCSC value of -999.
%%%         In this case, one approach is to manually exclude these cells from the
%%%         mask by editing 'mask-final.csv' and using the ammended mask in
%%%         further analyses.
%%%
%%%
%%%      DEVELOPER
%%%
%%%         Oren Civier (orenciv@gmail.com)
%%%         https://www.swinburne.edu.au/research/our-research/access-our-research/find-a-researcher-or-supervisor/researcher-profile/?id=ocivier
%%%
%%%      CITATIONS
%%%
%%%         When using mfcsc, authors should cite:
%%%
%%%             Civier O, Sourty M, Calamante F (2023) MFCSC: Novel method to calculate mismatch between functional and structural brain connectomes, 
%%%             and its application for detecting hemispheric functional specialisations. Scientific Reports https://doi.org/10.1038/s41598-022-17213-z
%%%
%%%             Rubinov M, Sporns O (2010) Complex network measures of brain connectivity: Uses and interpretations. NeuroImage 52:1059-69.%
%%%
%%%      ACKNOWLEDGMENTS
%%%
%%%         National Health and Medical Research Council of Australia (grant numbers APP1091593 andAPP1117724)
%%%         Australian Research Council (grant number DP170101815)
%%%         National Imaging Facility (NIF), a National Collaborative Research Infrastructure Strategy (NCRIS) capability at Swinburne Neuroimaging, 
%%%         Swinburne University of Technology.
%%%


% To compile for MCR, run:
% cd /home/ocivier/mfcsc
% mcc -m mfcsc.m -o mfcsc -R -nojvm

% To test MCR, run:
% cd /fred/oz120/oren/mfcsc_input/Data
% !/home/ocivier/mfcsc/mfcsc FC_SC_list.txt Functional_Connectomes_From_SCFC_FS Structural_Connectomes_From_SCFC_FS ../output

% OC - Mac/Linux example - remove from final version
% mfcsc('/Users/ocivier/Dropbox/Documents/fernando/fcsc/Data/FC_SC_list.txt', '/Users/ocivier/Dropbox/Documents/fernando/FCSC/Data/Functional_Connectomes_From_SCFC_FS', '/Users/ocivier/Dropbox/Documents/fernando/FCSC/Data/Structural_Connectomes_From_SCFC_FS', '/Users/ocivier/Dropbox/Documents/fernando/FCSC/Code/MFCSC_code/Output',false,false,false,false,'/Users/ocivier/Dropbox/Documents/fernando/phase between structural and functional/Code/2017_01_15_BCT');
% OC - Windows example - remove from final version
% mfcsc('C:\Users\ocivier\mfcsc-data\FC_SC_list.txt', 'C:\Users\ocivier\mfcsc-data\Functional_Connectomes_From_SCFC_FS', 'C:\Users\ocivier\mfcsc-data\Structural_Connectomes_From_SCFC_FS', 'C:\Users\ocivier\mfcsc-output',false,false,false,false)

% OC - remove from final version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The values below are for 50 HPC datasets, 4 runs, no z-score and no fisher, using power 2, LAR, Levenberg
%a =      0.4114;
%b =     0.09267;
%c =     -0.3789;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OC - remove from final version
%dbstop if error;

INCLUDED_BCT_SUBDIR = '2017_01_15_BCT';

% OC - I think the below can be removed
% add current folder to find other functions in package
% addpath '.'

if ~exist('FC_SC_LIST') || ~exist('FC_INPUT_DIR') || ~exist('SC_INPUT_DIR') || ~exist('OUTPUT_DIR')
    fprintf('\nmfcsc: ERROR: not enough input arguments!\n');
    fprintf('\nUSAGE -\n');
    fprintf('Matlab''s IDE (MATLAB license required): mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dir)\n');
    fprintf('Command line (MATLAB license required): matlab -batch "mfcsc(FC_SC_LIST,FC_INPUT_DIR,SC_INPUT_DIR,OUTPUT_DIR,not_in_mask_value,is_contra,is_keep_neg_fc,is_symmetrical,is_figures,bct_dir)"\n');
    fprintf('Neurodesk (no MATLAB license required): mfcsc FC_SC_LIST FC_INPUT_DIR SC_INPUT_DIR OUTPUT_DIR not_in_mask_value is_contra is_keep_neg_fc is_symmetrical is_figures\n');
    fprintf('\nNotice first 4 arguments (in caps) are mandatory\n'); 
    return;
end

if ~exist('not_in_mask_value')
    not_in_mask_value = '-99';
end
global exclude_value;
exclude_value = str2num(not_in_mask_value);
% set flags to false if omitted
if ~exist('is_contra')
    is_contra = false;
elseif isdeployed
    is_contra = str2logical(is_contra);
end
if ~exist('is_keep_neg_fc')
    is_keep_neg_fc = false;
elseif isdeployed
    is_keep_neg_fc = str2logical(is_keep_neg_fc);
end
if ~exist('is_symmetrical')
    is_symmetrical = false;
elseif isdeployed
    is_symmetrical = str2logical(is_symmetrical);
end
if ~exist('is_figures')
    is_figures = false;
elseif isdeployed
    is_figures = str2logical(is_figures);
end
if ~exist('bct_dir')
    % Get the folder of the mfcsc.m script
    [folder,~,~] = fileparts(mfilename('fullpath'));
    bct_dir = [folder filesep INCLUDED_BCT_SUBDIR];
%     if ~exist(bct_dir,'dir')
%         'folder does not exist'
%         BCT_DIR
%         BCT_DIR = 'built-in';
%     end
end

% make flag global, as it is used by the connectome_write function
global is_symmetrical_global;
is_symmetrical_global = is_symmetrical;

% print all input arguments
fprintf('\nInput arguments (or if not provided, defaults)\n\n');
fprintf('FC_SC_LIST:\t%s\n',FC_SC_LIST); 
fprintf('FC_INPUT_DIR:\t%s\n',FC_INPUT_DIR); 
fprintf('SC_INPUT_DIR:\t%s\n',SC_INPUT_DIR); 
fprintf('OUTPUT_DIR:\t%s\n',OUTPUT_DIR); 
fprintf('not_in_mask_value:\t%s\n',not_in_mask_value);
LogicalStr = {'false', 'true'};
fprintf('is_contra:\t%s\n',LogicalStr{is_contra+1}); 
fprintf('is_keep_neg_fc:\t%s\n',LogicalStr{is_keep_neg_fc+1}); 
fprintf('is_symmetrical:\t%s\n',LogicalStr{is_symmetrical+1}); 
fprintf('is_figures:\t%s\n',LogicalStr{is_figures+1}); 
if ~isdeployed % BCT_DIR is built in if code is deployed
    fprintf('bct_dir:\t%s\n',bct_dir); 
end
fprintf('\n'); 

% add the BCT to Matlab's path
if ~isdeployed % no need to add BCT_DIR if code is deployed
    if exist(bct_dir,'dir')
        addpath(bct_dir);
    else
        fprintf('\nmfcsc: Brain Connetivity Toolbox is not in bct_dir!\n');
        return;
    end
end

% enable calculating mask. Disabling it not tested; do at your own risk!
CALCULATE_MASK = true;

% enable the calculation of the mFCSC matrices. Can be disabled if one is interested in the mask only       
CALCULATE_SCFC = true;

% save the average matrices. Can be discabled if desired
SAVE_AVERAGES = true;

GROUP_FOLDER = 'group_connectomes';

%%% This should point to the BCT package
%%% For compatability, it is best to use the version with which MFCSC was
%%% tested (15/1/2017). It is provided as part of the MFCSC package

% Set general constants
MINIMUM_TO_SHOW = -1000;
PLOT_NEW = false;
PLOTS = 3;

% Set plots constants
LEFT = -0.25;
RIGHT = 1;
BINS = round((RIGHT - LEFT) / 0.1 * 2);
HEIGHT = 600;

SUBJ_TO_SHOW = 5;
%% N = 84;   %% DON'T USE ANYMORE, AS MATRICES MIGHT BE EXTENDED
SAVE_FC = false; %%% also need to remove comments from section
SHOW_FC = true;

% go over subjects
clear fc_all;
clear sc_all;

% OC - should check if can allow either tab or space (and multiples of
% them), basically like AWK
opts = detectImportOptions(FC_SC_LIST,'Delimiter','\t');
opts.VariableNames = {'FC','SC'};
opts.DataLines= [1 Inf];
filenames = readtable(FC_SC_LIST,opts);

numSubjs = height(filenames);

fprintf(['List includes ' num2str(numSubjs) ' pairs of FC and SC connectomes\n']);

for i=1:numSubjs
    % read structural and functional matrices
    sc = dlmread([SC_INPUT_DIR filesep filenames.SC{i}]);
    %sc = remove_subcortical_and_extend(sc);

    fc = dlmread([FC_INPUT_DIR filesep filenames.FC{i}]);
        
        %%% OR INSTEAD, CODE TO TURN MAT FILES INTO TXT FILES
        %[~,filename,~]=fileparts(filenames.FC{i});
        %load([FC_INPUT_DIR '/' filename '.mat']);
        % fc = sFC_mat;
        % OC - writematrix is recommended for compatability, but since cannot
        % specify precision and default is less than what I need, I do not use it
        %dlmwrite([FC_INPUT_DIR '/' filenames.FC{i} '.txt'],fc,'precision', 32);

    %load(['./' subjs{i} '/' subjs{i} '_sFC_meanFC.mat']);
    %fc = meanFC;

    fc = remove_lower_half(fc);
    sc = remove_lower_half(sc);

    % add to a big matrix that summarizes everything
    fc_all(i,:,:) = fc;
    sc_all(i,:,:) = sc;

end

N = size(fc_all,2);
if N ~= size(sc_all,2)
    error('FC and SC have different sizes');
end

% calculate mean matrices
fc_avg = squeeze(mean(fc_all));
sc_avg = squeeze(mean(sc_all));

if ~exist(OUTPUT_DIR, 'dir')
       mkdir(OUTPUT_DIR);
end
if ~exist([OUTPUT_DIR filesep GROUP_FOLDER], 'dir')
       mkdir([OUTPUT_DIR filesep GROUP_FOLDER]);
end


%if SAVE_AVG

% OC - Save averages, in original matrix form
% OC - check if really necessary to save these?
connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'FC_average-unmasked.csv'],fc_avg);
connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'SC_average-unmasked.csv'],sc_avg);    

fprintf('The connectomes have %d regions (including both hemispheres)\n',N);
fprintf('Each connectome has %d connections (including both ipsilateral and contralateral connections)\n\n',N*N/2-N/2);

num_of_removed_connections = 0;
if ~is_keep_neg_fc
    I = fc_avg < 0;
    num_of_removed_connections = sum(sum(I));
    fprintf('Out of %d connections, in %d connections FC is not negative\n\n',N*N/2-N/2,N*N/2-N/2-num_of_removed_connections);
    fc_avg(I) = I(I) * NaN;
    sc_avg(I) = I(I) * NaN;
end

% turn into vectors, and sort to prepare for createFit_fcsc_LAR
sc_avg_vec_sort = sort(sc_avg(:));
fc_avg_vec_sort = sort(fc_avg(:));

% remove NaNs to prevent warning message in createFit_fcsc_LAR
sc_avg_vec_sort = sc_avg_vec_sort(~isnan(sc_avg_vec_sort));
fc_avg_vec_sort = fc_avg_vec_sort(~isnan(fc_avg_vec_sort));

% Calculate model using cuvefit toolbox
[fcscfit, ~] = createFit_fcsc_LAR(sc_avg_vec_sort, fc_avg_vec_sort,is_figures);

fprintf('Transform SC using the below model (transformed_SC = fcscfit(SC)):\n');
disp(fcscfit);
fprintf('\n');

if SAVE_AVERAGES
    
    sc_avg_nonzero = sc_avg;
    sc_avg_nonzero(sc_avg == 0) = nan;

    transformed_sc_avg = (sc_avg_nonzero.^fcscfit.b)*fcscfit.a+fcscfit.c;

    % If IGNORE_ZERO_SC then turn cells that are not in the mask,
    % but still NaN into -999, otherwise turn them into 0
    if any(isnan(transformed_sc_avg))
        fprintf('mfsc: WARNING. Average SC has connections that are zero. Set them in transformed average SC to -999\n\n');
        transformed_sc_avg(isnan(transformed_sc_avg)) = -999;
    end

    connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'transformed_SC_average-unmasked.csv'],transformed_sc_avg);

end    
       
if CALCULATE_MASK
    % input must be a mat! Notice that retursn a symmetric mat!
    direct_is_shortest_mask_sym = calculate_direct_is_shortest_mask_sym(sc_avg,fcscfit.a,fcscfit.b,fcscfit.c,is_figures);
    % print out the number of connections in the mask assuming it's
    % symmeric
    fprintf('Out of %d connections, in %d connections the path length of the direct structural connection (1/transformed_SC) is shorter than any other indirect path between the two regions\n\n',N*N/2-N/2-num_of_removed_connections,sum(sum(direct_is_shortest_mask_sym))/2);
    connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'mask-direct_SC_is_shortest_path.csv'],direct_is_shortest_mask_sym);
    
    N = length(direct_is_shortest_mask_sym);
    
    connection_type_mask_sym = false(N,N);
    if is_contra
        % turn on contralateral in upper triangle
        connection_type_mask_sym(1:N/2,N/2+1:N) = true;
        % turn on contralateral in lower triangle
        connection_type_mask_sym(N/2+1:N,1:N/2) = true;
        postfix = 'contralateral';
    else
        % turn on ipsilateral (affects both sides of the triangle)
        connection_type_mask_sym(1:N/2,1:N/2) = true;
        connection_type_mask_sym(N/2+1:N,N/2+1:N) = true;
        postfix = 'ipsilateral';
    end

    connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'mask-' postfix '.csv'],connection_type_mask_sym);    
    
    combined_mask_sym = direct_is_shortest_mask_sym & connection_type_mask_sym;
    
    connectome_write([OUTPUT_DIR filesep 'mask.csv'],combined_mask_sym);
    
    fprintf('Out of %d connections, only calculating mFCSC for the %d %s connections\n\n',sum(sum(direct_is_shortest_mask_sym))/2,sum(sum(combined_mask_sym))/2,postfix);
else
    combined_mask_sym = ones(N,N);
end
 
        % OC - need to find how to dump it into an HTML or image file?
        
% Not tested. Use at your own risk
if is_figures

    figure;
    
    for i=1:SUBJ_TO_SHOW

        subplot(PLOTS,SUBJ_TO_SHOW,i);
        if i==1
            title('avg (ex. FC negatives)');
            fc = fc_avg;
            sc = sc_avg;
        else
            title(filenames.SC{i});
            fc = squeeze(fc_all(i,:,:));
            sc = squeeze(sc_all(i,:,:));
        end
        hold on;


        % plot #1: imagesc

        subplot(PLOTS,SUBJ_TO_SHOW,i);
        imagesc(fc);
        %axis([0 size(fc,1) 0 size(fc,2)]);
        if i==1
            ylabel('4 run mean fc');
        end

        % plot #2: FC
        subplot(PLOTS,SUBJ_TO_SHOW,SUBJ_TO_SHOW+i);
        hist([fc(:); LEFT; RIGHT],BINS, 'FaceColor', 'y');
        hold on;

        %% remove from both FC and SC all connections with SC == 0
        sc_nonzero = sc;
        sc_nonzero(sc == 0) = nan;

        fc_nonzero_sc = fc;
        fc_nonzero_sc(sc == 0) = nan;

        hist([fc_nonzero_sc(:); LEFT; RIGHT],BINS, 'FaceColor', 'g');

        plot([0 0],[0 HEIGHT],'k');
        if i==1
            ylabel(sprintf('fc w/wo. sc zeros\n(yellow vs green)'));  %  num2str(nansum(fc(:)))
        end
        axis([LEFT RIGHT 0 HEIGHT]);

        % plot #3: transformed SC and original FC
        subplot(PLOTS,SUBJ_TO_SHOW,SUBJ_TO_SHOW*2+i);

        transformed_sc = (sc_nonzero.^fcscfit.b)*fcscfit.a+fcscfit.c;

        % If IGNORE_ZERO_SC then turn cells that are not in the mask,
        % but still NaN into -999, otherwise turn them into 0
        if any(isnan(transformed_sc))
            ['mfsc: WARNING. transformed average SC has a SC cell that is zero and within the mask. Set to -999']
            transformed_sc(isnan(transformed_sc)) = -999;
        end

        connectome_write([OUTPUT_DIR filesep GROUP_FOLDER filesep 'transformed_SC_avg.csv'],transformed_sc);

        % if negatives FC are not taken into account in fitting, don't show
        % them
        if ~is_keep_neg_fc && i==1
            I = fc_nonzero_sc < 0;
            fc_nonzero_sc(I) = I(I) * NaN;
            transformed_sc(I) = I(I) * NaN;
        end 
        

        histogram([transformed_sc(transformed_sc < 1); LEFT; RIGHT],BINS,'FaceColor','m');
        hold on;
        histogram([fc_nonzero_sc(fc_nonzero_sc < 1); LEFT; RIGHT],BINS,'FaceColor','g');
        plot([0 0],[0 HEIGHT],'k');
        xlabel(['fc:' num2str(round(sum(fc_nonzero_sc(fc_nonzero_sc >= 0)))) ...  % only greater than zero, like van den Huevel 2017
                ' sc:' num2str(round(sum(transformed_sc(transformed_sc >= 0))))]);  % only greater than zero, so will be fair
        if i==1
            ylabel(sprintf('fc vs trans. sc\nboth wo. sc zeros\n(green vs magneta)'));
        end
        axis([LEFT RIGHT 0 HEIGHT]);

    end
 
end

if CALCULATE_SCFC  % formerly, this part was in "calculate_all_laterality_reg.m"

    %%% set these two if want to see specific connection
    %roi1 = 
    %roi2 =
    
    rs = [];
    intersects = [];
    slopes = [];
    mismatches = [];
    L_gt_Rs = [];
    FC_L_gt_Rs = [];
    order_i = 0;
    
    %figure;
    
    for i=1:numSubjs
        %for i= [4, 14, 19, 46]  %[*13, 14 , *15, 19] %35]   %19 - both have large mismatch %length(subjs)
        
        order_i = order_i + 1;
        
        % calculate SCFC for each subject
        
        [~,filename_SC,~] = fileparts(filenames.SC{i});
        
        is_plot_reg = false;
        is_plot_example = false;
        
        if is_figures % only do it if user requested figures to be printed out
            % Participants for Fig. 5
            if any(strcmp(filename_SC,'131217'))
                is_plot_example = true;
            % Participants for Fig. 2 
            elseif any(strcmp(filename_SC,{'130316','133928','149337','211720'})) %'130316',
                is_plot_reg = true;
            end
        end
        
        % OC - check what FC_PREFIX is for
        FC_PREFIX = 'FC';
        % calculate_mismatch_for_subject(sc,fc,a,b,c,PREFIX,subj_name,order_i)
        [r,intersect,slope,mismatch,L_gt_R,FC_L_gt_R] = calculate_mismatch_for_subject(sc_all(i,:,:),fc_all(i,:,:),fcscfit.a,fcscfit.b,fcscfit.c,FC_PREFIX,filenames.SC{i},order_i,combined_mask_sym,is_plot_reg,is_plot_example);

        rs = [rs r];
        intersects = [intersects intersect];
        slopes = [slopes slope];
        mismatches(i,:,:,:) = mismatch;
        L_gt_Rs(i,:,:,:) = L_gt_R;
        FC_L_gt_Rs(i,:,:,:) = FC_L_gt_R;

        % turn all cells that are not in the mask to ignore, including sc=0
        % cells that might have been set to 0
        mismatch(~combined_mask_sym) = exclude_value; 
        % If IGNORE_ZERO_SC then turn cells that are not in the mask,
        % but still NaN into -999, otherwise turn them into 0
        if any(isnan(mismatch))
            [~,filename_SC,~] = fileparts(filenames.SC{i});
            ['mfsc: WARNING. ' filename_SC ' has a SC cell that is zero and within the mask. Set to -999']
            mismatch(isnan(mismatch)) = -999;
        end
        
        [~,filename_FC,~]=fileparts(filenames.FC{i});
        [~,filename_SC,~]=fileparts(filenames.SC{i});

        % mirror the matrix
        connectome_write([OUTPUT_DIR filesep 'mFCSC-' filename_FC '-' filename_SC '-masked.csv'],mismatch);
        
        clear r intersect slope mismatch L_gt_R;
        order_i = order_i + 1;
        
    end

end

fprintf('MFCSC finished successfuly\n\n');

if isdeployed
    exit 
end
