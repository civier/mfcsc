function [return_r,return_intersect,return_slope, mismatch, left_gt_right, FC_left_gt_right] = calculate_mismatch_for_subject(sc,fc,a,b,c,PREFIX,subj_name,order_i, mask, IS_PLOT_REG,IS_PLOT_EXAMPLE)
% OC - do not run Roshenthal anymore, so all laterality calculations are
% wrong!

global fig_2_order;

% formerly was named SUBJ_1_plot_connectome_90_mismatch_FS_reg

dbstop if error;

sc = squeeze(sc);
fc = squeeze(fc);

%a =      0.4114;
%b =     0.09267;
%c =     -0.3789;

DO_EXTRA = false;
IS_PLOT = false;
SCALING = true;
FISHER = false;
SAVE_FILE_FOR_STATS = false;
%PREFIX = 'ALL'; %'LRFCgtRLSC' (if only wants to consider when LRFC > RLSC)

THREE_SUBJ = false;

sc_THRESH = 10; %10000; %70000/100;
NORM_FUNC = 'norm_dist'; %'log';

%sc = dlmread([SC_DIR '/' SC_PREFIX subjs{i} '.csv']);
%fc = dlmread([FC_DIR '/' FC_PREFIX subjs{i} '.csv']);

%%%% JUST SO CAN VERIFY LABELS ARE CORRECT %%% 
i = 1; % 1 for FS
BASE = ['.' filesep];
label_files = {[BASE 'MRtrix_FS_FIRST_just_labels.txt'],[BASE 'AAL_just_labels.txt']};


% OC - instead of the actual lables (commented out), just use serial
% numbers

%fid = fopen(label_files{i});
%labels = textscan(fid,'%s','delimiter','\n');
%labels = labels{1};
%fclose(fid);
labels = 1:length(sc);

%%%%%%%%%%%%%%%%%%%%%%%%%%

% OC - do not run Roshenthal anymore, so all laterality calculations are
% wrong!
%[sc,fc,labels] = prepare_sc_fc_Rosenthal(sc,fc,labels);

if size(fc) ~= size(sc)
    error('index matrices have different sizes')
end
N = length(fc);

% turn diagonal to NaN
for i=1:N
    fc(i,i) = nan;
    sc(i,i) = nan;
end

% nan lower part of functional and structual matrices
% and transfer value to other half if necessary
% This is important because the switch from MRtrix ordering to Rosenthal
% ordering might move some values to the lower triangle
for i=1:N
    for j=i:N
        if fc(j,i) > 0
            fc(i,j) = fc(j,i);
        end
        fc(j,i) = nan;

        if sc(j,i) > 0
            sc(i,j) = sc(j,i);
        end
        sc(j,i) = nan;
    end
end

% turn them into vectors for later processing
fc_vec = fc(:);
sc_vec = sc(:);

% Used to check correct number of NaNs -- but not practical if some columns
% are empty
%if sum(isnan(fc_vec)) ~= N*N/2+N/2
%    error('bad NANs in fc');
%end
%if sum(isnan(sc_vec)) ~= N*N/2+N/2
%    error('bad NANs in sc');
%end

% turn connections with 0 connectivity to NaNs
% It doesn't matter, because in any case the direct is nont the shortest
% (direct is 0)
% if using 2-power law, need to remove to get rid from left edge
zero_idx = (sc_vec == 0);
% turn structural with 0 SC to nan
sc_vec(zero_idx) = nan;
% turn corresponding functional to nan
fc_vec(zero_idx) = nan;

% Used to check the same number of NaNs
if sum(isnan(fc_vec)) ~= sum(isnan(sc_vec))
    error('bad NANs in sc or fc');
end

%%%%% OUTLIERS %%%%%%
if SCALING
    
    normalized_fc = fc_vec;

    %%% MAKE SURE THAT LOWER PART OF MATRIX IS NAN, OTHERWISE WILL GET SOME
    %%% VALUES
    normalized_sc = (sc_vec.^b).*a+c;
        
    sc_stronger_than_fc_vec = normalized_sc - normalized_fc;
else
    sc_stronger_than_fc_vec = nanrank(sc_vec) - nanrank(fc_vec);
end
sc_stronger_than_fc = reshape(sc_stronger_than_fc_vec,N,N);

% calculate the linear regression line of the subject
    %load ind_eq_dir_ipsi_mask.mat
    %ind_eq_dir_ipsi = ones(size(sc));
    
    
    normalized_sc_mat = reshape(normalized_sc,N,N);
    % excluding all connections that are not in mask
    normalized_sc_mat(~mask) = NaN;
    
    normalized_fc_mat = reshape(normalized_fc,N,N);
    % excluding all connections that are not in mask
	normalized_fc_mat(~mask) = NaN;
    
    p = polyfit(normalized_sc_mat(~isnan(normalized_sc_mat)),normalized_fc_mat(~isnan(normalized_fc_mat)),1);
    return_r = corr(normalized_sc_mat(~isnan(normalized_sc_mat)),normalized_fc_mat(~isnan(normalized_fc_mat)));
    return_intersect = p(2);
    return_slope =  p(1);

% calculate mismatch matrix for the subject
    mismatch = normalized_fc_mat - (return_slope .* normalized_sc_mat + return_intersect); % calculate mismatch

% calculate (mismatch_left - mismatch_right) for the subject and save it to a .txt file    
    left_gt_right = mismatch(1:N/2,1:N/2) - mismatch(N/2+(1:N/2),N/2+(1:N/2)); % automatically takes care of NaN, because operation between two NaNs is NaN    
    % BETTER NOT SAVE, BECAUSE ROSENTHAL ORDER - dlmwrite([PREFIX '_' subj_name '.txt'],left_gt_right); %mirror_mat((-a) - (-b)));    

% calculate (FC_left - FC_right) and save it to a .txt file    
    FC_left_gt_right = normalized_fc_mat(1:N/2,1:N/2) - normalized_fc_mat(N/2+(1:N/2),N/2+(1:N/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FROM HERE BELOW IT IS OPTIONAL!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot linear regression per subject
if IS_PLOT_REG || IS_PLOT_EXAMPLE

%    if IS_PLOT_REG
%        if isempty(fig_2_order)
%               fig_2_order = 1;
%        end
%
%        subplot(2,2,fig_2_order);
%        fig_2_order = fig_2_order + 1;
%    end
    figure;
    
    grayColor = [.7 .7 .7];
    plot(normalized_sc_mat(:),normalized_fc_mat(:),'d','MarkerSize',3,'Color',grayColor,'MarkerFaceColor', grayColor);
    set(gca,'linewidth',1);
    hold on;
    plot([0 1],return_slope*[0 1]+return_intersect,'-k','LineWidth', 1);
    xlabel('SC^{norm}');
    ylabel('FC');
    title(['subject ' subj_name]);
    
    % make font of tics texr larger
    ax=gca;
    ax.FontSize = 14;
    % this is required because the larger font causes for matlab to reduce
    % the number of ticks
    xticks(0:0.1:1)

    %title(['r = ' num2str(return_r) ' intersect = ' num2str(return_intersect) ' slope = ' num2str(return_slope)]);

    if IS_PLOT_EXAMPLE
        % highlight mismatches of a specific connection
        %       L1 L2  R1  R2
        %figure;
        rois= [ 7  17 49+7 49+17
                26 28 49+26 49+28
                9  36 36+7  49+9
                1  8  49+1  49+8
                3  10 49+3  49+10];
        % OLD - MARKER_SHAPE = ['o','s','d','p','h'];
        
        MARKER_SHAPE = ['s','h','o','p','d'];
%Now circle SHOULD BE square
%Now square SHOULD BE hexagram
%Now diamond SHOULD BE circles
%Now hexagram SHOULD BE diamond
    elseif IS_PLOT_REG
        %figure;
        rois= [ 17  27 49+17 49+27];
        MARKER_SHAPE = ['d'];
    end     

    for roi_i=1:size(rois,1)
    
        for side=0:1

            roi1=rois(roi_i,side*2+1);
            roi2=rois(roi_i,side*2+2);

            % notice that smaller region in order need to be first, as it is an
            % upper triangular

            % THE REGIONS USED IN FIG 2

            %%%%%%%%%%%%%%%%%%%%%
            %%% TESTING_FIG2 - SET ROI1 and ROI2 below to draw a specific connection
            %%%                on all participants;
            %%%                * Notice that ROI1 shold be smaller than
            %%%                ROI2 both when side=0 and when side=1!
            %%%                * The lables are MRtrix3 Desikan labels.
            %%%                * Notice that it won't be exactly like the article
            %%%                because the linear regression was done on all
            %%%                connections where direct SC is shortest path,
            %%%                even if it's only in one hemisphere.
            %%%                To replicate the article 100%, should use
            %%%                /Users/ocivier/Dropbox/Data_Analyzed/Connectome/test-SCFC_PNAS/SUBJ_1_plot_connectome_90_mismatch_FS_reg.m 
            %%%%%%%%%%%%%%%%%%%%%

            % The connection showed in Fig. 2
            %roi1 = side*49 + 17;  % parsopercularis
            %roi2 = side*49 + 27;  % superiorfrontal

            % Square
            %roi1 = side*49 + 7;   % Inferiorparietal
            %roi2 = side*49 + 17;  % parsopercularis

            % Six point star
            %roi1 = side*49 + 26;  % Rostralmiddlefrontal
            %roi2 = side*49 + 28;  % Superiorparietal

            % Circle
            %if side == 0
            %    roi1 = side*49 + 9;  % isthmuscingulate
            %    roi2 = 36 + side*7;  % Thalamus proper (doesn't work, because NaN)
            %else % need to shift between roi1 and roi2 becuase smaller should be first
            %    roi2 = side*49 + 9;  % isthmuscingulate
            %    roi1 = 36 + side*7;  % Thalamus proper (doesn't work, because NaN)
            %end

            % Five point star
            %roi1 = side*49 + 1;  % inferior temporal
            %roi2 = side*49 + 8;  % bankssts

            % Diamond
            %roi1 = side*49 + 3;  % caudal middle frontal
            %roi2 = side*49 + 10;  % lateral occipital

            
            if side == 0
                color = 'r';
                hemi = 'L';
            else
                color = 'b';
                hemi = 'R';
            end
            % OC - remove the below if want to draw in colos
            color = 'k';
            
            % plot line first, so shape will overlay the line
            plot([normalized_sc_mat(roi1,roi2) normalized_sc_mat(roi1,roi2)], ...
                 [normalized_fc_mat(roi1,roi2) (return_slope .* normalized_sc_mat(roi1,roi2) + return_intersect)],color,'LineWidth',2);

            % plot a shape filled in white. A bit larger if plotting the
            % example
            plot(normalized_sc_mat(roi1,roi2),normalized_fc_mat(roi1,roi2),[MARKER_SHAPE(roi_i) color],'MarkerFaceColor', 'w', 'MarkerSize', 12+IS_PLOT_EXAMPLE*2, 'LineWidth', 1);
            
            % plot either L or R above diamond. Make it closer and a bit larget if plotting
            % the multi pairs example
            text(normalized_sc_mat(roi1,roi2),normalized_fc_mat(roi1,roi2)+0.05 - IS_PLOT_EXAMPLE*0.03,hemi,'HorizontalAlignment','center','FontSize',10+IS_PLOT_EXAMPLE,'FontWeight','bold');
            
            axis ([0 1 0 1]);
            axis square;
        end

    end
        %%%

        %figure; imagesc(mismatch); title('mismatch');
        %figure; plot(mismatch(1:N/2,1:N/2),mismatch(N/2+(1:N/2),N/2+(1:N/2)),'*');
        %hold on; plot([-0.5 0.5],[-0.5 0.5],'k');
        %axis([-1 1 -1 1]);    
        %figure; imagesc(left_gt_right); title('left_gt_right');        
end

if DO_EXTRA


    % plot matrices
    if IS_PLOT
        figure;
        hold on
        imagesc(fc)
        plot([N/2,N/2]+0.5,[0,N]+0.5,'k-'); plot([0,N]+0.5,[N/2,N/2]+0.5,'k-');

        axis xy;
        title('functional')
        colorbar;
        axis tight

        figure;
        imagesc(log(sc))
        hold on
        plot([N/2,N/2]+0.5,[0,N]+0.5,'k-'); plot([0,N]+0.5,[N/2,N/2]+0.5,'k-');
        axis xy
        colorbar;
        colormap(gray)
        title('scructural (log)')
    end

    % plot the differnce between functional and structural ranking
    if IS_PLOT
        figure;
        imagesc_equal_colorbar(-sc_stronger_than_fc) % minus to have more FC in red!
        hold on
        plot([N/2,N/2]+0.5,[0,N]+0.5,'k-'); plot([0,N]+0.5,[N/2,N/2]+0.5,'k-');

        colorbar;
        axis xy;
        axis square
        title('diff between sc and fc;   red = fc>sc;  blue = sc>fc');
    end

    % plot least difference in ranking between hemispheres
    if IS_PLOT
        figure;
        imagesc(N*N/2 - abs(sc_stronger_than_fc));
        colormap(gray)
        hold on
        plot([N/2,N/2]+0.5,[0,N]+0.5,'k-'); plot([0,N]+0.5,[N/2,N/2]+0.5,'k-');

        colorbar;
        axis xy;
        axis square
        title('match between sc and fc;   white = high;   black = low');
    end

    % plot cases where the mismatch in structural is lager
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if IS_PLOT
        figure; subplot(1,2,1)
    end
    if SCALING
        sc_rank_vec = normalized_sc;
        fc_rank_vec = normalized_fc;
    else
        sc_rank_vec = nanrank(sc_vec);
        fc_rank_vec = nanrank(fc_vec);
    end

    sc_rank = reshape(sc_rank_vec,N,N);
    fc_rank = reshape(fc_rank_vec,N,N);

    % check if mismatch in structural is greater than functional
    sc_mis = abs(sc_rank(1:N/2,1:N/2)-sc_rank(N/2+(1:N/2),N/2+(1:N/2))) > ...
             abs(fc_rank(1:N/2,1:N/2)-fc_rank(N/2+(1:N/2),N/2+(1:N/2)));

    % plot cases where structural is greater
    mat_sc = (sc_stronger_than_fc(1:N/2,1:N/2)-sc_stronger_than_fc(N/2+(1:N/2),N/2+(1:N/2))) .* sc_mis;
    %mat_sc(N/2,1) = -max(mat_sc(:)); % to make sure that colorbar will be around zero
    %mat_sc(N/2-1,1) = -min(mat_sc(:)); % "
    if IS_PLOT
        imagesc_equal_colorbar(mat_sc);
    end
    %colorbar;
    if IS_PLOT
        title(sprintf('different sc, similar fc\n red = more sc in left,  blue = more sc in right'));
        axis xy;
        axis square

    % plot cases where functional is greater

        subplot(1,2,2)
    end
    mat_fc = (-sc_stronger_than_fc(1:N/2,1:N/2)+sc_stronger_than_fc(N/2+(1:N/2),N/2+(1:N/2))) .*~sc_mis;
    %mat_fc(N/2,1) = -max(mat_fc(:));  % to make sure that colorbar will be around zero
    %mat_fc(N/2-1,1) = -min(mat_fc(:));  % "
    if IS_PLOT
        imagesc_equal_colorbar(mat_fc);
    %colorbar;
        title(sprintf('different fc, similar sc\n red = more fc in left,  blue = more fc in right'));


        axis xy;
        axis square
    end

    % plot fc right vs. left and sc right vs. left (disregard if the other is
    % similar or not)

    if IS_PLOT
        figure;
        subplot(1,2,1);
    end
    mat_sc_abs = sc_rank(1:N/2,1:N/2)-sc_rank(N/2+1:N,N/2+1:N)
    %mat_sc_abs(N/2,1) = -max(mat_sc_abs(:));  % to make sure that colorbar will be around zero
    %mat_sc_abs(N/2-1,1) = -min(mat_sc_abs(:));  % "

    %colorbar;
    if IS_PLOT
        imagesc_equal_colorbar(mat_sc_abs);
        title(sprintf('sc left vs. right\n red = more sc in left,  blue = more sc in right'));
        axis xy;
        axis square
    end

    % plot cases where functional is greater
    if IS_PLOT
        subplot(1,2,2);
    end
    mat_fc_abs = fc_rank(1:N/2,1:N/2)-fc_rank(N/2+1:N,N/2+1:N);
    %mat_fc_abs(N/2,1) = -max(mat_fc_abs(:));  % to make sure that colorbar will be around zero
    %mat_fc_abs(N/2-1,1) = -min(mat_fc_abs(:));  % "

    if IS_PLOT
        imagesc_equal_colorbar(mat_fc_abs);
        %colorbar;
        title(sprintf('fc left vs. right\n  red = more fc in left, blue = more fc in right'));

        axis xy;    
        axis square
    end



    % scatter plot comparing left and right
    if IS_PLOT
        figure;
    end
    a = -sc_stronger_than_fc(1:N/2,1:N/2);
    b = -sc_stronger_than_fc(N/2+(1:N/2),N/2+(1:N/2));
    if IS_PLOT
        plot(a(mat_fc_abs > 0.4),b(mat_fc_abs > 0.4),'.r');
        hold on;
        plot(a(mat_fc_abs < 0.4),b(mat_fc_abs < 0.4),'.b');
        xlabel('sc > fc       LEFT HEMI       fc > sc');
        ylabel('sc > fc       RIGHT HEMI      fc > sc');
        hold on;
        plot([-1,1],[0,0],'k'); plot([0,0],[-1,1],'k');
        plot([-1,1],[-1,1],'k');
        title([subj_name ' red - fc > 0.4, blue - fc < 0.4']);
    end

    % do below if only consider when FC L > FC R and FC L - FC R > SC R - SC L
    if strcmp(PREFIX,'LRFCgtRLSC')
        mat_norm_fc = reshape(normalized_fc,N,N)
        mat_norm_sc = reshape(normalized_sc,N,N);
        for  i=1:size(a,1)
            for j=1:size(a,2)
                if mat_norm_fc(i,j) < mat_norm_fc(i+45,j+45) || ...
                    mat_norm_fc(i,j) - mat_norm_fc(i+45,j+45) < mat_norm_sc(i+45,j+45) - mat_norm_sc(i,j)
                    a(i,j) = 0;
                    b(i,j) = 0;
                    zerod(i,j) = true; % about half (546) of the ~1000 (1012 - 22) connections in the triangle are zerod
                end
            end
        end
    end

    if SAVE_FILE_FOR_STATS
        dlmwrite([PREFIX '_' subj_name '.txt'],mirror_mat((-a) - (-b)));
    end
    save([PREFIX '_' subj_name],'mat_sc_abs','mat_fc_abs','normalized_fc','normalized_sc');

    if IS_PLOT
        plot_sc_vs_fc_laterality_FS
    end

    end



    % 
    % 
    % 
    % figure;
    % mat_fc_comp = mat_fc_abs - mat_sc_abs;
    % %mat_fc_comp(N/2,1) = -max(mat_fc_comp(:));  % to make sure that colorbar will be around zero
    % %mat_fc_comp(N/2-1,1) = -min(mat_fc_comp(:));  % "
    % imagesc_equal_colorbar(mat_fc_comp);
    % 
    % %colorbar;
    % 
    % title(sprintf('(fc left vs. right) minus (sc left vs. right)\n  red = more fc in left, blue = more fc in right'));
    % 
    % axis xy;
    % axis square
    % 
    % if SCALING
    %     BOUND = 0.2;
    % else
    %     BOUND = 400;
    % end
    % 
    % % PLOT BELOW BOUND
    % figure; subplot(121)
    % mat_fc_BOUND = mat_fc
    % mat_fc_BOUND(mat_fc_BOUND > BOUND | mat_fc_BOUND < -BOUND) = 0
    % imagesc(mat_fc_BOUND)
    % title('fc relative up to BOUND')
    % axis xy
    % colorbar
    % 
    % subplot(122)
    % mat_fc_abs_BOUND = mat_fc_abs
    % mat_fc_abs_BOUND(mat_fc_abs_BOUND > BOUND | mat_fc_abs_BOUND < -BOUND) = 0
    % imagesc(mat_fc_abs_BOUND)
    % title('fc abs up to BOUND')
    % axis xy
    % colorbar
    % 
    % % PLOT ABOVE BOUND
    % figure; subplot(121)
    % mat_fc_BOUND = mat_fc
    % mat_fc_BOUND(mat_fc_BOUND < BOUND & mat_fc_BOUND > -BOUND) = 0
    % imagesc(mat_fc_BOUND)
    % title('fc relative above BOUND')
    % axis xy
    % colorbar
    % 
    % subplot(122)
    % mat_fc_abs_BOUND = mat_fc_abs
    % mat_fc_abs_BOUND(mat_fc_abs_BOUND < BOUND & mat_fc_abs_BOUND > -BOUND) = 0
    % imagesc(mat_fc_abs_BOUND)
    % title('fc abs above BOUND')
    % axis xy
    % colorbar
    % 
    % %subplot(133)
    % %imagesc(mat_sc_abs < -BOUND)
    % %title('sc the other direction gt BOUND')
    % %axis xy
    % %colorbar
    % 
    % %%% SHOW ALL OCCURENCES OF OPPOSITE DIRECTION ABOVE SOME EXTENT
    % figure;
    % thresh_opp = 300;
    % opp = zeros(45,45);
    % opp((mat_fc_abs > thresh_opp) & (mat_sc_abs < -thresh_opp)) = 1;
    % opp((mat_fc_abs < -thresh_opp) & (mat_sc_abs > thresh_opp)) = -1;
    % imagesc(opp)
    % title(['opposite direction, both above ' num2str(thresh_opp) '; red - fc left sc right; blue - fc right sc left'])
    % axis xy
    % 
    % %%% DRAWING CHANGE IN DISTRIBUTIONS FROM ABSOLUTE TO RELATIVE %%%
    % figure; colorhist(mat_sc_abs(:),100,'r');
    % hold on; colorhist(mat_sc(:),100,'g'); title('sc - from abs to relative')
    % 
    % figure; colorhist(mat_fc_abs(:),100,'r'); title('sc')
    % hold on; colorhist(mat_fc(:),100,'g'); title('fc - from abs to relative')
    % 
    % 
    % % find regions with comparable big change both in sc and fc
    % (abs(mat_sc) < 100) & (abs(mat_fc) < 100) & abs(sc(1:45,1:45) - sc(46:90,46:90)) > 10000
    % % for example, 6 to 35 have comparable change of almost 1000 places, with
    % % advantage to left!
    % 
    % %diausgdiuhsaidh
    % 
    % % remove NaNs to enable correlation
    % sc_vec_nnan = sc_vec(~isnan(sc_vec));
    % fc_vec_nnan = fc_vec(~isnan(fc_vec));
    % 
    % % normalize
    % figure;
    % subplot(1,2,1);
    % hold on;
    % hist(sc_vec_nnan,10);
    % sc_vec_nnan = eval([NORM_FUNC '(sc_vec_nnan);']);
    % subplot(1,2,2);
    % hist(sc_vec_nnan,10);
    % 
    % corr_fc_sc = corr(fc_vec_nnan,sc_vec_nnan)
    % Edges = length(fc_vec_nnan)
    % figure; plot(sc_vec_nnan,fc_vec_nnan,'*')
    % title(['half of ' num2str(N) 'x' num2str(N) ' connections without zero SC and ' NORM_FUNC ': r=' num2str(corr_fc_sc) ' n=' num2str(Edges)])
    % 
    % THRESH_IDX = find(sc_vec_nnan >= eval([NORM_FUNC '(sc_THRESH)']));
    % corr_fc_sc = corr(fc_vec_nnan(THRESH_IDX),sc_vec_nnan(THRESH_IDX))
    % Edges_Thresh = length(fc_vec_nnan(THRESH_IDX))
    % figure; plot(sc_vec_nnan(THRESH_IDX),fc_vec_nnan(THRESH_IDX),'*')
    % title(['half of only connections stronger than ' num2str(sc_THRESH) ': r=' num2str(corr_fc_sc) ' n=' num2str(Edges_Thresh)])
    % 
    % end
    
end