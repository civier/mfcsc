function direct_is_shortest_mask = calculate_direct_is_shortest_mask_sym(sc_avg,a,b,c,is_figures)
% OC - do not run Roshenthal anymore!

% Return mask in original order of data, NOT rosenthal!

% Original is 
% /Users/ocivier/Dropbox/Data_Analyzed/Connectome/SCFC_4_runs-WITH_ORIGINAL_FC/direct_is_shortest_mask

dbstop if error;

OPTION = 4;
SAVE_MASK = false;

% if OPTION == 4
%     % OPTION 4 - avg of all normalized connectomes
%     [normalized_sc_avg,~] = whole_brain_normalized_average('ROSE');
%     W = prepare_for_bct(normalized_sc_avg);    
%     W = W - min(min(W));     % for the avg of the 50 subjects, it's -0.2308
% elseif OPTION == 3
%     % OPTION 3 - normalization I did for abstract (but shift above 0)
%     load('ROSE_189450.mat');
%     W = prepare_for_bct(reshape(normalized_sc,84,84));
%     W = W - min(min(W)); % shift values so all weights are above 0. Usually around -0.2X
%     %W = W .* (W>0);  % eliminate negatives
% elseif OPTION == 1
%     % OPTION 1 - reading original connectome
%     %W = prepare_for_bct(dlmread('129028.csv'));
% elseif OPTION == 2
%     % OPTION 2 - coarse normalization (^0.25)
%     %W = W.^0.25; % normalize
% end

%load('SC_avg_50_subjects.mat');
%load('FC_avg_50_subjects.mat');
%fc_avg_MRtrix = fc_avg;

%[sc_avg_rosenthal,~] = prepare_sc_fc_Rosenthal(sc_avg,sc_avg);

  %     a =      0.4114;
  %     b =     0.09267;
  %     c =     -0.3789;

% calculate normalize matrix
W = (sc_avg.^b).*a + c;

% OC - not required
% [Orig_transformed_sc,fc_avg] = prepare_sc_fc_Rosenthal(W,fc_avg_MRtrix);

% zero out negatives
W = W .* (W>=0);

% make symmetrical for BCT
W = mirror_mat(W);

L_MRtrix = W;                                                      % connection-length matrix
A = W > 0;                                                  % adjacency matrix
L_MRtrix(A) = 1 ./ L_MRtrix(A);

ind_L_MRtrix = distance_wei(L_MRtrix);

% transfer to rosenthal for visualization and calculation of two
% hemispheres!
% OC - not necessary to transfer it to Rosenthal order
%[ind_L,L] = prepare_sc_fc_Rosenthal(ind_L_MRtrix,L_MRtrix);
% OC - instead, just copy to new names
ind_L = ind_L_MRtrix;
L = L_MRtrix;

%M = charpath(distance_wei(L),0,0);

% Find strongest connection = 23,21
[I, val] = find(W == max(max(W)));

if is_figures % plot only if is_figures is on

    figure;
    
    subplot(231);
    imagesc(W); colorbar; title('Original matrix');
    subplot(232);
    
    %imagesc(L); colorbar; title('distance of direct paths'); - TOTALLY DIFF
    %SCALE THAN INDIRECT DISTANCE MATRIX

    mat_to_plot = L .* (L<max(max(ind_L(~isinf(ind_L)))));
    mat_to_plot(1,1) = 0; % add zero in case lowest value is not 0
    imagesc(mat_to_plot); colorbar; title('Lengths of direct paths. Large values are zeroed out'); % drop long-distance connections to make scale identical to indirect distance matrix

    subplot(233);
    mat_to_plot = ind_L;
    mat_to_plot(1,1) = 0; % add zero in case lowest value is not 0
    imagesc(mat_to_plot); colorbar; title('Lengths of shortest paths, either direct or not');

    subplot(236);
    imagesc(ind_L == L); title(['Length of direct == Length of shortest.' num2str(sum(sum(ind_L == L))/2-length(W)) ' cases out of ' num2str((length(W)^2)/2)]); colorbar; % remove length(W) beacuse diagonal does not count

end

%%% ONLY 3% OF SHORTEST PATHS ARE DIRECT IN REGULAR CONNECTOM
%%% WITH NORMALIZATION, IT IS ABOUT 25% (the weakest is 0.3656, on scale of 0 to 2)
%%%                                     (the ones where direct ~ indirect
%%%                                      can go up to 0.8858 on that scale)

%%% FOR AVERAGE SC:
%%%       weakest direct == indirect is 0.3283
%%%    Observations: 1) More likely for ipsilateral connections
%%%                  2) Homotopic nodes are more likely to have
%%%                     direct==indirect than non-homotopic
                     
N = size(ind_L,1);

ind_eq_dir = (ind_L == L);

% OC - zero the diagonal, because should not be incldued in any calculation
% anyway
ind_eq_dir(logical(eye(length(ind_eq_dir))))=0;

% NOTICE THAT IT IS IN THE ORIGINAL ORDERING
direct_is_shortest_mask = ind_eq_dir;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%  BELOW IS EXTRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







figure; 
plot(Orig_transformed_sc,fc_avg,'.b');
hold on;
plot(Orig_transformed_sc(ind_eq_dir),fc_avg(ind_eq_dir),'.r');
title(['all: ' num2str(corr(Orig_transformed_sc(~isnan(fc_avg) & ~isnan(Orig_transformed_sc)),fc_avg(~isnan(fc_avg) & ~isnan(Orig_transformed_sc)))) ... 
       '; only direct=shortest (red): ' num2str(corr(Orig_transformed_sc(ind_eq_dir),fc_avg(ind_eq_dir)))]);
xlabel('transformed SC');
ylabel('FC');
figure;
imagesc(ind_eq_dir);
title('direct=shortest mask');

figure;
ind_eq_dir_ipsi = ind_eq_dir;
ind_eq_dir_ipsi(1:N/2,(N/2+1):N) = 0;
plot(Orig_transformed_sc(ind_eq_dir),fc_avg(ind_eq_dir),'.r');
hold on;
plot(Orig_transformed_sc(ind_eq_dir_ipsi),fc_avg(ind_eq_dir_ipsi),'.g');
xlabel('transformed SC');
ylabel('FC');
title(['all: ' num2str(corr(Orig_transformed_sc(ind_eq_dir),fc_avg(ind_eq_dir))) ... 
       '; only ipsilateral (green): ' num2str(corr(Orig_transformed_sc(ind_eq_dir_ipsi),fc_avg(ind_eq_dir_ipsi)))]);
save('inq_eq_dir_ipsi_mask.mat','ind_eq_dir_ipsi');

figure;
diff = (ind_eq_dir_ipsi(1:42,1:42) ~= ind_eq_dir_ipsi((1:42)+42,(1:42)+42)); % calculate connections where direct=shortest only on one side
imagesc(diff);
title(['diff between hemis, total ' num2str(sum(sum(diff))) ' out of ~' num2str(sum(sum(ind_eq_dir_ipsi))/2)]);

figure;
sum(sum(ind_eq_dir)) / 2
% 936 out of 3486, about 25%
ind_eq_dir_in_both_hemi = ind_eq_dir(1:N/2,1:N/2) .* ind_eq_dir(N/2+(1:N/2),N/2+(1:N/2));
'direct is shortest in both hemi, out of 861'
sum(sum(ind_eq_dir_in_both_hemi)) / 2
% 304 of 861, about 30% (so quite consistent, but also because that no
% inter-hemisphere connections)

subplot(235);
imagesc(ind_eq_dir_in_both_hemi); colorbar;
title(['direct in shortest path in both hemispheres.' num2str(sum(sum(ind_eq_dir_in_both_hemi)) / 2) ' our of 861']);

if SAVE_MASK
    dlmwrite('4_RUNS_EX_ABS_MASKED_DIRECT_SHORTEST_MASK.txt',ind_eq_dir_in_both_hemi);
end

% make direct_is_shortest 84x84
ind_eq_dir_in_both_hemi_big(1:42,1:42) = ind_eq_dir_in_both_hemi;
ind_eq_dir_in_both_hemi_big(42+(1:42),42+(1:42)) = ind_eq_dir_in_both_hemi;

MRtrix_ind_eq_dir_in_both_hemi_big = prepare_sc_fc_MRtrix(ind_eq_dir_in_both_hemi_big,ind_eq_dir_in_both_hemi_big);

MRtrix_RgtL = dlmread('BIN_MRtrix_ROSE_MASKED_none_output_ISE_0uncorrected_pvalue_t1_BEN_BIN_SCFC.csv');
MRtrix_LgtR = dlmread('BIN_MRtrix_ROSE_MASKED_none_output_ISE_0uncorrected_pvalue_t2_BEN_BIN_SCFC.csv');


% USE DIRECT IS SHORTEST TO MASK

figure;
subplot(231);
imagesc(MRtrix_LgtR);
title(['L gt R' num2str(sum(MRtrix_LgtR(:))/2)]);
subplot(232);
imagesc(MRtrix_RgtL);
title(['R gt L' num2str(sum(MRtrix_RgtL(:))/2)]);
subplot(233);
imagesc(MRtrix_ind_eq_dir_in_both_hemi_big);
title('direct is shortest in both hemis');

MRtrix_LgtR_MASKED = MRtrix_LgtR .* MRtrix_ind_eq_dir_in_both_hemi_big;
MRtrix_RgtL_MASKED = MRtrix_RgtL .* MRtrix_ind_eq_dir_in_both_hemi_big;

subplot(234);
imagesc(MRtrix_LgtR_MASKED);
title(['L gt R - MASKED' num2str(sum(MRtrix_LgtR_MASKED(:))/2)]);
subplot(235);
imagesc(MRtrix_RgtL_MASKED);
title(['R gt L - MASKED' num2str(sum(MRtrix_RgtL_MASKED(:))/2)]);
c