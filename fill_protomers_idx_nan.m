%%%%%%
%%% Identify the outlier POs (fill the nan values) cluster IDs after intial
%%% PCA assessment based 1) POs relative distance to the cluster center
%%% (measurement by z-score, and 2) POs positions in the protomer time
%%% series (not used by default).
%%% Input: 
%%%        1: cluster_idx_col: % column in protomer_idx file to be
%%%           addressed (recording intial PCA assessment)
%%%        2: protomers_idx: recording PO infos (tipj IDs and cluster IDs)
%%%        3: score: score values for calculation POs relative distances to
%%%           cluster centers
%%%        4: exclu_sel: excluded POs 
%%%        5: out_sel: outlier POs 
%%%        6: idx0: baseline cluster ID. 
%%%           Default: 0 for OFS-IFS clustering. 1 for IFSo-IFSc cluster.
%%%        7. min_mean_zscore_diff: minimal zscore difference. Default: 0 
%%%           (note that choice of 0 will accept all IDs assigned at this 
%%%           stage, therefore no time-serie comparisons will be made)
%%%        8. Oligomer molecular symmetry value
%%% Output:
%%%        1: protomers_idx: updated protomer_idx file
%%%%%%
function protomers_idx = fill_protomers_idx_nan(cluster_idx_col, protomers_idx, score, exclu_sel, out_sel, idx0, min_mean_zscore_diff, nf)

%%% fill the gap 1 - compare zscore
%%% purpose:
% identify the outlier POs according to their relative distance to the
% cluster center (measured by z-score)

%%%
[~, Ncol] = size(score);
protomers_idx_out = zeros(sum(out_sel), 1);
protomers_up_stats = score(~exclu_sel, :);
protomers_out_stats = protomers_up_stats(out_sel, :);
num_states = max(protomers_idx(:, cluster_idx_col));

% record the mean of std of the POs in all clusters (reference points for
% identifying outlier POs
protomers_ref_stats = zeros(num_states, Ncol, 2);
for i = 1:num_states
    sel = protomers_idx(:, cluster_idx_col) == i + idx0;    
    protomers_sel_stats_i = score(sel, :);

    % calculate the mean and std of the POs in the two clusters
    protomers_ref_stats(i, :, 1) = mean(protomers_sel_stats_i(:, :), 1);
    protomers_ref_stats(i, :, 2) = std(protomers_sel_stats_i(:, :), 1);
end

% calculate z-scores of all outlier POs (distances to the cluster centers)
protomers_out_zscores = zeros(sum(out_sel), num_states, Ncol);
for p = 1:sum(out_sel)
    protomers_out_stats_p = protomers_out_stats(p, :);
    for i = 1:num_states
        zscores_p_i = (protomers_out_stats_p-protomers_ref_stats(i, :, 1))./protomers_ref_stats(i, :, 2);
        protomers_out_zscores(p, i, :) = zscores_p_i;
    end
end
protomers_out_zscore_mean = sqrt(sum(protomers_out_zscores.^2, 3));
protomers_out_zscore_mean = protomers_out_zscore_mean./(Ncol-1);

% find the closest culster
[protomers_out_zscore_min1, protomers_out_zscore_idx] = min(protomers_out_zscore_mean, [], 2);
protomers_out_zscore_mean(protomers_out_zscore_mean==protomers_out_zscore_min1) = Inf;
[protomers_out_zscore_min2, ~] = min(protomers_out_zscore_mean, [], 2);

protomers_out_zscore_diff = protomers_out_zscore_min2-protomers_out_zscore_min1;
sel1 = protomers_out_zscore_diff > min_mean_zscore_diff;  %% mean zscore difference;
protomers_idx_out(sel1) = protomers_out_zscore_idx(sel1) + idx0;

%%% update the PO IDs (for outlier POs) in the following column in the
% the protomer_idx file
update_idx = cluster_idx_col + 1;
protomers_idx1 = protomers_idx(~exclu_sel, cluster_idx_col);
protomers_idx1(out_sel) = protomers_idx_out;
protomers_idx(exclu_sel, update_idx) = 1;
protomers_idx(~exclu_sel, update_idx) = protomers_idx1;



%%% fill the gap 2 - compare time
%%% purpose:
% identify the outlier POs according to closest neighbors in the time
% series (usually the PO of the same transporter domain one time-step
% before)

cluster_idx_col = update_idx;
% update_idx = cluster_idx_col + 1;    % assign IDs to the next column
protomers_idx_out = protomers_idx(:, cluster_idx_col);
state_time_threshold = 2;
for n = 1:nf
    % find the PO time trace
    sel_n = protomers_idx(:, 2) == n;
    protomers_idx_n = protomers_idx(sel_n, :);
    protomers_idx_opr = protomers_idx_n(:, [1 cluster_idx_col]);
    protomers_idx_opr_t = protomers_idx_opr(:, 1);
    protomers_idx_out1 = protomers_idx_opr(:, 2);

    % find the un-determined POs
    opr_sel = protomers_idx_out1 == 0;
    exclude_sel = protomers_idx_out1 <= 1 | isnan(protomers_idx_out1);  % excluded state

    % assign the unidentified PO state to its closest neighbot in the time
    % series (usually the one before)
    protomers_idx_opr_t_close = abs(protomers_idx_opr_t' - protomers_idx_opr_t);
    protomers_idx_opr_t_thresh = protomers_idx_opr_t_close <= state_time_threshold;
    protomers_idx_opr_t_close = protomers_idx_opr_t_close.*protomers_idx_opr_t_thresh;
    protomers_idx_opr_t_close = protomers_idx_opr_t_close(opr_sel, :);
    protomers_idx_opr_t_close(:, exclude_sel) = 0;
    protomers_idx_opr_t_close(protomers_idx_opr_t_close==0) = state_time_threshold + 1;
    [protomers_idx_opr_t_close_min, protomers_idx_opr_t_close_min_idx] = min(protomers_idx_opr_t_close, [], 2);
    protomers_idx_out_n = protomers_idx_out1(protomers_idx_opr_t_close_min_idx);
    protomers_idx_out_n(protomers_idx_opr_t_close_min == state_time_threshold + 1) = nan;
    protomers_idx_out1(opr_sel) = protomers_idx_out_n;
    protomers_idx_out(sel_n) = protomers_idx_out1;
end
protomers_idx(:, update_idx) = protomers_idx_out;
end