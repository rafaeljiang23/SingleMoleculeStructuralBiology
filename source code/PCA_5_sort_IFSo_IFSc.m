close all

%% Sorting IFSc and IFSo
%%% sort IFS and OFS POs
% IFS open (IFSo) POs are assigned an ID of 2
% IFS closed (IFSc) POs are assigned an ID of 3

%%% select image characteristics (parameter)
cluster_idx_col = 4;    % column recording the final 2-state PO cluster IDs
parameter_idx = 4:11;   % use all parameters, note that 1-3 are not normalized (protomer dependent) and thus excluded
down_sel = protomers_idx(:, cluster_idx_col) == 1;    % select all IFS POs (protomer_idx column 5 from PCA_4 should be 1 in this case)

protomers_up_stats = protomers_sel_stats(~down_sel, :);
X = protomers_up_stats(:, parameter_idx);
[protomer_tot, ~] = size(protomers_sel_stats);

%%% remove outliers
% outliers should be removed for better PCA results. The ID of the outlier
% POs is later assigned according to their relative distance to the cluster
% centers (higher priority) and their positions in the time series.
[X, out_sel] = rmoutliers(X, 'median', 'ThresholdFactor', 3);

%%% PCA
[coeff, score, latent] = pca(X, 'VariableWeights','variance');


%%% display result
figure();
hold on
scatter(score(:, 1), score(:, 2));
hold off

%% clusters 
%%% Gaussian Mixture Models (GMM) for clustering
% note that the clustering result may slightly vary due to the initial
% condition variation for the 'fitgmdist' function.
input = score(:, 1:2);     % sort POs in the pc1-pc2 space
group_proposed = 2;    % IFSo and IFSc at this stage 

gmfit = fitgmdist(input, group_proposed, 'Replicates', 5, 'Start', 'randSample');
gmfit_idx = cluster(gmfit,input); % Cluster index
cluster_idx = gmfit_idx;
%% display cluster result and the energy landscape
%%% calculate energy landscape
num_bin = 600;
[EE, XX, YY] = energylanscape(num_bin, score(:, 1:2), 15, 2);

figure()
hold on
contour(XX, YY, EE, 50);
for i = 1:max(cluster_idx)
    cluster_idx_sel = cluster_idx == i;
    sum(cluster_idx_sel)
    scatter(score(cluster_idx_sel, 1), score(cluster_idx_sel, 2), "filled", "MarkerFaceAlpha", 0.3);
    xlabel("pca 1")
    ylabel("pca 2")
end
hold off
legend

%% cluster idx update
%%% the cluster result (PO ID) is recorded in the protomer_idx file
% the initial two state (IFSo vs IFSc) ID is recorded at column 6. The OFS
% POs are assigned to 1 automatically at this stage.
update_idx = cluster_idx_col + 1;   

% note that the IFSo cluster should be labeled as group 1 and the 
% IFS cluster group 2 using the following codes if necessary. A value of 1
% is added to all IFS POs to adjust their IDs to 2 and 3 respectively
% cluster_idx(cluster_idx == 1) = 5;     % temporary storage
% cluster_idx(cluster_idx == 2) = 1;
% cluster_idx(cluster_idx == 5) = 2;

cluster_idx_all = zeros(sum(~down_sel), 1);
cluster_idx_all(~out_sel) = cluster_idx + 1;   % add 1 to the IFS PO IDs
cluster_idx_all(out_sel) = nan;
protomers_idx(down_sel, update_idx) = 1;
protomers_idx(~down_sel, update_idx) = cluster_idx_all;


%% Identify the outliers
%%% assign ID values to the outlier POs with nan

cluster_idx_col = update_idx;   % column recording the initial cluster IDs

% input all POs
score2 = protomers_sel_stats(:, parameter_idx);
% score2 = normalize(score2)* inv(coeff);  % project all POs onto the pc1-pc2 space 

%%% fill the nan values
nf = max(protomers_idx(:, 2));
protomers_idx = fill_protomers_idx_nan(cluster_idx_col, protomers_idx, score2, down_sel, out_sel, 1, 0, nf);

%%% protomers_idx file at this stage N x 6
%%% N: POs
%%% columes: 1: frame ti. 2: original protomer pj
%%%          3-4: two-state IDs: ID-1: OFS. ID-2: IFS.
%%%          3: initial cluster ID from PCA
%%%          4: final cluster ID (see fill_protomers_idx_nan.m function)
%%%          5-6: three-state IDs: ID-1: OFS. ID-2: IFSo. ID-3: IFSc.
%%%          5: initial cluster ID from PCA
%%%          6: final cluster ID (see fill_protomers_idx_nan.m function)

%% display IFSo and IFSc clusters
% display_idx = 6;
% for i = 1:max(cluster_idx_all)
%     MIJ.createImage(protomers_nf(:, :, protomers_idx(:, display_idx)==i));
%     MIJ.createImage(mean(protomers_nf(:, :, protomers_idx(:, display_idx)==i), 3));
% end

