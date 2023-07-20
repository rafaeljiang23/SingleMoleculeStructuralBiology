close all

%% Sorting IFS and OFS
%%% sort IFS and OFS POs
% OFS POs are assigned an ID of 1
% IFS POs are assigned an ID of 2

%%% select image characteristics (parameter)
parameter_idx = [7 11];   % only consider height and volume for sorting OFS and IFS
X = protomers_sel_stats(:, parameter_idx);
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
input = score;
group_proposed = 2;    % IFS and OFS at this stage 

gmfit = fitgmdist(input, group_proposed, 'Replicates', 5, 'Start', 'randSample');
gmfit_idx = cluster(gmfit,input); % Cluster index
cluster_idx = gmfit_idx;
%% display cluster result and the energy landscape
%%% calculate energy landscape
num_bin = 600;
[EE, XX, YY] = energylanscape(num_bin, score, 15, 2);

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
% the initial two state (IFS vs OFS) ID is recorded at column 3.
update_idx = 3;   

% note that the OFS cluster should be labeled as group 1 and the 
% IFS cluster group 2 using the following codes if necessary
cluster_idx(cluster_idx == 1) = 3;     % temporary storage
cluster_idx(cluster_idx == 2) = 1;
cluster_idx(cluster_idx == 5) = 2;

cluster_idx_all = zeros(protomer_tot, 1);
cluster_idx_all(~out_sel) = cluster_idx;
cluster_idx_all(out_sel) = nan;    % all omitted point is assigned to nan
protomers_idx(:, update_idx) = cluster_idx_all;

%% Identify the outliers
%%% assign ID values to the outlier POs with nan

cluster_idx_col = update_idx;   % column recording the initial cluster IDs

% down_sel records all POs in the OFS state, and not created at this stage
% a psedo-file is created for the operation

down_sel = zeros(protomer_tot, 1);
down_sel = down_sel == 1;

% input all POs
X2 = protomers_sel_stats(:, parameter_idx);
% score2 = normalize(X2)* inv(coeff);   % project all POs onto the pc1-pc2 space 

%%% fill the nan values
nf = max(protomers_idx(:, 2));
protomers_idx = fill_protomers_idx_nan(update_idx, protomers_idx, score2, down_sel, out_sel, 0, 0, nf);

%%% protomers_idx file at this stage N x 4
%%% N: POs
%%% columes: 1: frame ti. 2: original protomer pj
%%%          3-4: two-state IDs
%%%          3: initial ID from PCA
%%%          4: final ID (see fill_protomers_idx_nan.m function)
%% display IFS and OFS clusters
% display_idx = 4;
% for i = 1:max(cluster_idx_all)
%     MIJ.createImage(protomers_nf(:, :, protomers_idx(:, display_idx)==i));
%     MIJ.createImage(mean(protomers_nf(:, :, protomers_idx(:, display_idx)==i), 3));
% end

