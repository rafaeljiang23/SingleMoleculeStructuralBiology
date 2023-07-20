close all

%% organize LAFM detections
%%% purpose: assign PO tipj ID and cluster ID to the detections for LAFM
%%% map constructions

%% step 1: assign PO object ID (tipj) to the detections
[~, ~, d3] = size(data);
%%%
LAFM_detections_protomers = organize_protomers_local_max(d3, protomers_stat, local_max);
% LAFM_detections_protomers matrix, size N x 6
% N: LAFM detections
% columns: 1: PO object ID (pj)
%          2-6 are columns 1-5 from local_max matrix:
%          2-3: local maxima coordiates x and y
%          4: local maxima height value
%          5: local maxima prominence (significance)
%          6: local maxima frame
% note that this matrix does not dependent on the cluster ID

%% step 2: assign PO cluster ID to the detections
%%% purpose: allow the LAFM map to be constructed with user-defined cluster ID

%%% define useful parameter
% colume in protomers_idx that corresponds to final PO cluster id
cluster_idx_col = 6;     % three-state clustering results

%%%
LAFM_detections_species = organize_species_local_max(cluster_idx_col, d3, protomers_idx, LAFM_detections_protomers);
% LAFM_detections_species matrix, size N x 7
% N: LAFM detections
% columns: 1: PO cluster ID (OFS:1, IFSo:2, IFSc: 3)
%          2-7 are columns 1-6 from LAFM_detections_protomers matrix:
%          2: local maxima PO object ID (pj)
%          3-4: local maxima coordiates x and y
%          5: local maxima height value
%          6: local maxima prominence (significance)
%          7: local maxima frame
% note that this matrix does not dependent on the cluster ID


%% functions
%%%%%%
%%% Grow LAFM_detections_protomers matrix for LAFM map construction
%%% Input: 
%%%        1. d3: total frame number
%%%        2. P_stat: protomers_stat file, recording pixel info
%%%        3. M_max: local maxima matrix, recording local maxima info
%%% Output:
%%%        1. M_pro: LAFM_detections_protomers matrix 
%%%%%%
function M_pro = organize_protomers_local_max(d3, P_stat, M_max)
seed = RandStream('mlfg6331_64'); 
[num_local_max, ~] = size(M_max);
M_pro = zeros(num_local_max, 1);
for t = 1: d3
    t
    P_stat_t = P_stat(P_stat(:, 1) == t, :);    % submatrix of P_stat with time == t
    sel_t = M_max(:, 5) == t;
    M_max_t = M_max(sel_t, :);
    M_prot_t = M_pro(sel_t);     % submatrix of M_protomers with time == t
    for i = 1 : sum(sel_t)
        coor = M_max_t(i, 1:2);         % coor
        idx = (P_stat_t(:, 3) == coor(2)) & (P_stat_t(:, 4) == coor(1));
        if sum(idx) == 1     % not boundary pixels
           M_prot_t(i, 1) = P_stat_t(idx , 2);   % update protomers_id
        elseif sum(idx) > 1    
           % boundary pixels are randomly assigned to a group 
           % note that these pixels are usually noise
           P_stat_t1 = P_stat_t(idx , 2);
           prandom = randsample(seed, P_stat_t1, 1);
           M_prot_t(i, 1) = prandom;   % update protomers_id
        end
    end
    M_pro(sel_t) = M_prot_t;    % update M_protomers
end
M_pro = [M_pro M_max];
end

%%%%%%
%%% Assign PO cluster ID to the detections for LAFM map construction
%%% Input: 
%%%        1. P_idx_col: column in protomers_idx that corresponds to final PO cluster id
%%%        2. d3: total frame number
%%%        3. P_idx: protomers_idx file
%%%        4. M_pro: LAFM_detections_protomers matrix 
%%% Output:
%%%        1. M_spe: LAFM_detections_species matrix 
%%%%%%
function M_spe = organize_species_local_max(P_idx_col, d3, P_idx, M_pro)
nf = max(P_idx(:, 2));
[num_local_max, ~] = size(M_pro);
M_spe = zeros(num_local_max, 1);
for t = 1:d3
    t
    for p = 1:nf
        % find the PO object tipj in M_pro
        sel_t_p = (M_pro(:, 1) == p) & (M_pro(:, 6) == t);
        % find the PO object tipj in P_idx
        sel_t_p_idx = (P_idx(:, 1) == t) & (P_idx(:, 2) == p);
        % assign PO cluster ID
        M_spe(sel_t_p) = P_idx(sel_t_p_idx, P_idx_col);
    end
end
M_spe = [M_spe M_pro];
end