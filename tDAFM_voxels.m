%%%%%%
%%% Allocate detections into the voxels volume space
%%% Input: 
%%%        1: b: discretized detection 3D coordinates
%%%        2: bin_num: total bin numbers
%%%        3: showprogress: a flag to show progress. Default: 1
%%% Output:
%%%        1: voxels: the voxels volume space of detections
%%%%%%

function voxels = tDAFM_voxels(b, bin_num, showprogress)
[d1, d2, d3] = size(b);
voxels = zeros(d1, d2, bin_num);

edges = {1:1:d1 1:1:d2};
[CC, RR, ~] = meshgrid(1:d2, 1:d1, 1:d3);
for k = 1:bin_num
    if showprogress
        k + "/" + bin_num
    end
    sel_k = b == k;
    if sum(sel_k(:)) > 0
        CC_k = CC(sel_k);
        RR_k = RR(sel_k);
        voxels(:, :, k) = hist3([RR_k CC_k], edges);
    end
end
end