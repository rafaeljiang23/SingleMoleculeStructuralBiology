%%%%%%
%%% Create energy landscape
%%% Input: 
%%%        1: num_bins: number of histogram bins for binning PCA scores 
%%%        2: scores: PCA scores
%%%        3: sigma: sigma value for gaussian bluring the landscape map
%%%        4: mincount: background count. Default: 2
%%% Output:
%%%        1: EE: energy (2D matrix)
%%%        2: XX: positions in axis 1 (2D matrix)
%%%        3: YY: positions in axis 2 (2D matrix)
%%%%%%

function [EE, XX, YY] = energylanscape(num_bin, score, sigma, mincount)
xbins_edg = linspace(min(score(:, 1))-std(score(:, 1)), max(score(:, 1))+std(score(:, 1)), num_bin);
ybins_edg = linspace(min(score(:, 2))-std(score(:, 2)), max(score(:, 2))+std(score(:, 2)), num_bin);
h = histogram2(score(:, 1), score(:, 2), xbins_edg, ybins_edg, "Visible","off");
xbins_ctr = 0.5.*(xbins_edg(1:end-1) + xbins_edg(2:end));
ybins_ctr = 0.5.*(ybins_edg(1:end-1) + ybins_edg(2:end));
[XX, YY] = meshgrid(xbins_ctr,ybins_ctr);
XX = XX';
YY = YY';
count = h.BinCounts;
count = imgaussfilt(count, sigma);
p = count * (sigma^2*2*pi);
p(p<mincount) = mincount;
EE = -log(p);
EE = EE - max(EE, [], "all");
end