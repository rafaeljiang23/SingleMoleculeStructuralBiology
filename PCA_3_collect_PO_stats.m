close all

%% collect statistics
%%% define useful parameters
hisbins = linspace(0, 4,200) ;   % for calculating image entropy
sigma3 = 30;  % estimated size of a protomer. Unit: pixel
lower_limit = -1; % lower height boundary for PO pixel vlaues
upper_limit = 3; % upper height boundary for PO pixel vlaues

%%%
% PO stats boundary filter
radius_filter1 = zeros(lx, ly);
radius_in1 = radius_in;
radius_ot1 = floor(lx/2)-1;

for x = 1:lx
    for y = 1:ly
        radius_temp = (x - xc)^2 + (y - yc)^2;
        radius_filter1(x, y) = radius_temp > radius_in1 ^2 & radius_temp < radius_ot1 ^2;
    end
end

% combine POs in the same frame (should give back the cropped raw data)
% subjected to later operations
protomers_combine = zeros(lx, ly, d3);
for t = 1 : d3
    for i = 1 : nf
    protomers_combine(:, :, t) = protomers_combine(:, :, t) + protomers(:, :, (t-1)*nf+i) - protomers(1, 1, (t-1)*nf+i);
    end
end


%%% create PO statistics file
% PO pixel stats: N x M matrix for PCA 
% N: POs
% M: PO image charateristics: (cols 1-5 from PO COM stats)
%           1-2: PO COM cartisian coordinates x and y.
%           3-4: PO COM polar coordinates theta and rho.
%           5: PO COM normalized theta
%           6: PO ROI histogram kurtosis. 
%           7: PO ROI delta height.
%           8: PO ROI mean gradient.
%           9: PO ROI histogram skewness 
%           8: PO ROI image entropy  
%           9: PO ROI delta volume 
protomers_sel_stats = zeros(lz, 20);
protomers_sel_stats(:, 1:5) = protomers_sel_locs;  


%%% create PO ROI file
% pixels within PO ROI are subjected to PCA
protomers_roi = protomers;
% protomers_roi2 = protomers;

[lx, ly, lz] = size(protomers);
for p = 1:lz
    p
    % get PO objects
    pt = floor((p-1)/nf)+1;
    obj_frame = protomers_combine(:, :, pt);
    obj_frame(obj_frame<lower_limit) = lower_limit;
    obj_frame(obj_frame>upper_limit) = upper_limit;
    obj_frame = reshape(obj_frame, lx ,ly);
    % get PO area
    obj_sel = protomers_loc(:, :, p);
    obj_sel = reshape(obj_sel, lx, ly);
    obj_sel = imgaussfilt(double(obj_sel), sigma3)*(sigma3^2*2*pi);
    obj_sel = obj_sel > 0.68 & radius_filter1;
    num_pixel = sum(obj_sel, "all");
    % get PO pixels
    X = obj_frame .*(obj_sel);
    protomers_roi(:, :, p) = X;
    X1 = obj_frame(obj_sel);
    X1 = rmoutliers(X1, 'median', 'ThresholdFactor', 5); %remove outliers
    background = quantile(X1, 0.05, "all");   % background value 
    [Gmag,Gdir] = imgradient(obj_frame,'sobel');   % PO image gradient
    %%% update PO stats
    protomers_sel_stats(p, 6) = kurtosis(X1);  % PO histogram kurtosis
    protomers_sel_stats(p, 7) = quantile(X1, 0.9, "all") - background;  % PO delta height 
    protomers_sel_stats(p, 8) = mean(sum(Gmag(obj_sel)));  % PO mean gradient   
    protomers_sel_stats(p, 9) = skewness(X1);  % PO histogram skewness  
    protomers_sel_stats(p, 10) = entropy2(X1, hisbins);    % PO image entropy  
    protomers_sel_stats(p, 11) = sum(X1, "all") - num_pixel*background;  % PO delta volume   

end

%% display PO ROI
protomers_roi_combine = zeros(lx, ly, d3);
for t = 1 : d3
    for i = 1 : nf
    protomers_roi_combine(:, :, t) = protomers_roi_combine(:, :, t) + protomers_roi(:, :, (t-1)*nf+i);
    end
end

% MIJ.createImage(protomers_roi);
% MIJ.createImage(protomers_roi_combine);

%% functions
%%%%%%
%%% Calculating image entropy with user defined histogram bins
%%% Input: 
%%%        1: X: input image 
%%%        2: bins: histogram bins
%%% Output:
%%%        1: e: image entropy
%%%%%%
function e = entropy2(X, bins)
count = histcounts(X, bins);
p = count./sum(count);
p = p(p>0);
e = -sum(p.*log2(p));
end