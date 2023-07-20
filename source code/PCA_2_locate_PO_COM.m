close all

%% locate protomer observation (PO) center of mass (COM) in the PO objects
%%% define useful parameters
sigma1 = 10;      % sigma value for gauss filtering, purpose: denoise the raw data for locating PO COM. Unit: pixel
% sigma2 = 5;

radius_in = 15;     % radius for the inner boundary of the PO COM. Unit: pixel
radius_ot = 50;     % radius for the outer boundary of the PO COM. Unit: pixel
thresh_up = 2.8;     % height boundary of the PO object. Unit: pixel

%%% 
[d1, d2, d3] = size(data);
[lx, ly, lz] = size(protomers);
[xx, yy] = meshgrid(1:lx, 1:ly);
xx = xx(:);
yy = yy(:);
xc = lx/2;
yc = ly/2;

% PO COM boundary filter
radius_filter = zeros(lx, ly);
for x = 1:lx
    for y = 1:ly
        radius_temp = (x - xc)^2 + (y - yc)^2;
        radius_filter(x, y) = radius_temp > radius_in ^2 & radius_temp < radius_ot ^2;
    end
end

%%% create PO location files
% PO COM stats (location stats): N x M matrix
% N: POs
% M: PO COM stats: 1-2: cartisian coordinates x and y. 3-4: polar
% coordinates theta and rho. 5: normalized theta
protomers_sel_locs = zeros(lz, 5);

protomers_loc = zeros(lx, ly, d3);     % PO COM in PO objects
%%% locate POs in the PO objects
for p = 1 : lz
    p
    obj_frame = protomers(:, :, p);
    obj_frame = reshape(obj_frame, lx ,ly);
    % set out-ranged pixels to the min value of the image
    obj_frame(~radius_filter) = min(obj_frame(:)); 
    obj_frame(obj_frame > thresh_up) = min(obj_frame(:));  

    % select valid height values
    X = obj_frame(obj_frame> min(obj_frame(:)));
    thresh = quantile(X, 0.75);   % 75% quantile of the valid height values
    obj_frame_bw = obj_frame > thresh;   % possible PO pixels: height values > 75% quantile 

    % connect possible PO pixels, raw PO area
    bw = imgaussfilt(double(obj_frame_bw), sigma1);
    bw = bw * (sigma1^2*2*pi);
    bw = bw > 0.68;

    % find the PO center of mass (COM)
    D = bwdist(~bw);
    M = islocalmax(D, 1) & islocalmax(D, 2);
    max_value = max(D(M), [], "all");
    M_sel = D == max_value;
    if sum(M_sel, "all") > 1
        % if multiple PO COM found, the PO COM should be the one with the
        % larget height in the PO object
        max_value2 = max(obj_frame(M_sel), [], "all");
        M_sel = obj_frame == max_value2;
    end
    protomers_loc(:, :, p)  = M_sel;

    % cartesian coordinates of the PO COM
    M_sel_x = xx(M_sel);
    M_sel_y = yy(M_sel);
    % polar coordinates of PO COM
    [M_sel_theta, M_sel_rho] = cart2pol(M_sel_x - xc, M_sel_y - yc);

    % update PO location stats
    protomers_sel_locs(p, 1) = M_sel_x(1) - xc;
    protomers_sel_locs(p, 2) = M_sel_y(1) - yc;
    protomers_sel_locs(p, 3) = M_sel_theta(1);
    protomers_sel_locs(p, 4) = M_sel_rho(1);
end

%%% normalize PO theta (angle in the polar coordinates)
protomers_sel_locs(:, 5) = rem(pi + protomers_sel_locs(:, 3), 2*pi/nf);
protomers_sel_locs(:, 5) = protomers_sel_locs(:, 5) - mean(protomers_sel_locs(:, 5));


%% display PO objects
protomers_loc_dis = imgaussfilt(protomers_loc, 5);
protomers_loc_combine_dis = zeros(lx, ly, d3);
protomers_loc_combine = zeros(lx, ly, d3);

for t = 1 : d3
    for i = 1 : nf
    protomers_loc_combine_dis(:, :, t) = protomers_loc_combine_dis(:, :, t) + protomers_loc_dis(:, :, (t-1)*nf+i);
    protomers_loc_combine(:, :, t) = protomers_loc_combine(:, :, t) + protomers_loc(:, :, (t-1)*nf+i);
    end
end

% MIJ.createImage(protomers_loc_combine_dis);