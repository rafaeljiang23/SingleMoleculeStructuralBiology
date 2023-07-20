close all

%% Construct (3D) LAFM map
%% select POs for LAFM construcction

%%% Note that at this stage, user should define a 'LAFM_sel' variable 
%%% according to the puspose of LAFM map construction. 
%%% POs can be selected using various riteria: 
%%% examples are: 1: cluster IDs. 
%%%               2: positions (measured by PCA score) on the pc1-pc2 space
%%%               3: activity (measured by their distance to an event on 
%%%                  the time-series (time-state traces)

%%% here the example to construct IFSo LAFM map is given:
% selecting POs that has a specific cluster ID
LAFM_sel = LAFM_detections_species(:, 1) == 2;  % cluster ID=2 corresponds
                                                % to IFSo POs

%% LAFM detection matrix
%%% allocate the selected LAFM detections into a particle stack 
% dimensions X-by-Y-by-N
% X and Y are real space plane dimensions
% N is the AFM frames

LAFM_input = LAFM_detections_species(LAFM_sel, :);
z_min = min(LAFM_input(:, 5));
[d1, d2, d3] = size(data);
detections = zeros(d1, d2, d3);
detections = detections + z_min;
for i = 1:sum(LAFM_sel)
    LAFM_input_i = LAFM_input(i, :);
    x = LAFM_input_i(3);
    y = LAFM_input_i(4);
    t = LAFM_input_i(7);
    h = LAFM_input_i(5);
    detections(y, x, t) = h;
end
detections(detections == z_min) = nan;
%% crop bad pixels
%%% defube useful parameters
xy_radius = 70;    % effective (3D) LAFM map radius
%%%

[XX,YY,ZZ] = meshgrid(1:d2, 1:d1,1:d3);
sel_radius = (XX - d1/2).^2 + (YY-d2/2) .^2 < xy_radius.^2;

detections = detections .* sel_radius;
detections(~sel_radius) = nan;
%% voxels construction
%%% construct a 3D volume space, and allocate detections into this space
% dimensions X-by-Y-by-H
% X and Y are real space plane dimensions
% H is real space height dimension

%%% define useful parameters 
nf = 3;   % oligomer molecular symmetry
z_min = -1.5;   % maximum z detection. Unit: nm
z_max = 1.5;    % minimum z detection. Unit: nm 
resolution_z = 0.01;  % target resolution in z dimension. Unit: nm/pixel; 
resolution_xy = 0.1;   % resolution in xy dimension (input data). Unit: nm/pixel
%%%

b = floor((detections - z_min)./resolution_z);
bin_num = floor((z_max - z_min)./resolution_z);
%%% allocate detections into the voxels volume space
voxels = tDAFM_voxels(b, bin_num, 1);

%%% apply molecular symmetry to the voxels volume space
voxels_nf = zeros(d1, d2, bin_num);
for i = 1:nf
    angle = (i-1)*(360/nf);
    voxels_nf = voxels_nf + imrotate(voxels, angle, "nearest", "crop");
end
voxels = voxels_nf;
voxels(:, :, 1:10) = 0;   % background noise

%% convolution with a point spreading function

%%% define useful parameters
sigma_xy = 1.4; 
sigma_z = 1.4;
%%%

sigma_xy = sigma_xy / (resolution_xy * 10);
sigma_z = sigma_z / (resolution_z * 10);

%%% customized gaussian point spreading function 
h = make_3D_LAFM_kernel1e(sigma_xy, sigma_z, resolution_xy, resolution_z);  %% shape: gauss z; psf: gauss xyz

%%% apply the psf to the voxels volume splace
% voxels_h is the detection density volume space
voxels_h = imfilter(voxels, h);

%%% apply molecular symmetry to the detection density volume space (voxels_h)
voxels_hs = voxels_h;
if nf > 1
for i = 2:nf
    angle = (i-1)*(360/nf);
    voxels_hs = voxels_hs + imrotate(voxels_h, angle, "bicubic", "crop");
end
end
voxels_hs = voxels_hs./sum(voxels_hs(:));

% normalize the density to z-scores
voxels_hsz = (voxels_hs - mean(voxels_hs(:)))./std(voxels_hs(:));

%% surface normalization
%%% purpose: to find the most likely surface of the detection density
%%% volume space (voxels_hsz). The normalized detecition density volume 
%%% space (voxels_hsn) should have a value 0-1 for each x-y cylinder, where
%%% a value of 1 represents the most likely hight at this x-y position.

voxels_hsn = voxels_hs./max(voxels_hs, [], 3);
voxels_hsn(voxels_hsn < 0) = 0;
voxels_hsn(isnan(voxels_hsn)) = 0;
final_radius_ratio = 0.99;
[d11, d22, d33] = size(voxels_hsn);
crop_radius = xy_radius * final_radius_ratio;

[X2, Y2, Z2] = meshgrid(1:d11, 1:d22, 1:d33);
sel_radius_crop = (X2-d11/2).^2 + (Y2-d22/2).^2 <crop_radius.^2;

voxels_hsn = voxels_hsn .*sel_radius_crop;
voxels_hsn(~sel_radius_crop) = min(voxels_hsn(:));

%% Display (3D) LAFM construction result
%%% NOTE that these volume space (saved as tiff) can be directly 
%%% visualized in ChimeraX. Measurements are also directly made on
%%% these volume space data.

% MIJ.createImage(voxels);
% MIJ.createImage(voxels_hsz);
% MIJ.createImage(voxels_hsn);