close all;

%% create protomer observation (PO) object files
%%% define useful parameters
nf = 3;    % oligomer molecular symmetry value
radius = 70;  % radius for creating POs. For determining pixels to be
              % included in PO objects. Unit: pixel 
ang0 = 60;   % angle zero for creating POs. For determining boundries 
             % between neighboring PO objects in the oligomer. Unit: degree 
ang_tolerance_deg = 30;    % angle allowance for boundary pixels. Pixels
                           % within this allowance near the boundary are 
                           % counted for both object. Unit: degree
[d1, d2, d3] = size(data);

%%% PO files
protomers = zeros(d2, d1, nf*d3);    % PO objects
protomers_nf = zeros(d2, d1, nf*d3);     % symmetrized PO objects
protomers_stat = zeros(d2*d1*d3, 4);     % record PO object info 
                                         % ti, pj, x, y 
                                         % (for LAFM constructions)
idx = 0;  
for i = 1 : d3
    i
    data_frame = data(:, :, i);
    [prot, prot_nf, stat] = pick_protomer(data_frame, nf, ang0, radius, ang_tolerance_deg);
    [icr, ~] = size(stat);
    protomers_stat(idx+1: idx+icr, 2:4) = stat(:, :);
    protomers_stat(idx+1: idx+icr, 1) = i;
    protomers(:, :, nf*(i-1)+1: nf*(i-1)+nf) = prot(:, :, 1:nf);
    protomers_nf(:, :, nf*(i-1)+1: nf*(i-1)+nf) = prot_nf(:, :, 1:nf);
    idx = idx + icr;
end
protomers_stat = protomers_stat(1: idx, :);
%%% croping the black pixels
protomers = stack_crop(protomers, radius);
protomers_nf = stack_crop(protomers_nf, radius);

%% creating PO index
%%% PO index file records the protomer identity in PCA
% PO index tipj characterizes protomer j in frame i 
protomer_tot = d3 * nf;
protomers_idx = zeros(protomer_tot, 10);
idx = 0:protomer_tot-1;

protomers_idx(:, 1) = floor(idx/nf) + 1;    % PO index ti
protomers_idx(:, 2) = rem(idx, nf) + 1;   % PO index pj
%% display PO objects
% MIJ.createImage(protomers)
% MIJ.createImage(protomers_combine)
% MIJ.createImage(protomers_nf)
%% functions

%%%%%%
%%% Pick PO objects from the oligomers
%%% Input: 
%%%        1: oligo: input oligomer object (single frame)
%%%        2: nf: oligomer molecular symmetry value
%%%        3: ang0_deg: angle zero for creating POs. For determining boundries between neighboring PO objects in the oligomer. Unit: degree 
%%%        4. radius: radius for creating POs. For determining pixels to be included in PO objects. Unit: pixel 
%%%        5. ang_tolerance_deg: angle allowance for boundry pixels. Boundary pixels are counted for both objects. Unit: degree
%%% Output:
%%%        1: prot: protomer object file (PO)
%%%        2: nf_prot: symmetrized protomer object file
%%%        3: stat: PO object info file, recording ti, pj, x, y values of PO pixels
%%%%%%
function [prot, nf_prot, stat] = pick_protomer(oligo, nf, ang0_deg, radius, ang_tolerance_deg)     
nf_rad = deg2rad(360/nf);
ang0_rad = deg2rad(ang0_deg);
ang_tolerance_rad = deg2rad(ang_tolerance_deg);
[ly, lx] = size(oligo);
% find oligomer center
x_c = lx/2;
y_c = ly/2;
% create PO file
prot = zeros(ly, lx, nf);
prot = prot + min(oligo, [], "all");
stat = zeros(ly*lx, 3);    % record PO object ID, x and y
% allocate pixels to the PO objects
idx = 0;
for x = 1 : lx
    for y = 1 : ly
        if (x - x_c)^2 + (y - y_c)^2 < radius ^ 2
            idx = idx + 1;
            stat(idx, 2:3) = [x, y];   % cartesian coor
            [theta_rad, ~] = cart2pol(x-x_c, y-y_c);
            theta_rad =  theta_rad + 2*pi;
            inc = 0;
            idx0 = 0;
            % find protomer id
            while theta_rad > -ang0_rad + nf_rad 
                theta_rad = theta_rad - nf_rad;
                inc = inc + 1;
                if abs(theta_rad + ang0_rad - nf_rad) < ang_tolerance_rad   
                    % pixel near protomer boundary
                    idx0 = idx0 + 1;
                    idx1 = idx + idx0;
                    prot_id_boundary = rem(inc, nf) + 1;
                    stat(idx1, 2:3) = [x, y];
                    stat(idx1, 1)  = prot_id_boundary;    % PO pj

                end
            end
            prot_id = rem(inc, nf) + 1;
            stat(idx, 1)  = prot_id;    % PO pj
            % allocate pixels
            prot(x, y, prot_id) = oligo(x, y);
            idx = idx + idx0;
        end
    end
end
stat = stat(1: idx, :);
% create symmetrized PO file
nf_prot = stack_nfold_fill(prot, nf);
end


%%%%%%
%%% Create symmetrized objects
%%% Input: 
%%%        1: in: input object 
%%%        2: nf: symmetry value
%%% Output:
%%%        1: out: symmetrized object
%%%%%%
function out = stack_nfold_fill(in, nf)
[ly, lx, lz] = size(in);
out = zeros(ly, lx, lz);
for t = 1:lz
    for i = 1:nf
        nf_deg = (i - 1) * (360/nf);
        in0 = in(:, :, t);
        in0 = reshape(in0, [ly lx]);
        inr = imrotate(in0, nf_deg, 'bicubic', 'crop');
        out(:, :, t) = out(:, :, t) + inr;
    end
    out(:, :, t) = out(:, :, t) - (nf-1)* min(in0, [], "all");
end
end

%%%%%%
%%% Crop out-ranged pixels 
%%% Input: 
%%%        1: in: input object
%%%        2: radius: effective range (circle radius)
%%% Output:
%%%        1: out: cropped object
%%%%%%
function out = stack_crop(in, radius)
[ly, lx, lz] = size(in);
x_c = lx/2;
y_c = ly/2;
out = zeros(2*radius+1, 2*radius+1, lz);
for t = 1:lz
        in0 = in(:, :, t);
        in0 = reshape(in0, [ly lx]);
        out(:, :, t) = imcrop(in0, [floor(y_c)-radius+1 floor(x_c)-radius+1 2*radius 2*radius]);
end
end