close all;


%% create LAFM detections (3D local maxima matrix) 
%%% purpose: find local maxima from the raw data as LAFM detections
%%% define useful parameters
radius = 70;  % effective PO object radius, Unit: pixel

%%%
[d2, d1, d3] = size(data);
detect_tot = 0;    % total detection count

local_max = zeros(d1*d2*d3, 5);    
% local maxima matrix, size N x 5
% N: LAFM detections
% columns: recording local maxima info
%          1-2: local maxima coordiates x and y
%          3: local maxima height value
%          4: local maxima prominence (significance)
%          5: local maxima frame

for i = 1 : d3
    i
    data_frame = data(:, :, i);
    %find local maximum
    [data_local_max_c, prominence_c] = islocalmax(data_frame, 1 , 'MinProminence', 0);
    [data_local_max_r, prominence_r] = islocalmax(data_frame, 2 , 'MinProminence', 0);
    data_local_max_rc = data_local_max_c & data_local_max_r;   %local maximum positions
    min_prominence_rc = min(prominence_c, prominence_r);    %local maximum prominence
    data_local_max_height = data_frame .* (data_local_max_rc > 0);     %local maximum height
    % grow local maxima matrix
    local_max = grow_local_max_matrix(local_max, d2, d1, detect_tot, data_local_max_rc, data_local_max_height, min_prominence_rc, i, radius);
    detect_tot = detect_tot + sum(data_local_max_rc > 0, "all");
end
local_max = local_max(local_max(:, 5)>0, :);


%% functions
%%%%%%
%%% Grow local maxima matrix for LAFM map construction
%%% Input: 
%%%        1. M_max: local maxima matrix 
%%%        2-3. lx, ly: data dimension
%%%        4. c: local maxima counter
%%%        5. M_local_max: local maxima positions in the data frame
%%%        6. M_local_height: local maxima height values
%%%        7. M_prominence: local maxima prominence values
%%%        8. f: current frame
%%%        9. radius: effective PO object radius, used to eliminate
%%%           out-ranged local maxima
%%% Output:
%%%        1: M_max: updated local maxima matrix 
%%%%%%
function M_max = grow_local_max_matrix(M_max, lx, ly, c, M_local_max, M_local_height, M_prominence, f, radius)
x_c = lx/2;
y_c = ly/2;
i = 1;
for x = 1 : lx
    for y = 1 : ly
        if (x - x_c)^2 + (y - y_c)^2 < radius ^ 2
            if M_local_max(y, x) > 0
                M_max(c + i, 1) = x;   %local max x
                M_max(c + i, 2) = y;   %local max y
                M_max(c + i, 3) = M_local_height(y, x);  %local max height
                M_max(c + i, 4) = M_prominence(y, x);  %local max prominence
                M_max(c + i, 5) = f;  %local max frame
                i = i + 1;
            end
        end
    end
end
end

