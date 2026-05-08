function [xo, fo, gridMask] = newGrid2d_v2(xi, fi, th, Nd, iw)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    xi,fi,th,Nd,iw
end

arguments (Output)
    xo, fo, gridMask
end
N = size(fi,1); %Length of one side of square grid
fo = zeros(size(fi));

%Initalize data storage matrices
a_data = zeros(N,N,Nd);
H_data = zeros(size(a_data));
V_data = zeros(size(a_data));
D_data = zeros(size(a_data));

data = fi;

for idecomp =1:Nd
    % fprintf("Decomp level: %i\n", idecomp)
    Ndim = N / (2^idecomp-1);
    [data, H_data(:,:,idecomp), V_data(:,:,idecomp), D_data(:,:,idecomp)]...
        = getApproxAndDetails(data,'haar', idecomp,N);
end

% Stores whether gridpoint is turned ON (1) or OFF (0)
% Initalize Mask
gridMask = false(N, N);

% 1. Keep the grid at the coarsest user-defined scale
step_coarse = 2^Nd;
gridMask(1:step_coarse:end, 1:step_coarse:end) = true;

% 2. Check detail coefficients via logical indexing
for idecomp = 1:Nd
    step_decomp = 2^idecomp;
    
    % Sum the coefficients for the entire matrix at once
    coefSum = abs(H_data(:,:,idecomp)) + abs(V_data(:,:,idecomp)) + abs(D_data(:,:,idecomp));
    
    % Find all points that exceed the threshold
    threshold_mask = coefSum > th;
    
    % Create a boolean grid that isolates the correct scale
    scale_grid = false(N, N);
    scale_grid(1:step_decomp:end, 1:step_decomp:end) = true;
    
    % Turn ON points that are BOTH on the correct scale AND exceed the threshold
    gridMask(threshold_mask & scale_grid) = true;
end

plotGridmask =0;
if plotGridmask
    % gridMask
    figure; imagesc(gridMask);
    colormap gray
end
xo=1;
end