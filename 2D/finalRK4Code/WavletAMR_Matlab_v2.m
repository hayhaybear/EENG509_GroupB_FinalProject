clc, clear all, close all
dwtmode('sp0');
%Following paper describing 1D WOFD
%% Input Variables
N = 2^6;
fi = ones(N+1,1);
% xi = ones(N+1,1);
xi = linspace(-30,30, length(fi));
mid = round(length(fi)/2);
fi(mid:end) = 0;
% fi = exp(-1 * (xi*0.25).^2);
Nd = 1; % Maximum number of wavelet coefficients
iw =1; % Width of refinement stencil (decides how far forward and backward to look when deciding if grid point at location 0 should be used
th =0.01; % Needs to be tuned finely to see the grid refinement
L=4;
%{
%% Testing newGrid2d_v2
xi = 1;
fi = rand(N,N);
xi = linspace(-5,5, N); yi = linspace(-5,5, N);
[XGrid, YGrid] = meshgrid(xi,yi);
fi = 20*exp(-1 * (XGrid.^2 + YGrid.^2));
disp("Max")
max(fi(:))
fi(N/2:end,:) = fi(N/2:end,:) +20;

Nd =3;
figure;surf(fi);
[xo, fo, gMask] = newGrid2d_v2(xi, fi, th, Nd, iw);
fprintf("Uncompressed Numpts: %i\t Compressed Numpts: %i\n", N^2, sum(gMask(:)));