clc, clear all, close all
dwtmode('sp0');
%Following paper describing 1D WOFD
%% Input Variables
N = 2^4;
fi = ones(N+1,1);
% xi = ones(N+1,1);
xi = linspace(-30,30, length(fi));
mid = round(length(fi)/2);
fi(mid:end) = 0;
% fi = exp(-1 * (xi*0.25).^2);
Nd = 1; % Maximum number of wavelet coefficients
iw =0; % Width of refinement stencil (decides how far forward and backward to look when deciding if grid point at location 0 should be used
th =0.00; % Needs to be tuned finely to see the grid refinement
L=4;

%% Program
[xo, fo, No] = newgrExample(xi, fi,L, N, th, Nd, iw);

figure;
plot(xi, fi);
hold on; plot(xo(1:No),fo(1:No),'--o'); hold off;
legend("Input", "Output");
disp("Length xi: " + length(xi));
disp("Length xo: " + No);