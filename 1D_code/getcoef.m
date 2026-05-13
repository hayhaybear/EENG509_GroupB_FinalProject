function [h, g] = getcoef(L)
% Equivalent MATLAB version of the Fortran subroutine getcoef

Lmax = 8;

% Preallocate
h = zeros(Lmax,1);
g = zeros(Lmax,1);

% Coefficients (Daubechies-4)
h(1) = 0.482962913145;
h(2) = 0.836516303738;
h(3) = 0.224143868042;
h(4) = -0.129409522551;

% h(4) = 0.482962913145;
% h(3) = 0.836516303738;
% h(2) = 0.224143868042;
% h(1) = -0.129409522551;


% Normalize
for i = 1:L
    h(i) = h(i) / sqrt(2.0);
end

% Compute high-pass filter coefficients
for i = 1:L
    g(i) = (-1)^(i-1) * h(L - i + 1);
end

end