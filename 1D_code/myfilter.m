function [low, high] = myfilter(Extdata, h, g, N, L)
% Equivalent MATLAB version of the Fortran subroutine filter
% (renamed to avoid conflict with MATLAB built-in "filter")

% disp("Extdata length " + length(Extdata));
Nmax = 260;
Lmax = 8;

% Preallocate
low  = zeros(Nmax/2 + Lmax, 1);
high = zeros(Nmax/2 + Lmax, 1);

% Initialize (not strictly needed in MATLAB, but kept for clarity)
for i = 1:(Nmax/2 + Lmax)
    low(i)  = 0.0;
    high(i) = 0.0;
end

% Main filtering loop

for i = 1:(N/2 + (L-2)/2)
    for j = 1:L
        
        % Adjusted indexing for MATLAB (1-based)
        ij = 2*(i-1) + j - (L-2);
        
        idx = ij + 2;  % matches Extdata(ij+2) in Fortran
        
        % Boundary safety (important in MATLAB)
        if idx >= 1 && idx <= length(Extdata)
            low(i)  = low(i)  + h(j) * Extdata(idx);
            high(i) = high(i) + g(j) * Extdata(idx);
        end
        
    end
end

% [c, s] = wavedec(Extdata, 1, 'db2');
% length(low)
% s
% lowMat = c(1:length(c)/2);
% highMat = c(length(c)/2 + 1:end);

% low = lowMat;
% high = highMat;

% approx = appcoef(c,s,"db2")
% 
% s
% figure; plot(lowMat);
% hold on; plot(low,'--'); 
% plot(approx,':'); hold off;
% title("Low Coefficients");
% legend("Matlab", "Fortran")
% 
% figure; plot(highMat);
% hold on; plot(high, '--'); hold off;
% title("High Coeficients")
% legend("Matlab", "Fortran")




end