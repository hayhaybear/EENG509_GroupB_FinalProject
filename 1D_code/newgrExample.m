function [xo, fo, No] = newgrExample(xi, fi, L, N, th, Nd, iw)
% Equivalent MATLAB version of the Fortran subroutine newgr

Nmax = 260;
Lmax = 8;
Ndmax = 8;

% Preallocate
xo = zeros(Nmax,1);
fo = zeros(Nmax,1);
h  = zeros(Lmax,1);
g  = zeros(Lmax,1);

HPF = zeros(Nmax/2 + Lmax, Ndmax);
data = zeros(Nmax,1);
Extdata = zeros(Nmax + Lmax,1);

% Get filter coefficients
[h, g] = getcoef(L); %--------------

% Copy data
data = fi;

% Decomposition loop
% Gets the coarsest scaling coefficients and stores them in data
% Stores the wavelet coefficieints in HPF
for idecomp = 1:Nd
    
    Ndim = N / (2^(idecomp-1));
    
    Extdata = constext(data, Ndim, L);
    
    [data, HPF(:, idecomp)] = myfilter(Extdata, h, g, Ndim, L);
end
% disp("length data: " + length(data))




% Grid selection
igrid = 0;

% Repeat for all points
for ipnt = 1:N
    
    % Base grid condition
    % Makes sure that the grid at the user defined coarsest scale is kept
    if abs(mod(ipnt-1, 2^Nd)) < 1e-5
        igrid = igrid + 1;
        xo(igrid) = xi(ipnt);
        fo(igrid) = fi(ipnt);
    end
    
    % Multilevel refinement
    % Repeat for all levels of decomposition
    for idecomp = 1:Nd
        
        n1 = abs(ipnt - 2^(idecomp-1) - 1);
        n2 = 2^(idecomp);
        
        % Make sure that we are looking at the correct scale
        if abs(mod(n1, n2)) < 1e-5
            
            indexl = 1 + round(n1 / n2);
            iflagpoint = 0;
            
            for iwiden = -iw:iw
                iindex = indexl + iwiden;
                
                % Boundary check
                if iindex <= 1 || iindex >= N / (2^(idecomp))
                    iindex = indexl;
                end
                
                if abs(HPF(iindex, idecomp)) > th
                    iflagpoint = iflagpoint + 1;
                end
            end
            
            if iflagpoint >= (iw + 1)
                igrid = igrid + 1;
                xo(igrid) = xi(ipnt);
                fo(igrid) = fi(ipnt);
            end
            
        end
    end
end

% Add final point
xo(igrid+1) = xi(N+1);
fo(igrid+1) = fi(N+1);

No = igrid + 1;

end