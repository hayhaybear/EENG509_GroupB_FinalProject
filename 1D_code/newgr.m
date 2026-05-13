function [xo, fo,No] = newgr(xi, fi, L, Nd, iw, xo, fo, No)
arguments (Input)
    xi, fi,L, Nd,iw;
end

arguments (Output)
    xo, fo, No;
end


Nmax = 260;
Lmax = 8;
Ndmax = 8;

% Guessing syntax
% Line 3
xi = zeros(Nmax,1);
fi = zeros()

% Lines 6 - 8
data = fi; % Copy fi into data


for idecomp = 1:Nd
Ndim = N/ (2^(idecomp-1));
% Call context
% Call filter
end

igrid