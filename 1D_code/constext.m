function Extdata = constext(data, N, L)
% Equivalent MATLAB version of the Fortran subroutine constext
% This code extends the data before giving it to the filter command
% This ensures that as the filter slides over the function values there is
% enough room at the end to calculated the filtering
Nmax = 260;
Lmax = 8;

% Preallocate
Extdata = zeros(Nmax + Lmax, 1);

halfL = L/2;

% Copy main data into extended array
for i = 1:N
    idx = halfL + i - 1;
    Extdata(idx) = data(i);
end

% Left extension (constant using first value)
for i = 1:(L-3)
    idx = halfL - i;
    if idx >= 1
        Extdata(idx) = data(1);
    end
end

% Right extension (constant using last value)
for i = 1:(L-1)
    idx = N + halfL + i - 1;
    if idx <= length(Extdata)
        Extdata(idx) = data(N);
    end
end

end