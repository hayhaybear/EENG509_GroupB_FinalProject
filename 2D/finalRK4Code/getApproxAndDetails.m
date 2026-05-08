function [A1,H1,V1,D1] = getApproxAndDetails(filterData,waveletFilter,decompNum,N)
% Function calculates one level wavelet transform and returns approximation
% (A1) and detail coefficients(H1,V1,D1)
arguments (Input)
    filterData
    waveletFilter
    decompNum % Used to determine upsampling
    N
end

arguments (Output)
    A1, H1,V1,D1
end
H1 = zeros(N,N); V1 = zeros(N,N); D1 = zeros(N,N);

[c,s] = wavedec2(filterData,1,waveletFilter);
[h1,v1,d1] = detcoef2('all', c,s,1);
A1 = appcoef2(c,s,waveletFilter,1);

H1(1:2^decompNum:end,1:2^decompNum:end) = h1;
V1(1:2^decompNum:end,1:2^decompNum:end) = v1;
D1(1:2^decompNum:end,1:2^decompNum:end) = d1;


% A1=1;H1=1;V1=1;D1 =1;
end