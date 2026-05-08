%% Test run

% Parameters Case One 
L_x = 20000000;
dx = 156000;
dt = 3050; 
numSteps = 5000;
g_prime = 0.049;
beta = 2 * 10^-11;
lambda = 1 * 10^4; 
f = 7.2921*10^-5; 
l_x = 667000; 
H = 40;
delta_h = 60;

% Plot Points Grid
x = -L_x/2:dx:L_x/2;

% Initial Conditions
u = zeros(size(x)); 
h = H + delta_h*(exp(-(x.^2) / l_x^2));
plot(x,h); 

% Try Finite Difference
for t = 1:numSteps
    hu = u .* h; 
    old_h = h; 

    % RK4 was implemented from https://www.cfm.brown.edu/people/dobrush/am33/Matlab/ch3/RK4.html
    % and debugging from chatgpt

    [dhdt_1, dhudt_1] = derivativeStuff(h, hu, dx, g_prime, lambda);
    k1_h  = dt .* dhdt_1;
    k1_hu = dt .* dhudt_1;
    
    h_mid1  = h  + 0.5 .* k1_h;
    hu_mid1 = hu + 0.5 .* k1_hu;

    [dhdt_2, dhudt_2] = derivativeStuff(h_mid1, hu_mid1, dx, g_prime, lambda);
    k2_h  = dt .* dhdt_2;
    k2_hu = dt .* dhudt_2;
    
    h_mid2  = h  + 0.5 .* k2_h;
    hu_mid2 = hu + 0.5 .* k2_hu;

    [dhdt_3, dhudt_3] = derivativeStuff(h_mid2, hu_mid2, dx, g_prime, lambda);
    k3_h  = dt .* dhdt_3;
    k3_hu = dt .* dhudt_3;
    
    h_end  = h  + k3_h;
    hu_end = hu + k3_hu;
    
    [dhdt_4, dhudt_4] = derivativeStuff(h_end, hu_end, dx, g_prime, lambda);
    k4_h  = dt .* dhdt_4;
    k4_hu = dt .* dhudt_4;
    
    h  = h  + (1/6) .* (k1_h  + 2.*k2_h  + 2.*k3_h  + k4_h);
    hu = hu + (1/6) .* (k1_hu + 2.*k2_hu + 2.*k3_hu + k4_hu);
    
    u = hu ./ h; 
    
    h(1:2) = h(3); 
    h(end-1:end) = h(end-2);
    u(1:2) = 0; 
    u(end-1:end) = 0;
    hu = h .* u; 

    plot(x,h);
    ylim([20,100])
    pause(0.00005)
    disp(t)
end


function [dhdt, dhudt] = derivativeStuff(h, hu, dx, g_prime, lambda)
    u = hu ./ h; 
    
    dh_dx = SpatialX(h, dx);
    du2h_dx = SpatialX(u.*u.*h, dx);
    dhu_dx = SpatialX(hu, dx);
    laplacian_uh = Laplacian(u.*h, dx);
    
    dhdt  = -dhu_dx; 
    dhudt = -(g_prime .* h .* dh_dx) + (lambda .* laplacian_uh)- du2h_dx; 
end

% Fourth order derivative functions
function dfdx = SpatialX(F, dx)
    dfdx = zeros(size(F)); 
    dfdx(3:end-2) = (-F(5:end) + 8.*F(4:end-1) - 8.*F(2:end-3) + F(1:end-4)) ./ (12*dx);
end

function laplaceResult = Laplacian(F, dx)
    xPart = zeros(size(F)); 
    xPart(3:end-2) = (-F(5:end)+16.*F(4:end-1)-30.*F(3:end-2)+16.*F(2:end-3)-F(1:end-4))./(12*dx^2);
    laplaceResult = xPart; 
end