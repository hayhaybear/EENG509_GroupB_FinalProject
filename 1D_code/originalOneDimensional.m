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
    uh = u .* h; 
    old_h = h; 
    old_uh = uh;

    dhu_dt = -g_prime*h.*SpatialX(h,dx) - SpatialX(u.*u.*h, dx) + lambda*Laplacian(h.*u, dx);
    dh_dt = -SpatialX(h.*u,dx); 

    h = old_h + dt*dh_dt;
    uh = old_uh + dt*dhu_dt;
    u = uh ./ h;

    % boundary conditions
    h(1) = h(2);
    h(end) = h(end-1);

    u(1) = 0;
    u(end) = 0; 

    plot(x,h);
    ylim([20,100])
    pause(0.00005)
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