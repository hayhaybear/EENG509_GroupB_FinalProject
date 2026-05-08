%% Test run

% Parameters Case One 
L_x = 20000000;
L_y = 8000000; 
dx = 156000;
dy = 95000; 
dt = 3050; 
numSteps = 1000;
g_prime = 0.049;
beta = 2 * 10^-11; 
lambda = 1 * 10^-4;
f = 7.2921*10^-5; 
l_x = 667000; 
l_y = 334000; 
H = 40;
delta_h = 60;

% Mesh Grid
x = -L_x/2:dx:L_x/2;
y = -L_y/2:dy:L_y/2;
[X, Y] = meshgrid(x,y);
length_x = length(x);
length_y = length(y); 

% Initial Conditions
u = zeros(size(X)); 
v = zeros(size(X)); 

h = H + delta_h*(exp(-((X.^2) / l_x^2) - (Y.^2) / l_y^2));

% Plot the intial h
surf(X,Y, h)
shading interp
title('Initial Condition Test')

% Try Finite Difference
for t = 1:numSteps
    hu = h.*u;
    hv = h.*v; 

    h_old = h; 
    u_old = u; 
    v_old = v; 
 
    hu_old = hu; 
    hv_old = hv;

    rightSide19 = (f.*h.*v) - (g_prime.*h.*SpatialX(h,dx)) + (lambda.*Laplacian(u.*h,dx,dy)) ...
                  -(SpatialX(u.*u.*h, dx)) - (SpatialY(u.*v.*h,dy)); 
    hu = hu_old + dt.*rightSide19;

    rightSide20 = -(f.*h.*u) - (g_prime.*h.*SpatialY(h,dy)) + (lambda.*Laplacian(v.*h,dx,dy))...
                  - (SpatialX((u.*v.*h), dx)) - (SpatialY((v.*v.*h),dy));
    hv = hv_old + dt.*rightSide20; 

    dhdt = -(SpatialX(hu_old, dx)) - (SpatialY(hv_old, dy));

    h = h_old + dt.*dhdt;
    h(h < 0.1) = 0.1; 
    u = hu ./ h; 
    v = hv ./ h;

    % Boundary Conditions
    u(:, 1) = 0;  
    u(:, end) = 0;  
    h(:, 1) = h(:, 2);       
    h(:, end) = h(:, end-1);   
    v(:, 1) = v(:, 2);
    v(:, end) = v(:, end-1);
    
    v(1, :) = 0;  
    v(end, :) = 0;  
    h(1, :) = h(2, :);       
    h(end, :) = h(end-1, :);   
    u(1, :) = u(2, :);
    u(end, :) = u(end-1, :);

    if (mod(t,10) == 0)
        surf(X,Y, h);
        zlim([20 120]);
        pause(0.1); 
    end

end


% Second order derivative functions
function dfdx = SpatialX(F, dx)
    dfdx = zeros(size(F)); 
    dfdx(:,2:end-1) = (F(:,3:end) - F(:,1:end-2)) / (2*dx);
end

function dfdy = SpatialY(F, dy)
    dfdy = zeros(size(F)); 
    dfdy(2:end-1,:) = (F(3:end,:) - F(1:end-2,:)) / (2*dy);
end

function laplaceResult = Laplacian(F, dx, dy)
    xPart = zeros(size(F)); 
    yPart = zeros(size(F)); 

    xPart(:,2:end-1) = (F(:,3:end) - 2*F(:,2:end-1) + F(:,1:end-2)) / (dx^2);
    yPart(2:end-1,:) = (F(3:end,:) - 2*F(2:end-1,:) + F(1:end-2,:)) / (dy^2);

    laplaceResult = xPart + yPart;
end