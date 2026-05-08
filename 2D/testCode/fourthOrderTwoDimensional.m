% Original in paper but too large for 8GB of Macbook
L_x = 20000000;
L_y = 8000000;
dx = 156000;
dy = 95000;
dt = 3050; 

% Good enough I think
L_x = 3000*255;
L_y = 3000*255;
dx = 3000;
dy = 3000; 
dt = 600; 

numSteps = 5000;

g_prime = 0.049;  
lambda = 1 * 10^3;

x = -L_x/2:dx:L_x/2;
y = -L_y/2:dy:L_y/2;
[X, Y] = meshgrid(x,y);

% This allows for Corioulis to change as Y changes stopping that problem
% from earlier
beta = 2e-11;
f_0 = 0; 
f = f_0 + beta .* Y;

H = 40;
delta_h = 60;
l_x = 20000; 
l_y = 20000; 

h = H + delta_h * exp(-(X.^2/l_x^2 + Y.^2/l_y^2));
hu = zeros(size(h));
hv = zeros(size(h));

% This is from Google Gemini 
figure('Color', 'w');
hPlot = surf(X/1e3, Y/1e3, h);
shading interp; lighting gouraud; camlight;
xlabel('x (km)'); ylabel('y (km)'); zlabel('Height (m)');
colorbar;
zlim([30 80]); 
clim([30 80]);

timerVal = tic; 
for t = 1:numSteps
% I used Google Gemini to refine and find errors in my original adaptation
% of the RK4 finite difference method

    % k1
    [dhdt_1, dhudt_1, dhvdt_1] = derivativeStuff(h, hu, hv, dx, dy, g_prime, f, lambda);
    [h, hu, hv] = apply_boundaries_optimized(h, hu, hv); % Base boundaries
    k1_h  = dt .* dhdt_1;
    k1_hu = dt .* dhudt_1;
    k1_hv = dt .* dhvdt_1; 
    
    % k2
    h_mid1  = h  + 0.5 .* k1_h;
    hu_mid1 = hu + 0.5 .* k1_hu;
    hv_mid1 = hv + 0.5 .* k1_hv;
    [h_mid1, hu_mid1, hv_mid1] = apply_boundaries_optimized(h_mid1, hu_mid1, hv_mid1); % Mid1 boundaries
    [dhdt_2, dhudt_2, dhvdt_2] = derivativeStuff(h_mid1, hu_mid1, hv_mid1, dx, dy, g_prime, f, lambda);
    k2_h  = dt .* dhdt_2;
    k2_hu = dt .* dhudt_2;
    k2_hv = dt .* dhvdt_2;
    
    % k3
    h_mid2  = h  + 0.5 .* k2_h;
    hu_mid2 = hu + 0.5 .* k2_hu;
    hv_mid2 = hv + 0.5 .* k2_hv;
    [h_mid2, hu_mid2, hv_mid2] = apply_boundaries_optimized(h_mid2, hu_mid2, hv_mid2); % Mid2 boundaries
    [dhdt_3, dhudt_3, dhvdt_3] = derivativeStuff(h_mid2, hu_mid2, hv_mid2, dx, dy, g_prime, f, lambda);
    k3_h  = dt .* dhdt_3;
    k3_hu = dt .* dhudt_3;
    k3_hv = dt .* dhvdt_3;
    
    % k4
    h_end  = h  + k3_h;
    hu_end = hu + k3_hu;
    hv_end = hv + k3_hv;
    [h_end, hu_end, hv_end] = apply_boundaries_optimized(h_end, hu_end, hv_end); % End boundaries
    [dhdt_4, dhudt_4, dhvdt_4] = derivativeStuff(h_end, hu_end, hv_end, dx, dy, g_prime, f, lambda);
    k4_h  = dt .* dhdt_4;
    k4_hu = dt .* dhudt_4;
    k4_hv = dt .* dhvdt_4;
    
    % combine everything
    h  = h  + (1/6) .* (k1_h  + 2.*k2_h  + 2.*k3_h  + k4_h);
    hu = hu + (1/6) .* (k1_hu + 2.*k2_hu + 2.*k3_hu + k4_hu);
    hv = hv + (1/6) .* (k1_hv + 2.*k2_hv + 2.*k3_hv + k4_hv);
    
    % Final boundary application for the timestep
    [h, hu, hv] = apply_boundaries_optimized(h, hu, hv);
    
    threshold = 0.01; 
    
    hit_top    = any(abs(h(3, :) - H) > threshold);
    hit_bottom = any(abs(h(end-2, :) - H) > threshold);
    hit_left   = any(abs(h(:, 3) - H) > threshold);
    hit_right  = any(abs(h(:, end-2) - H) > threshold);
    
    if hit_top || hit_bottom || hit_left || hit_right
        fprintf(['Boundary hit detected at timestep %d at time ' ...
            '%d. Ending program.\n'], t, toc(timerVal));
        break; 
    end
    
    if mod(t, 2) == 0
        set(hPlot, 'ZData', h);
        title(['Time Step: ', num2str(t)]);
        drawnow limitrate;
    end

    if t == 40
        clim([30 50]);
    end
end

baselineH = h; 

% derivativeStuff is from Google Gemini but refined from 1D
function [dhdt, dhudt, dhvdt] = derivativeStuff(h, hu, hv, dx, dy, g_prime, f, lambda)
    safe = max(h, 0.1); 
    u = hu ./ safe; 
    v = hv ./ safe; 
    
    dh_dx = SpatialX(h, dx);
    dh_dy = SpatialY(h, dy);
    
    dhu_dx = SpatialX(hu, dx);
    dhv_dy = SpatialY(hv, dy);
    
    du2h_dx = SpatialX(hu .* u, dx);
    duvh_dy = SpatialY(hu .* v, dy);
    
    dv2h_dy = SpatialY(hv .* v, dy);
    duvh_dx = SpatialX(hv .* u, dx);
    
    laplacian_hu = Laplacian(hu, dx, dy);
    laplacian_hv = Laplacian(hv, dx, dy);
    
    dhdt  = -dhu_dx - dhv_dy; 
    
    dhudt = -(g_prime .* h .* dh_dx) - du2h_dx - duvh_dy + (f .* hv) + (lambda .* laplacian_hu); 
    
    dhvdt = -(g_prime .* h .* dh_dy) - dv2h_dy - duvh_dx - (f .* hu) + (lambda .* laplacian_hv); 
end


function dfdx = SpatialX(F, dx)
    dfdx = zeros(size(F)); 
    dfdx(:, 3:end-2) = (-F(:, 5:end) + 8.*F(:, 4:end-1) - 8.*F(:, 2:end-3) + F(:, 1:end-4)) ./ (12*dx);
end

function dfdy = SpatialY(F, dy)
    dfdy = zeros(size(F)); 
    dfdy(3:end-2, :) = (-F(5:end, :) + 8.*F(4:end-1, :) - 8.*F(2:end-3, :) + F(1:end-4, :)) ./ (12*dy);
end

function laplaceResult = Laplacian(F, dx, dy)
    xPart = zeros(size(F)); 
    yPart = zeros(size(F)); 
    
    xPart(:, 3:end-2) = (-F(:, 5:end) + 16.*F(:, 4:end-1) - 30.*F(:, 3:end-2) + 16.*F(:, 2:end-3) - F(:, 1:end-4)) ./ (12*dx^2);
    
    yPart(3:end-2, :) = (-F(5:end, :) + 16.*F(4:end-1, :) - 30.*F(3:end-2, :) + 16.*F(2:end-3, :) - F(1:end-4, :)) ./ (12*dy^2);
    
    laplaceResult = xPart + yPart; 
end

function [h, hu, hv] = apply_boundaries_optimized(h, hu, hv)
    h(:, 1:2) = repmat(h(:, 3), 1, 2);
    h(:, end-1:end) = repmat(h(:, end-2), 1, 2);
    hu(:, 1:2) = 0; hu(:, end-1:end) = 0;
    hv(:, 1:2) = 0; hv(:, end-1:end) = 0;
    
    h(1:2, :) = repmat(h(3, :), 2, 1);
    h(end-1:end, :) = repmat(h(end-2, :), 2, 1);
    hu(1:2, :) = 0; hu(end-1:end, :) = 0;
    hv(1:2, :) = 0; hv(end-1:end, :) = 0;
end






