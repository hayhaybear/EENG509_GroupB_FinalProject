% Google Gemini was used for most of the plotting logic and it was
% prompted to make functions out of most the logic from the original code
% for readability purposes

% Grid Parameters
L_x = 3000*255;
L_y = 3000*255;
dx = 3000;
dy = 3000; 
dt = 600/3; 
numSteps = 5000;
remeshNum = 5; 

% Parameters
g_prime = 0.049;  
lambda = 1 * 10^3;
x = -L_x/2:dx:L_x/2;
y = -L_y/2:dy:L_y/2;
[X, Y] = meshgrid(x,y);
beta = 2e-11;
beta = 0; 
f_0 = 0; 
f = f_0 + beta .* Y;

% Initial Conditions
H = 40;
delta_h = 60;
l_x = 20000; 
l_y = 20000; 
h = H + delta_h * exp(-(X.^2/l_x^2 + Y.^2/l_y^2));
hu = zeros(size(h));
hv = zeros(size(h));
figure('Color', 'w');
filter = zeros(size(h));
activeIDX = []; 

% These plots was from Google Gemini
subplot(1, 2, 1);
hPlot = surf(X/1e3, Y/1e3, h);
shading interp; lighting gouraud; camlight;
title_h = title(['Water Height, Time: 0']);
zlim([30 80]); clim([30 80]);

subplot(1, 2, 2);
maskPlot = imagesc(x/1e3, y/1e3, filter);
axis xy; 
title('Wavelet Refinement Mask');
colormap(gca, 'parula');

timerVal = tic; 
for t = 1:numSteps
    if mod(t, remeshNum) == 1 || t == 1
        [~,~,grid_mask] = newGrid2d_v2(0, h, 0.5, 4, 0); 
        dialateMask = imdilate(grid_mask, ones(5,5));
    
        safeZone = false(size(h));
        safeZone(3:end-2, 3:end-2) = 1; 
        filter = dialateMask & safeZone;
    
        activeIDX = find(filter);
    end
    
    [h, hu, hv] = rk4_step(h, hu, hv, dt, dx, dy, g_prime, f, lambda, activeIDX);

    % when it hits the wall, just end simulation
    threshold = 0.01; 
    
    hit_top    = any(abs(h(3, :) - H) > threshold);
    hit_bottom = any(abs(h(end-2, :) - H) > threshold);
    hit_left   = any(abs(h(:, 3) - H) > threshold);
    hit_right  = any(abs(h(:, end-2) - H) > threshold);
    
    if hit_top || hit_bottom || hit_left || hit_right
        fprintf(['Boundary hit detected at timestep %d at time ' ...
            '%d.\n'], t, toc(timerVal));
        break; 
    end
    
    % This plot is from Google Gemini
    if mod(t, 2) == 0
        set(hPlot, 'ZData', h);
        
        rgb_mask = repmat(1 - filter, 1, 1, 3); 
        
        set(maskPlot, 'CData', rgb_mask);
        
        set(title_h, 'String', ['Water Height, Time: ', num2str(t)]);
        drawnow limitrate;
    end
end

% CITE THIS IN PAPER WHEN YOU WRITE
% http://www.cfm.brown.edu/people/dobrush/am33/Matlab/ch3/rk4.html

function [h, hu, hv] = rk4_step(h, hu, hv, dt, dx, dy, g_prime, f, lambda, active_indexes)
    % k1
    [dhdt_1, dhudt_1, dhvdt_1] = derivativeStuff(h, hu, hv, dx, dy, g_prime, f, lambda, active_indexes);
    [h, hu, hv] = apply_boundaries(h, hu, hv);
    k1_h  = dt .* dhdt_1;
    k1_hu = dt .* dhudt_1;
    k1_hv = dt .* dhvdt_1; 
    
    % k2
    h2 = h + 0.5.*k1_h; hu2 = hu + 0.5.*k1_hu; hv2 = hv + 0.5.*k1_hv;
    [h2, hu2, hv2] = apply_boundaries(h2, hu2, hv2); % Apply to h2, not h
    [dhdt_2, dhudt_2, dhvdt_2] = derivativeStuff(h2, hu2, hv2, dx, dy, g_prime, f, lambda, active_indexes);
    k2_h  = dt .* dhdt_2;
    k2_hu = dt .* dhudt_2;
    k2_hv = dt .* dhvdt_2;
    
    % k3
    h3 = h + 0.5.*k2_h; hu3 = hu + 0.5.*k2_hu; hv3 = hv + 0.5.*k2_hv;
    [h3, hu3, hv3] = apply_boundaries(h3, hu3, hv3); % Apply to h3, not h
    [dhdt_3, dhudt_3, dhvdt_3] = derivativeStuff(h3, hu3, hv3, dx, dy, g_prime, f, lambda, active_indexes);
    k3_h  = dt .* dhdt_3;
    k3_hu = dt .* dhudt_3;
    k3_hv = dt .* dhvdt_3;
    
    % k4
    h4 = h + k3_h; hu4 = hu + k3_hu; hv4 = hv + k3_hv;
    [h4, hu4, hv4] = apply_boundaries(h4, hu4, hv4); % Apply to h4, not h
    [dhdt_4, dhudt_4, dhvdt_4] = derivativeStuff(h4, hu4, hv4, dx, dy, g_prime, f, lambda, active_indexes);
    k4_h  = dt .* dhdt_4;
    k4_hu = dt .* dhudt_4;
    k4_hv = dt .* dhvdt_4;
    
    % Combine
    h  = h  + (1/6) .* (k1_h  + 2.*k2_h  + 2.*k3_h  + k4_h);
    hu = hu + (1/6) .* (k1_hu + 2.*k2_hu + 2.*k3_hu + k4_hu);
    hv = hv + (1/6) .* (k1_hv + 2.*k2_hv + 2.*k3_hv + k4_hv);
    
    % Final boundary application
    [h, hu, hv] = apply_boundaries(h, hu, hv);
end

% Extracted your exact boundary logic into a helper function
function [h, hu, hv] = apply_boundaries(h, hu, hv)
    h(:, 1:2) = repmat(h(:, 3), 1, 2);
    h(:, end-1:end) = repmat(h(:, end-2), 1, 2);
    hu(:, 1:2) = 0; hu(:, end-1:end) = 0;
    hv(:, 1:2) = 0; hv(:, end-1:end) = 0;
    
    h(1:2, :) = repmat(h(3, :), 2, 1);
    h(end-1:end, :) = repmat(h(end-2, :), 2, 1);
    hu(1:2, :) = 0; hu(end-1:end, :) = 0;
    hv(1:2, :) = 0; hv(end-1:end, :) = 0;
end

function [dhdt, dhudt, dhvdt] = derivativeStuff(h, hu, hv, dx, dy, g_prime, f, lambda, active_indexes)
    safe = max(h, 0.1); 
    u = hu ./ safe; 
    v = hv ./ safe; 
    
    dh_dx = SpatialX(h, dx, active_indexes);
    dh_dy = SpatialY(h, dy, active_indexes);
    dhu_dx = SpatialX(hu, dx, active_indexes);
    dhv_dy = SpatialY(hv, dy, active_indexes);
    du2h_dx = SpatialX(hu .* u, dx, active_indexes);
    duvh_dy = SpatialY(hu .* v, dy, active_indexes);
    dv2h_dy = SpatialY(hv .* v, dy, active_indexes);
    duvh_dx = SpatialX(hv .* u, dx, active_indexes);
    
    laplacian_hu = Laplacian(hu, dx, dy, active_indexes);
    laplacian_hv = Laplacian(hv, dx, dy, active_indexes);
    
    dhdt  = -dhu_dx - dhv_dy; 
    dhudt = -(g_prime .* h .* dh_dx) - du2h_dx - duvh_dy + (f .* hv) + (lambda .* laplacian_hu); 
    dhvdt = -(g_prime .* h .* dh_dy) - dv2h_dy - duvh_dx - (f .* hu) + (lambda .* laplacian_hv); 
end

function dfdx = SpatialX(F, dx, active_indexes)
    dfdx = zeros(size(F)); 
    N = size(F, 1); 
  
    dfdx(active_indexes) = (-F(active_indexes + 2*N) + ...
                            8.*F(active_indexes + N) - ...
                            8.*F(active_indexes - N) + ...
                            F(active_indexes - 2*N)) ./ (12*dx);
end

function dfdy = SpatialY(F, dy, active_indexes)
    dfdy = zeros(size(F)); 
    dfdy(active_indexes) = (-F(active_indexes+2) + ...
                           8.*F(active_indexes+1) - ...
                           8.*F(active_indexes-1) + ...
                           F(active_indexes-2)) ./ (12*dy); 
end

function laplaceResult = Laplacian(F, dx, dy, active_indexes)
    xPart = zeros(size(F)); 
    yPart = zeros(size(F)); 
    N = size(F, 1);
    xPart(active_indexes) = (-F(active_indexes + 2*N) + ...
                            16.*F(active_indexes + N) - ...
                            30.*F(active_indexes) + ...
                            16.*F(active_indexes - N) - ...
                            F(active_indexes - 2*N)) ./ (12*dx^2);

    yPart(active_indexes) = (-F(active_indexes+2) + ...
                            16.*F(active_indexes+1) - ...
                            30.*F(active_indexes) + ...
                            16.*F(active_indexes-1) - ...
                            F(active_indexes-2)) ./ (12*dy^2); 

    laplaceResult = xPart + yPart; 
end