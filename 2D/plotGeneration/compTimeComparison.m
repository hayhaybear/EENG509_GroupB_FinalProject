% parameters for building graphs
timeStepsBeforeReMesh = [2 5 10 15 20 25];
epsilon = [0.05 0.1 0.2 0.5 1 2];
numSimulations = length(epsilon);
numSteps = 5000;
bigRMSE = zeros(1, numSimulations); 
bigTimeElapsed = zeros(1, numSimulations);
bigTimeElapsed_optimized = zeros(1, numSimulations);
guassian = 1;
pyramid = 0; 

for i = 1:numSimulations
    % Grid Parameters
    L_x = 3000*255;
    L_y = 3000*255;
    dx = 3000;
    dy = 3000; 
    dt = 600; 
    
    % Parameters
    g_prime = 0.049;  
    lambda = 1 * 10^3;
    x = -L_x/2:dx:L_x/2;
    y = -L_y/2:dy:L_y/2;
    [X, Y] = meshgrid(x,y);
    beta = 0; 
    f_0 = 0; 
    f = f_0 + beta .* Y;
    
    % Initial Conditions
    H = 40;
    delta_h = 60;
    l_x = 20000; 
    l_y = 20000; 
    
    if guassian == 1
        h = H + delta_h * exp(-(X.^2/l_x^2 + Y.^2/l_y^2));
    end

    if guassian == 0
        r2 = X.^2 + Y.^2;
        h = H + delta_h * (1 - r2/l_x^2) .* exp(-r2/l_x^2);
    end

    if pyramid == 1
        depression_depth = min(delta_h, H - 0.1);
        h = H - depression_depth * max(0, 1 - abs(X)/l_x) .* max(0, 1 - abs(Y)/l_y);
    end

    hu = zeros(size(h));
    hv = zeros(size(h));
    filtered_mask = zeros(size(h));
    active_indexes = []; 
    
    h_optimized = h;
    hu_optimized = hu;
    hv_optimized = hv; 
    
    timeNotOptimized = tic;
    for t = 1:numSteps
        % k1
        [dhdt_1, dhudt_1, dhvdt_1] = derivativeStuff(h, hu, hv, dx, dy, g_prime, f, lambda);
        [h, hu, hv] = apply_boundaries_optimized(h, hu, hv); 
        k1_h  = dt .* dhdt_1;
        k1_hu = dt .* dhudt_1;
        k1_hv = dt .* dhvdt_1; 
        
        % k2
        h_mid1  = h  + 0.5 .* k1_h;
        hu_mid1 = hu + 0.5 .* k1_hu;
        hv_mid1 = hv + 0.5 .* k1_hv;
        [h_mid1, hu_mid1, hv_mid1] = apply_boundaries_optimized(h_mid1, hu_mid1, hv_mid1); 
        [dhdt_2, dhudt_2, dhvdt_2] = derivativeStuff(h_mid1, hu_mid1, hv_mid1, dx, dy, g_prime, f, lambda);
        k2_h  = dt .* dhdt_2;
        k2_hu = dt .* dhudt_2;
        k2_hv = dt .* dhvdt_2;
        
        % k3
        h_mid2  = h  + 0.5 .* k2_h;
        hu_mid2 = hu + 0.5 .* k2_hu;
        hv_mid2 = hv + 0.5 .* k2_hv;
        [h_mid2, hu_mid2, hv_mid2] = apply_boundaries_optimized(h_mid2, hu_mid2, hv_mid2); 
        [dhdt_3, dhudt_3, dhvdt_3] = derivativeStuff(h_mid2, hu_mid2, hv_mid2, dx, dy, g_prime, f, lambda);
        k3_h  = dt .* dhdt_3;
        k3_hu = dt .* dhudt_3;
        k3_hv = dt .* dhvdt_3;
        
        % k4
        h_end  = h  + k3_h;
        hu_end = hu + k3_hu;
        hv_end = hv + k3_hv;
        [h_end, hu_end, hv_end] = apply_boundaries_optimized(h_end, hu_end, hv_end); 
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
            fprintf('Boundary hit detected at timestep %d (Remesh interval: %d, Case: %d)\n', t, timeStepsBeforeReMesh(i), i);
            break; 
        end
    end
    bigTimeElapsed(i) = toc(timeNotOptimized)

    timeOptimized = tic; 
    for t = 1:numSteps
        if mod(t, timeStepsBeforeReMesh(i)) == 1 || t == 1
            [~,~,grid_mask] = newGrid2d_v2(0, h_optimized, 0.5, 4, 0); 
            dialated_mask = imdilate(grid_mask, ones(5,5));
        
            safe_zone = false(size(h));
            safe_zone(3:end-2, 3:end-2) = true; 
            filtered_mask = dialated_mask & safe_zone;
        
            active_indexes = find(filtered_mask);
        end
    
        [h_optimized, hu_optimized, hv_optimized] = rk4_step_optimized(h_optimized, hu_optimized, hv_optimized, ...
                                                    dt, dx, dy, g_prime, f, lambda, active_indexes);
    
        threshold = 0.01; 
        
        hit_top    = any(abs(h_optimized(3, :) - H) > threshold);
        hit_bottom = any(abs(h_optimized(end-2, :) - H) > threshold);
        hit_left   = any(abs(h_optimized(:, 3) - H) > threshold);
        hit_right  = any(abs(h_optimized(:, end-2) - H) > threshold);
        
        if hit_top || hit_bottom || hit_left || hit_right
            fprintf('Boundary hit detected at timestep %d (Remesh interval: %d, Case: %d)\n', t, timeStepsBeforeReMesh(i), i);
            break; 
        end
    end 
    bigTimeElapsed_optimized(i) = toc(timeOptimized)

    bigRMSE(i) = sqrt(mean((h(:) - h_optimized(:)).^2)); 
end

% Plotting code was refined from Google Gemini
figure;

% Plot 1: RMSE
subplot(2,1,1);
plot(timeStepsBeforeReMesh, bigRMSE, '-o', 'LineWidth', 2);
title('RMSE vs. Re-Mesh Interval');
ylabel('RMSE');
grid on;

% Plot 2: Time Comparison
subplot(2,1,2);
plot(timeStepsBeforeReMesh, bigTimeElapsed, '-s', 'DisplayName', 'Standard Time');
hold on;
plot(timeStepsBeforeReMesh, bigTimeElapsed_optimized, '--^', 'DisplayName', 'Optimized Time');
title('Computation Time vs. Re-Mesh Interval');
xlabel('Time Steps Before Re-Mesh');
ylabel('Time (s)');
legend('Location', 'best');
grid on;

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

function [h, hu, hv] = rk4_step_optimized(h, hu, hv, dt, dx, dy, g_prime, f, lambda, active_indexes)
    % k1
    [dhdt_1, dhudt_1, dhvdt_1] = derivativeStuff_optimized(h, hu, hv, dx, dy, g_prime, f, lambda, active_indexes);
    [h, hu, hv] = apply_boundaries_optimized(h, hu, hv);
    k1_h  = dt .* dhdt_1;
    k1_hu = dt .* dhudt_1;
    k1_hv = dt .* dhvdt_1; 
    
    % k2
    h2 = h + 0.5.*k1_h; hu2 = hu + 0.5.*k1_hu; hv2 = hv + 0.5.*k1_hv;
    [h2, hu2, hv2] = apply_boundaries_optimized(h2, hu2, hv2); % Apply to h2, not h
    [dhdt_2, dhudt_2, dhvdt_2] = derivativeStuff_optimized(h2, hu2, hv2, dx, dy, g_prime, f, lambda, active_indexes);
    k2_h  = dt .* dhdt_2;
    k2_hu = dt .* dhudt_2;
    k2_hv = dt .* dhvdt_2;
    
    % k3
    h3 = h + 0.5.*k2_h; hu3 = hu + 0.5.*k2_hu; hv3 = hv + 0.5.*k2_hv;
    [h3, hu3, hv3] = apply_boundaries_optimized(h3, hu3, hv3); % Apply to h3, not h
    [dhdt_3, dhudt_3, dhvdt_3] = derivativeStuff_optimized(h3, hu3, hv3, dx, dy, g_prime, f, lambda, active_indexes);
    k3_h  = dt .* dhdt_3;
    k3_hu = dt .* dhudt_3;
    k3_hv = dt .* dhvdt_3;
    
    % k4
    h4 = h + k3_h; hu4 = hu + k3_hu; hv4 = hv + k3_hv;
    [h4, hu4, hv4] = apply_boundaries_optimized(h4, hu4, hv4); % Apply to h4, not h
    [dhdt_4, dhudt_4, dhvdt_4] = derivativeStuff_optimized(h4, hu4, hv4, dx, dy, g_prime, f, lambda, active_indexes);
    k4_h  = dt .* dhdt_4;
    k4_hu = dt .* dhudt_4;
    k4_hv = dt .* dhvdt_4;
    
    % Combine
    h  = h  + (1/6) .* (k1_h  + 2.*k2_h  + 2.*k3_h  + k4_h);
    hu = hu + (1/6) .* (k1_hu + 2.*k2_hu + 2.*k3_hu + k4_hu);
    hv = hv + (1/6) .* (k1_hv + 2.*k2_hv + 2.*k3_hv + k4_hv);
    
    % Final boundary application
    [h, hu, hv] = apply_boundaries_optimized(h, hu, hv);
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

function [dhdt, dhudt, dhvdt] = derivativeStuff_optimized(h, hu, hv, dx, dy, g_prime, f, lambda, active_indexes)
    safe = max(h, 0.1); 
    u = hu ./ safe; 
    v = hv ./ safe; 
    
    dh_dx = SpatialX_optimized(h, dx, active_indexes);
    dh_dy = SpatialY_optimized(h, dy, active_indexes);
    dhu_dx = SpatialX_optimized(hu, dx, active_indexes);
    dhv_dy = SpatialY_optimized(hv, dy, active_indexes);
    du2h_dx = SpatialX_optimized(hu .* u, dx, active_indexes);
    duvh_dy = SpatialY_optimized(hu .* v, dy, active_indexes);
    dv2h_dy = SpatialY_optimized(hv .* v, dy, active_indexes);
    duvh_dx = SpatialX_optimized(hv .* u, dx, active_indexes);
    
    laplacian_hu = Laplacian_optimized(hu, dx, dy, active_indexes);
    laplacian_hv = Laplacian_optimized(hv, dx, dy, active_indexes);
    
    dhdt  = -dhu_dx - dhv_dy; 
    dhudt = -(g_prime .* h .* dh_dx) - du2h_dx - duvh_dy + (f .* hv) + (lambda .* laplacian_hu); 
    dhvdt = -(g_prime .* h .* dh_dy) - dv2h_dy - duvh_dx - (f .* hu) + (lambda .* laplacian_hv); 
end

function dfdx = SpatialX_optimized(F, dx, active_indexes)
    dfdx = zeros(size(F)); 
    N = size(F, 1); 
  
    dfdx(active_indexes) = (-F(active_indexes + 2*N) + ...
                            8.*F(active_indexes + N) - ...
                            8.*F(active_indexes - N) + ...
                            F(active_indexes - 2*N)) ./ (12*dx);
end

function dfdy = SpatialY_optimized(F, dy, active_indexes)
    dfdy = zeros(size(F)); 
    dfdy(active_indexes) = (-F(active_indexes+2) + ...
                           8.*F(active_indexes+1) - ...
                           8.*F(active_indexes-1) + ...
                           F(active_indexes-2)) ./ (12*dy); 
end

function laplaceResult = Laplacian_optimized(F, dx, dy, active_indexes)
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