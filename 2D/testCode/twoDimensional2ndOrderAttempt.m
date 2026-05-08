%% WOFD-AHO Shallow Water Model - Case 1
clear; clc;

% --- Parameters Case One (from table) ---
L_x = 20000000;      % 20,000 km -> meters
L_y = 8000000;       % 8,000 km -> meters
dx = 156000;         % 156 km -> meters
dy = 95000;          % 95 km -> meters
dt = 3050;           % Time step (s)
g_prime = 0.049;     % Reduced gravity
beta = 2e-11;        % Beta parameter
lambda = 1e4;        % Horizontal viscosity (A)
f0 = 1e-4;           % Reference Coriolis frequency

% Domain setup
H = 40;              % Background height
delta_h = 60;        % Bump height
l_x = 667000;        % Initial scale x
l_y = 334000;        % Initial scale y
numSteps = 1000;

% --- Mesh Grid ---
x = -L_x/2:dx:L_x/2;
y = -L_y/2:dy:L_y/2;
[X, Y] = meshgrid(x, y);

% Define varying Coriolis parameter (Beta-plane: f = f0 + beta*y)
F_coriolis = f0 + beta .* Y;

% --- Initial Conditions ---
u = zeros(size(X)); 
v = zeros(size(X)); 
h = H + delta_h*(exp(-((X.^2) / l_x^2) - (Y.^2) / l_y^2));

% Plot the initial h (Your original style)
figure(1);
surf(X, Y, h);
shading interp;
title('Initial Condition Test');
xlabel('X (m)'); ylabel('Y (m)');

%% Stabilized WOFD-AHO Model
% ... (Keep your parameters and Initial Conditions the same) ...

for t = 1:numSteps
    % 1. Create a "smoothed" version of old values to prevent checkerboarding
    % This is a simple Lax-style filter
    h_smooth = h;
    h_smooth(2:end-1, 2:end-1) = 0.25 * (h(1:end-2, 2:end-1) + h(3:end, 2:end-1) + ...
                                         h(2:end-1, 1:end-2) + h(2:end-1, 3:end));
    
    hu_old = h .* u; 
    hv_old = h .* v;

    % 2. Calculate Right Hand Sides (using the smoothed height for pressure)
    % We use 'h_smooth' in the pressure gradient to stop the noise
    rhs_u = (F_coriolis .* hv_old) ...
            - (g_prime .* h .* SpatialX(h, dx)) ... 
            + (lambda .* Laplacian(hu_old, dx, dy)) ...
            - (SpatialX(u.*hu_old, dx) + SpatialY(v.*hu_old, dy));

    rhs_v = -(F_coriolis .* hu_old) ...
            - (g_prime .* h .* SpatialY(h, dy)) ...
            + (lambda .* Laplacian(hv_old, dx, dy)) ...
            - (SpatialX(u.*hv_old, dx) + SpatialY(v.*hv_old, dy));

    rhs_h = -(SpatialX(hu_old, dx) + SpatialY(hv_old, dy));

    % 3. Update using the "Lax" averaged values
    % Replacing h_old with h_smooth here is key for stability
    hu = hu_old + dt .* rhs_u;
    hv = hv_old + dt .* rhs_v;
    h  = h_smooth + dt .* rhs_h; % The smoothing happens here

    % 4. Prevent division by zero and update velocity
    h(h < 1) = 1; % Physical floor so depth doesn't go negative
    u = hu ./ h;
    v = hv ./ h;

    % 5. Boundary Conditions (Hard walls)
    u([1 end], :) = 0; u(:, [1 end]) = 0;
    v([1 end], :) = 0; v(:, [1 end]) = 0;
    % No-gradient height
    h(1,:) = h(2,:); h(end,:) = h(end-1,:);
    h(:,1) = h(:,2); h(:,end) = h(:,end-1);

    % Visualization
    if mod(t, 20) == 0
        surf(X, Y, h); shading interp; zlim([30 110]);
        title(['Step: ', num2str(t)]);
        drawnow;
    end
end

% --- Finite Difference Functions ---
function dfdx = SpatialX(F, dx)
    dfdx = zeros(size(F)); 
    dfdx(:, 2:end-1) = (F(:, 3:end) - F(:, 1:end-2)) / (2 * dx);
end

function dfdy = SpatialY(F, dy)
    dfdy = zeros(size(F)); 
    dfdy(2:end-1, :) = (F(3:end, :) - F(1:end-2, :)) / (2 * dy);
end

function laplaceResult = Laplacian(F, dx, dy)
    xPart = zeros(size(F)); 
    yPart = zeros(size(F)); 
    xPart(:, 2:end-1) = (F(:, 3:end) - 2*F(:, 2:end-1) + F(:, 1:end-2)) / (dx^2);
    yPart(2:end-1, :) = (F(3:end, :) - 2*F(2:end-1, :) + F(1:end-2, :)) / (dy^2);
    laplaceResult = xPart + yPart;
end