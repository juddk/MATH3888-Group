clear;
% Parameters
sigma = 0.6;            % sigmaatility
r = 0.2;                % Interest rate
K = 10;                 % K price
T = 1;                  % T time
NAS = 1000;             % Number of x steps
Le = 0.7;               % leland
rapm = 0.6;               % rapm

BS(sigma, r, K, T, NAS);
L(sigma, r, K, T, NAS, Le);
RAPM(sigma, r, K, T, NAS, rapm);

function BS(sigma, r, K, T, NAS)
    % Constants
    beta = 2 * r / sigma^2;
    tildeT = T * sigma^2 / 2;
    
    % Discretization
    x_min = -7;
    x_max = 1.5;
    x = linspace(x_min, x_max, NAS + 1); % Stock prices
    dx = (x_max - x_min) / NAS;
    NTauS = 15000; % Number of time steps
    dtau = tildeT / NTauS;
    tau = linspace(0, sigma^2*T/ 2, NTauS + 1);
    U = zeros(NAS + 1, NTauS + 1); % Option value array
    
    % Initialization
    U(:, 1) = max(1 - exp(-x), 0); % Vectorized initialization
    
    % Precompute constant coefficients
    dx2_inv = 1 / (dx^2); % Inverse of dx^2 to avoid recalculating in loop

    % Time stepping (numerical solution)
    for k = 2:(NTauS + 1)
        % Calculate U_x and U_xx
        U_x = (U(3:end, k - 1) - U(1:end-2, k - 1)) / (2 * dx); 
        U_xx = (U(3:end, k - 1) - 2 * U(2:end-1, k - 1) + U(1:end-2, k - 1)) * dx2_inv;
        
        % Update option values (vectorized for inner loop)
        U(2:end-1, k) = dtau * ((1 + beta) * U_x + U_xx) + U(2:end-1, k - 1);
        
        % Boundary condition at S=infinity
        U(end, k) = 1 - exp(-beta * (k - 1) * dtau - x(end));  
    end

    % Convert log prices to stock prices
    S = K * exp(x); 
    t = T - 2 * tau / sigma^2;
    % Compute V_num (numerical solution) in vectorized form
    V_num = S' .* U;

    % Limit to the range S in [0, 2K]
    idx = S <= 2 * K; 
    S = S(idx);     
    V_num = V_num(idx, :);

    % Plot the surface of the numerical solution
figure;
surf(S, t, V_num');
shading interp;
colorbar;
% Add labels for the X, Y, and Z axes
xlabel('Stock Price (S)');   % X-axis label
ylabel('Time (t)');          % Y-axis label
zlabel('Option Value (V)');  % Z-axis label

% Save the figure
print('C:\Users\14192\Desktop\group project\V_B', '-depsc');

    % Calculate analytic solution for the same stock prices and times
    V_analytic = zeros(size(V_num)); % Preallocate for analytic solution

    % Vectorized calculation for analytic solution
    for j = 1:length(t)
        if t(j) == T
            V_analytic(:, j) = max(S - K, 0); % Payoff at T (call)
        else
            d1 = (log(S / K) + (r + 0.5 * sigma^2) * (T - t(j))) ./ (sigma * sqrt(T - t(j)));
            d2 = d1 - sigma * sqrt(T - t(j));
            V_analytic(:, j) = S .* normcdf(d1) - K * exp(-r * (T - t(j))) .* normcdf(d2);
        end
    end
    % Compute and plot the absolute difference between numerical and analytic solutions
    figure;
    surf(S, t, abs(V_analytic - V_num)');
    shading interp;
    colorbar;
    view(2); 

    % Add labels for the X, Y, and Z axes
    xlabel('Stock Price (S)');   % X-axis label
    ylabel('Time (t)');          % Y-axis label
    print('C:\Users\14192\Desktop\group project\V_B-V', '-depsc');
end



function L(sigma, r, K, T, NAS, Le)
    % Constants
    beta = 2 * r / sigma^2;
    tildeT = T * sigma^2 / 2;
    
    % Discretization
    x_min = -7;
    x_max = 1.5;
    x = linspace(x_min, x_max, NAS + 1); % Stock prices
    dx = (x_max - x_min) / NAS;
    NTauS = 15*NAS; % Number of time steps
    dtau = tildeT / NTauS;
    tau = linspace(0, sigma^2 / 2, NTauS + 1);
    U = zeros(NAS + 1, NTauS + 1); % Option value array
    
    % Initialization
    U(:, 1) = max(1 - exp(-x), 0); % Vectorized initialization
    
    % Precompute constant coefficients
    dx2_inv = 1 / (dx^2); % Inverse of dx^2 to avoid recalculating in loop

    % Time stepping (numerical solution)
for k = 2:(NTauS + 1) 
    % Calculate U_x and U_xx
    U_x = (U(3:end, k - 1) - U(1:end-2, k - 1)) / (2 * dx); 
    U_xx = (U(3:end, k - 1) - 2 * U(2:end-1, k - 1) + U(1:end-2, k - 1)) * dx2_inv;
    UxUxx = U_x + U_xx;
    % Update option values (vectorized for inner loop)
    U(2:end-1, k) = dtau * ((1 + Le * sign(UxUxx)) .* UxUxx + beta * U_x) + U(2:end-1, k - 1);
    
    % Boundary condition at S=infinity
    U(end, k) = 1 - exp(-beta * (k - 1) * dtau - x(end));  
end


    % Convert log prices to stock prices
    S = K * exp(x); 
    t = T - 2 * tau / sigma^2;

    % Compute V_num (numerical solution) in vectorized form
    V_num = S' .* U;

    % Limit to the range S in [0, 2K]
    idx = S <= 2 * K; 
    S = S(idx);     
    V_num = V_num(idx, :);

    % Plot the surface of the numerical solution
    figure;
    surf(S, t, V_num');
    shading interp;
    % Add labels for the X, Y, and Z axes
    xlabel('Stock Price (S)');   % X-axis label
    ylabel('Time (t)');          % Y-axis label
    zlabel('Option Value (V)');  % Z-axis label
    colorbar;
    print('C:\Users\14192\Desktop\group project\V_L', '-depsc');
    % Calculate analytic solution for the same stock prices and times
    V_analytic = zeros(size(V_num)); % Preallocate for analytic solution

    % Vectorized calculation for analytic solution
    for j = 1:length(t)
        if t(j) == T
            V_analytic(:, j) = max(S - K, 0); % Payoff at T (call)
        else
            d1 = (log(S / K) + (r + 0.5 * sigma^2) * (T - t(j))) ./ (sigma * sqrt(T - t(j)));
            d2 = d1 - sigma * sqrt(T - t(j));
            V_analytic(:, j) = S .* normcdf(d1) - K * exp(-r * (T - t(j))) .* normcdf(d2);
        end
    end
    % Compute and plot the absolute difference between numerical and analytic solutions
    figure;
    surf(S, t, (V_num-V_analytic)');
    colorbar;
    view(2);
    % Add labels for the X, Y, and Z axes
xlabel('Stock Price (S)');   % X-axis label
ylabel('Time (t)');          % Y-axis label
zlabel('Option Value (V)');  % Z-axis label
    shading interp;
    print('C:\Users\14192\Desktop\group project\V_L-V', '-depsc');
end


function RAPM(sigma, r, K, T, NAS, rapm)
    % Constants
    beta = 2 * r / sigma^2;
    tildeT = T * sigma^2 / 2;
    
    % Discretization
    x_min = -7;
    x_max = 1.5;
    x = linspace(x_min, x_max, NAS + 1); % Stock prices
    dx = (x_max - x_min) / NAS;
    NTauS = 15*NAS; % Number of time steps
    dtau = tildeT / NTauS;
    tau = linspace(0, sigma^2 / 2, NTauS + 1);
    U = zeros(NAS + 1, NTauS + 1); % Option value array
    
    % Initialization
    U(:, 1) = max(1 - exp(-x), 0); % Vectorized initialization
    
    % Precompute constant coefficients
    dx2_inv = 1 / (dx^2); % Inverse of dx^2 to avoid recalculating in loop

    % Time stepping (numerical solution)
for k = 2:(NTauS + 1) 
    % Calculate U_x and U_xx
    U_x = (U(3:end, k - 1) - U(1:end-2, k - 1)) / (2 * dx); 
    U_xx = (U(3:end, k - 1) - 2 * U(2:end-1, k - 1) + U(1:end-2, k - 1)) * dx2_inv;
    UxUxx = U_x + U_xx;
    RAPM_term = abs(UxUxx).^(1/3);
    % Update option values (vectorized for inner loop)
    U(2:end-1, k) = dtau * ((1 + rapm * RAPM_term) .* UxUxx + beta * U_x) + U(2:end-1, k - 1);
    
    % Boundary condition at S=infinity
    U(end, k) = 1 - exp(-beta * (k - 1) * dtau - x(end));  
end


    % Convert log prices to stock prices
    S = K * exp(x); 
    t = T - 2 * tau / sigma^2;

    % Compute V_num (numerical solution) in vectorized form
    V_num = S' .* U;

    % Limit to the range S in [0, 2K]
    idx = S <= 2 * K; 
     S = S(idx);     
     V_num = V_num(idx, :);

    % Plot the surface of the numerical solution
    figure;
    surf(S, t, V_num');
    shading interp;
    colorbar;
    % Add labels for the X, Y, and Z axes
xlabel('Stock Price (S)');   % X-axis label
ylabel('Time (t)');          % Y-axis label
zlabel('Option Value (V)');  % Z-axis label
    print('C:\Users\14192\Desktop\group project\V_R', '-depsc');

    % Calculate analytic solution for the same stock prices and times
    V_analytic = zeros(size(V_num)); % Preallocate for analytic solution

    % Vectorized calculation for analytic solution
    for j = 1:length(t)
        if t(j) == T
            V_analytic(:, j) = max(S - K, 0); % Payoff at T (call)
        else
            d1 = (log(S / K) + (r + 0.5 * sigma^2) * (T - t(j))) ./ (sigma * sqrt(T - t(j)));
            d2 = d1 - sigma * sqrt(T - t(j));
            V_analytic(:, j) = S .* normcdf(d1) - K * exp(-r * (T - t(j))) .* normcdf(d2);
        end
    end
    % Compute and plot the absolute difference between numerical and analytic solutions
    figure;
    surf(S, t, (V_num-V_analytic)');
    shading interp;
    colorbar;
    view(2);
    % Add labels for the X, Y, and Z axes
    xlabel('Stock Price (S)');   % X-axis label
    ylabel('Time (t)');          % Y-axis label

    print('C:\Users\14192\Desktop\group project\V_R-V', '-depsc');
end

