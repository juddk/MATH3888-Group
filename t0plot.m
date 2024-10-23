clear;
% Parameters
sigma = 0.2;            % sigmaatility
r = 0.1;                % Interest rate
K = 100;                 % K price
T = 1;                  % T time
NAS = 1000;             % Number of x steps
Le_values = [0, 0.2, 0.4, 0.6]; % Array of different Le values
rapm_values = [0, 0.2, 0.4, 0.6]; % Array of different rapm values

% Call the function to plot the results
RAPM_multiple(sigma, r, K, T, NAS, rapm_values);
% Call the function to plot the results
L_multiple(sigma, r, K, T, NAS, Le_values);

% Function definition (placed at the bottom)
function L_multiple(sigma, r, K, T, NAS, Le_values)
    % Constants
    beta = 2 * r / sigma^2;
    tildeT = T * sigma^2 / 2;
    
    % Discretization
    x_min = -7;
    x_max = 1.5;
    x = linspace(x_min, x_max, NAS + 1); % Stock prices
    dx = (x_max - x_min) / NAS;
    NTauS = 15 * NAS; % Number of time steps
    dtau = tildeT / NTauS;
    tau = linspace(0, sigma^2 / 2, NTauS + 1);
    U = zeros(NAS + 1, NTauS + 1); % Option value array
    
    % Initialization
    U(:, 1) = max(1 - exp(-x), 0); % Vectorized initialization
    
    % Precompute constant coefficients
    dx2_inv = 1 / (dx^2); % Inverse of dx^2 to avoid recalculating in loop

    % Prepare the figure
    figure;
    hold on; % Keep the plot for multiple Le values

    % Loop over different Le values
    for i = 1:length(Le_values)
        Le = Le_values(i);
        
        % Reset U for each Le value
        U(:, 1) = max(1 - exp(-x), 0);
        
        % Time stepping (numerical solution)
        for k = 2:(NTauS + 1) 
            % Calculate U_x and U_xx
            U_x = (U(3:end, k - 1) - U(1:end-2, k - 1)) / (2 * dx); 
            U_xx = (U(3:end, k - 1) - 2 * U(2:end-1, k - 1) + U(1:end-2, k - 1)) * dx2_inv;
            UxUxx = U_x + U_xx;
            
            % Update option values (vectorized for inner loop)
            U(2:end-1, k) = dtau * ((1 + Le * abs(UxUxx)) .* UxUxx + beta * U_x) + U(2:end-1, k - 1);
            
            % Boundary condition at S=infinity
            U(end, k) = 1 - exp(-beta * (k - 1) * dtau - x(end));  
        end

        % Convert log prices to stock prices
        S = K * exp(x); 
        t = T - 2 * tau / sigma^2;

        % Compute V_num (numerical solution) in vectorized form
        V_num = S' .* U;

        % Limit to the range S in [0, 2K]
        idx = (S <= 1.2 * K)&(S >= 0.5*K); 
        S = S(idx);     
        V_num = V_num(idx, :);

        % Extract the values for t = 0
        [~, t_zero_idx] = min(abs(t)); % Find index closest to t = 0
        V_t0 = V_num(:, t_zero_idx);    % V values at t = 0
        
        % Plot V against S at t = 0 for this Le value
        plot(S, V_t0, 'LineWidth', 2, 'DisplayName', ['Le = ' num2str(Le)]);
    end
    
    % Add labels and legend
    xlabel('Stock Price (S)');
    ylabel('Option Value (V) at t = 0');
    legend show; % Show legend with Le values
    grid on;
    hold off;
    print('C:\Users\14192\Desktop\group project\diffLe', '-depsc');
end


% Function definition
function RAPM_multiple(sigma, r, K, T, NAS, rapm_values)
    % Constants
    beta = 2 * r / sigma^2;
    tildeT = T * sigma^2 / 2;
    
    % Discretization
    x_min = -7;
    x_max = 1.5;
    x = linspace(x_min, x_max, NAS + 1); % Stock prices
    dx = (x_max - x_min) / NAS;
    NTauS = 15 * NAS; % Number of time steps
    dtau = tildeT / NTauS;
    tau = linspace(0, sigma^2 / 2, NTauS + 1);
    U = zeros(NAS + 1, NTauS + 1); % Option value array
    
    % Prepare the figure
    figure;
    hold on; % Keep the plot for multiple rapm values

    % Loop over different rapm values
    for i = 1:length(rapm_values)
        rapm = rapm_values(i);
            % Precompute constant coefficients
    dx2_inv = 1 / (dx^2); % Inverse of dx^2 to avoid recalculating in loop
        % Reset U for each rapm value
        U(:, 1) = max(1 - exp(-x), 0); % Initial condition

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
        idx = (S <= 1.2 * K) & (S >= 0.5 * K); 
        S = S(idx);     
        V_num = V_num(idx, :);

        % Extract the values for t = 0
        [~, t_zero_idx] = min(abs(t)); % Find index closest to t = 0
        V_t0 = V_num(:, t_zero_idx);    % V values at t = 0
        
        % Adjust the line thickness dynamically based on rapm values
        line_thickness = 2 + i * 0.5;  % Example: thicker lines as i increases
        
        % Plot V against S at t = 0 for this rapm value
        plot(S, V_t0, 'LineWidth', 2, 'DisplayName', ['rapm = ' num2str(rapm)]);
    end
    
    % Add labels and legend
    xlabel('Stock Price (S)');
    ylabel('Option Value (V) at t = 0');
    legend show; % Show legend with rapm values
    grid on;
    hold off;
    print('C:\Users\14192\Desktop\group project\diffRapm', '-depsc');
end
