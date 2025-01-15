clear;
clc;

% Constants and parameters
a = 1e-4;                   % Crack length, m                                         
theta = 0.625;                 % Crack Geometry, b/a
b = a * theta;                   % Mineralization length, m
n = 1000;                    % Number of spatial nodes
T = 3;                      % Total time
M = 1000;                   % Number of time intervals
B = 0.7081374184;                    % Constant B in D_ij
sigma_0 = 2e7;             % far-field stress at depth of 500m
p_c = 1.89625e8;          % p_c, crystallization pressure
dt = T / M;                 % Time step size
sample_time_pick = 21;     % selected time number
% theta_values = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]; % Adjust number of theta values as needed
theta_values = [0.375, 0.625, 0.875]; % Adjust number of theta values as needed

fontSize_label = 20;
fontSize_title = 22;
fontSize_axis = 18;
fontName = 'Calibri';

% Variables
x = cos(pi * (1:n-1) / n);           % x_i values
s = cos(pi * (2*(1:n) - 1) / (2*n)); % s_j values
t = linspace(0, T, M+1);             % t_m values

% Preallocate arrays for results
F = zeros(n, M+1, length(theta_values)); % Initialize F(s_j, t_m)
delta = zeros(n-1, M+1, length(theta_values)); % Initialize delta(x_i, t_m)
sigma_yy = zeros(n-1, M+1, length(theta_values)); % Initialize sigma_yy(x_i, t_m)
kappa = zeros(M+1, length(theta_values)); % Initialize kappa(x_1, t_m)

F_far = zeros(n, M+1, length(theta_values)); % Initialize F(s_j, t_m)
delta_far = zeros(n-1, M+1, length(theta_values)); % Initialize delta(x_i, t_m)
sigma_yy_far = zeros(n-1, M+1, length(theta_values)); % Initialize sigma_yy(x_i, t_m)
kappa_far = zeros(M+1, length(theta_values)); % Initialize kappa(x_1, t_m)

% Compute P_ij for all i, j
P = zeros(n-1, n);
for i = 1:n-1
    for j = 1:n
        P(i, j) = (pi - acos(x(i))) / n + 2 / n * sum((sin(pi * (1:n-1)) - sin((1:n-1) * acos(x(i)))) ...
                  .* cos((1:n-1) * acos(s(j))) ./ (1:n-1));
    end
end

% Compute C_ij and D_ij for all i, j
C = zeros(n-1, n);
D = zeros(n-1, n);
for i = 1:n-1
    for j = 1:n
        C(i, j) = 1 / ((x(i) - s(j)) * n);
        D(i, j) = (B * M / T) * P(i, j);
    end
end

% Loop over different theta values
for k = 1:length(theta_values)
    theta = theta_values(k);

    % Solving the system iteratively for each time step
    for m = 2:M+1
        A = zeros(n, n);  % System matrix
        b = zeros(n, 1);  % Right-hand side vector

        A_far = zeros(n, n);  % System matrix
        b_far = zeros(n, 1);  % Right-hand side vector
        
        % Fill A and b based on the conditions
        for i = 1:n-1
            if abs(x(i)) < theta
                A(i, :) = (C(i, :) + D(i, :));
                b(i) = 2 + sum(D(i, :) .* F(:, m-1, k).');  % F(s_j, t_{m-1})

                A_far(i, :) = (C(i, :) + D(i, :));
                b_far(i) = 2 -(2*sigma_0/p_c) + sum(D(i, :) .* F_far(:, m-1, k).');  % F(s_j, t_{m-1})
            else
                A(i, :) = C(i, :);

                A_far(i, :) = C(i, :);
                b_far(i) =  - (2*sigma_0/p_c);
            end
        end
        A(n, :) = 1;  % Sum constraint
        A_far(n, :) = 1;  % Sum constraint
        
        % Solve the linear system A * f = b
        f = A \ b;
        F(:, m, k) = f;

        f_far = A_far \ b_far;
        F_far(:, m, k) = f_far;

        % Compute delta for each x_i and t_m
        delta(:, m, k) = P * F(:, m, k);

        delta_far(:, m, k) = P * F_far(:, m, k);

        % Compute sigma_yy for each x_i and t_m
        for i = 1:n-1
            sigma_yy(i, m, k) = (1/2) * 1/n * sum(F(:, m, k).' ./ (x(i) - s));
            sigma_yy_far(i, m, k) = (1/2) * 1/n * sum(F_far(:, m, k).' ./ (x(i) - s)) + (sigma_0/p_c);
        end
        
        % Compute kappa for x_1 and t_m
        kappa(m, k) = (sqrt(2)/4) * delta(1, m, k) / sqrt(1 - x(1));
        kappa_far(m, k) = (sqrt(2)/4) * delta_far(1, m, k) / sqrt(1 - x(1));
    end
end

% Plotting kappa against delta at x = 0 for different theta values
figure;
hold on;

% Generate colors based on the number of theta values
color_list = {'#0072BD', '#7E2F8E', '#A2142F'};  % Colors for each theta value

for k = 1:length(theta_values)
    color = color_list{k};  % Select color for the current line
    if theta_values(k) == 0.625
        plot(delta(n/2, :, k), kappa(:, k), 'Color', color, 'LineWidth', 1.5, 'DisplayName', ['\theta = ' num2str(theta_values(k))]);
        plot(delta_far(n/2, :, k), kappa_far(:, k),'--', 'Color', color, 'LineWidth', 1.5, 'DisplayName', ['\theta = ' num2str(theta_values(k))]);
    else
        plot(delta(n/2, :, k), kappa(:, k), 'Color', color, 'LineWidth', 1, 'DisplayName', ['\theta = ' num2str(theta_values(k))]); % black lines for other theta values
        plot(delta_far(n/2, :, k), kappa_far(:, k),'--', 'Color', color, 'LineWidth', 1, 'DisplayName', ['\theta = ' num2str(theta_values(k))]); % black lines for other theta values
    end    
end

hold off;
xlabel('H(0, T)', 'FontSize', fontSize_label,'FontName', fontName);
ylabel('$\hat{K}$', 'Interpreter', 'latex', 'FontSize', fontSize_label, 'FontWeight', 'bold','FontName', fontName);
legend('show', 'Location', 'eastoutside'); % Position the legend outside the plot
set(gca, 'FontSize', fontSize_axis, 'FontName', fontName);
set(gca,'box','on');
% xlim([0 1.01]); % Set vertical axis from 0 to 1
ylim([0 0.45]); % Set vertical axis from 0 to 1
