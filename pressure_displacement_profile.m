clear;
clc;

% MATLAB code to solve the provided linear system of equations and plot the results

% Constants and parameters
% a = 1e-4;                   % half Crack length, m                                         
theta = 0.625;                 % Crack Geometry, b/a
% b = a * theta;                   % Mineralization length, m
n = 1000;                    % Number of spatial nodes
T = 3;                      % Total time
M = 1000;                   % Number of time intervals
B = 0.7081374184;                    % Constant B in D_ij
% B = 2.8243;
dt = T / M;                 % Time step size
sample_time_pick = 25;     % selected time number

fontSize_label = 16;
fontSize_title = 22;
fontSize_axis = 16;
fontName = 'Calibri';  % Define the desired font

% Variables
x = cos(pi * (1:n-1) / n);           % x_i values
s = cos(pi * (2*(1:n) - 1) / (2*n)); % s_j values
t = linspace(0, T, M+1);             % t_m values
F = zeros(n, M+1);                   % Initialize F(s_j, t_m)
delta = zeros(n-1, M+1);             % Initialize delta(x_i, t_m)
sigma_yy = zeros(n-1, M+1);          % Initialize sigma_yy(x_i, t_m)
kappa = zeros(1, M+1);               % Initialize kappa(x_1, t_m)

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

% Solving the system iteratively for each time step
for m = 2:M+1
    A = zeros(n, n);  % System matrix
    b = zeros(n, 1);  % Right-hand side vector
    
    % Fill A and b based on the conditions
    for i = 1:n-1
        if abs(x(i)) < theta
            A(i, :) = (C(i, :) + D(i, :));
            b(i) = 2 + sum(D(i, :) .* F(:, m-1).');  % F(s_j, t_{m-1})
        else
            A(i, :) = C(i, :);
        end
    end
    A(n, :) = 1;  % Sum constraint
    
    % Solve the linear system A * f = b
    f = A \ b;
    F(:, m) = f;

    % Compute delta for each x_i and t_m
    delta(:, m) = P * F(:, m);

    % Compute sigma_yy for each x_i and t_m
    for i = 1:n-1
        sigma_yy(i, m) = (1/2) * 1/n * sum(F(:, m).' ./ (x(i) - s));
    end
    
    % Compute kappa for x_1 and t_m
    kappa(m) = (sqrt(2)/4) * (delta(1, m) / sqrt(1 - x(1)));
end

% Plotting at ten specific time steps
figure;
plotIndices = round(linspace(1, M+1, sample_time_pick));  % Indices of time steps to plot

colors = jet(length(plotIndices));  % Generate a colormap using 'jet' for rainbow colors

subplot(2, 1, 1);
hold on;
for i = 1:length(plotIndices)
    plot(x, delta(:, plotIndices(i)),'Color', colors(i,:),'LineWidth', 1.5);
end
hold off;
xlabel('X', 'FontSize', fontSize_label,'FontName', fontName);
ylabel('H(X,T)', 'FontSize', fontSize_label,'FontName', fontName);
% title('Non-dimensional displacement along the crack, \Delta(X, T)', 'FontSize', fontSize_title, 'FontWeight', 'bold');
colormap(jet); % Set colormap
colorbar('Ticks',linspace(0,1,5),'TickLabels', arrayfun(@(idx) sprintf('T=%.3f', t(plotIndices(idx))), linspace(1, sample_time_pick, 5), 'UniformOutput', false));
set(gca, 'FontSize', fontSize_axis,'FontName', fontName);
set(gca,'box','on');
% ylim([0,2])


subplot(2, 1, 2);
hold on;
for i = 1:length(plotIndices)
    plot(x, sigma_yy(:, plotIndices(i)),'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlabel('X', 'FontSize', fontSize_label,'FontName', fontName);
ylabel('P(X, T)', 'FontSize', fontSize_label,'FontName', fontName);
% title('Non-dimensional stress along the crack, P(X, T)', 'FontSize', fontSize_title, 'FontWeight', 'bold');
colormap(jet); % Set colormap
colorbar('Ticks',linspace(0,1,5),'TickLabels', arrayfun(@(idx) sprintf('T=%.3f', t(plotIndices(idx))), linspace(1, sample_time_pick, 5), 'UniformOutput', false));
set(gca, 'FontSize', fontSize_axis,'FontName', fontName);
set(gca,'box','on');



% subplot(3, 1, 3);
% hold on;
% plot(t(plotIndices),kappa(plotIndices),  '-', 'LineWidth', 1.5);
% hold off;
% xlabel('T', 'FontSize', fontSize_label, 'FontWeight', 'bold');
% ylabel('$\hat{K}$', 'Interpreter', 'latex', 'FontSize', fontSize_label, 'FontWeight', 'bold');
% title('Time to reach the stress intensity factor at crack tip in non-dimensional', 'FontSize', fontSize_title, 'FontWeight', 'bold');
% set(gca, 'FontSize', fontSize_axis);
% set(gca,'box','on');

