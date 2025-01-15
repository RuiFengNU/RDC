clc;
clear;

% Define the constant
thetas = [3/8, 5/8, 7/8];

fontSize_label = 20;
fontSize_title = 22;
fontSize_axis = 18;
fontName = 'Calibri';

% Define the range for xi
xi = linspace(0, 1, 1000); % Avoiding 0 to prevent division by zero

% Initialize figure
figure;

% Plot for each theta
hold on;

color_list = {'#0072BD', '#7E2F8E', '#A2142F'};  % Colors for each theta value

for k = 1:length(thetas)
    theta = thetas(k);
    color = color_list{k};  % Select color for the current line
    % Define the function hatK_crit
    hatK_crit = sqrt(1 ./ xi) .* (2 / pi) .* asin(theta * xi);
    % Plot the function
    plot(xi, hatK_crit,'Color', color, 'LineWidth', 1, 'DisplayName', ['\theta = ', num2str(theta)]);
end

% Add labels and legend
xlabel('$a_0/a$', 'Interpreter', 'latex', 'FontSize', fontSize_label, 'FontWeight', 'bold','FontName', fontName);
ylabel('$\hat{K}_{crit}$', 'Interpreter', 'latex', 'FontSize', fontSize_label, 'FontWeight', 'bold','FontName', fontName);
legend show;
set(gca, 'FontSize', fontSize_axis,'FontName', fontName);
set(gca,'box','on');
% grid on;

% Display the plot
hold off;
shg;
