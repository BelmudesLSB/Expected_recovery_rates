%%% This file creates an Mex file using the source code from C++ %%%

%% Clear all
clc; 
clear;

%% mex 

mex main.cpp economy.cpp initialization.cpp auxiliary.cpp 

%% Common parameters:
params.b_grid_size = 50;                % Number of points in the grid for the bond price.
params.b_grid_min = -0.6;               % Minimum value of the bond price.
params.b_grid_max = 0.00;               % Maximum value of the bond price.
params.y_grid_size = 6;                 % Number of points in the grid for the income.
params.y_default = 0.969;               % Maxsimum income under default.
params.beta = 0.953;                    % Discount factor.
params.gamma = 2;                       % Risk aversion.
params.r = 0.017;                       % Interest rate.
params.rho = 0.945;                     % Persistence of the income.
params.sigma = 0.025;                   % Standard deviation of the income.
params.theta = 0.282;                   % Probability of a re-entry.
params.max_iter = 2000;                 % Maximum number of iterations.
params.tol = 1e-7;                      % Tolerance for the convergence.
params.m = 3;                           % Number of standard deviations for the income grid.

% Benchmark model:
params_benchmark.b_grid_size = params.b_grid_size;       % Number of points in the grid for the bond price.
params_benchmark.b_grid_min = params.b_grid_min;         % Minimum value of the bond price.
params_benchmark.b_grid_max = params.b_grid_max;         % Maximum value of the bond price.
params_benchmark.y_grid_size = params.y_grid_size;       % Number of points in the grid for the income.
params_benchmark.y_default = params.y_default;           % Maximum income under default.
params_benchmark.beta = params.beta;                     % Discount factor.
params_benchmark.gamma = params.gamma;                   % Risk aversion.
params_benchmark.r = params.r;                           % Interest rate.
params_benchmark.rho = params.rho;                       % Persistence of the income.
params_benchmark.sigma = params.sigma;                   % Standard deviation of the income.
params_benchmark.theta = params.theta;                   % Probability of a re-entry.
params_benchmark.alpha = 0;                              % Recovery on defaulted debt.
params_benchmark.max_iter = params.max_iter;             % Maximum number of iterations.
params_benchmark.tol = params.tol;                       % Tolerance for the convergence.
params_benchmark.m = params.m;                           % Number of standard deviations for the income grid.

% Alpha low calibration:
params_low.b_grid_size = params.b_grid_size;       % Number of points in the grid for the bond price.
params_low.b_grid_min = params.b_grid_min;         % Minimum value of the bond price.
params_low.b_grid_max = params.b_grid_max;         % Maximum value of the bond price.
params_low.y_grid_size = params.y_grid_size;       % Number of points in the grid for the income.
params_low.y_default = params.y_default;           % Maximum income under default.
params_low.beta = params.beta;                     % Discount factor.
params_low.gamma = params.gamma;                   % Risk aversion.
params_low.r = params.r;                           % Interest rate.
params_low.rho = params.rho;                       % Persistence of the income.
params_low.sigma = params.sigma;                   % Standard deviation of the income.
params_low.theta = params.theta;                   % Probability of a re-entry.
params_low.alpha = 0.15;                           % Recovery on defaulted debt.
params_low.max_iter = params.max_iter;             % Maximum number of iterations.
params_low.tol = params.tol;                       % Tolerance for the convergence.
params_low.m = params.m;                           % Number of standard deviations for the income grid.

% Alpha high calibration:
params_high.b_grid_size = params.b_grid_size;      % Number of points in the grid for the bond price.
params_high.b_grid_min = params.b_grid_min;        % Minimum value of the bond price.
params_high.b_grid_max = params.b_grid_max;        % Maximum value of the bond price.
params_high.y_grid_size = params.y_grid_size;      % Number of points in the grid for the income.
params_high.y_default = params.y_default;          % Maximum income under default.
params_high.beta = params.beta;                    % Discount factor.
params_high.gamma = params.gamma;                  % Risk aversion.
params_high.r = params.r;                          % Interest rate.
params_high.rho = params.rho;                      % Persistence of the income.s
params_high.sigma = params.sigma;                  % Standard deviation of the income.
params_high.theta = params.theta;                  % Probability of a re-entry.
params_high.alpha = 0.30;                          % Recovery on defaulted debt.
params_high.max_iter = params.max_iter;            % Maximum number of iterations.
params_high.tol = params.tol;                     % Tolerance for the convergence.
params_high.m = params.m;                         % Number of standard deviations for the income grid.

%% Run code with both alphas;

tic;
arellano_solution_alpha_benchmark = main(params_benchmark);
toc;
tic;
arellano_solution_alpha_low = main(params_low);
toc;
tic;
arellano_solution_alpha_high = main(params_high);
toc;
save('Solution_alpha_benchmark', 'arellano_solution_alpha_benchmark')
save('Solution_alpha_low', 'arellano_solution_alpha_low')
save('Solution_alpha_high', 'arellano_solution_alpha_high')
save('Parameters_benchmark.mat', 'params')
save('Parameters_low', 'params_low')
save('Parameters_high', 'params_high')

%% Format variables:

% Exogenous:
Y_grid = arellano_solution_alpha_low.Y_grid;
Y_grid_default = arellano_solution_alpha_low.Y_grid_default;
B_grid = arellano_solution_alpha_low.B_grid;
P = reshape(arellano_solution_alpha_low.P, params.y_grid_size, params.y_grid_size)';

% For benchamark Calibration:
Q_benchmark = reshape(arellano_solution_alpha_benchmark.Q, params.b_grid_size, params.y_grid_size)';
B_policy_benchmark = reshape(arellano_solution_alpha_benchmark.B_policy, params.b_grid_size, params.y_grid_size)'+1; %MATLAB starts indexing at 1.
D_policy_benchmark = reshape(arellano_solution_alpha_benchmark.D_policy, params.b_grid_size, params.y_grid_size)';
B_policy_benchmark(D_policy_benchmark == 1) = nan;
V_r_benchmark = reshape(arellano_solution_alpha_benchmark.V_r, params.b_grid_size, params.y_grid_size)';
V_d_benchmark = reshape(arellano_solution_alpha_benchmark.V_d, params.b_grid_size, params.y_grid_size)';
V_benchmark = reshape(arellano_solution_alpha_benchmark.V, params.b_grid_size, params.y_grid_size)';

% For Alpha Low Calibration:
Q_low = reshape(arellano_solution_alpha_low.Q, params.b_grid_size, params.y_grid_size)';
B_policy_low = reshape(arellano_solution_alpha_low.B_policy, params.b_grid_size, params.y_grid_size)'+1; %MATLAB starts indexing at 1.
D_policy_low = reshape(arellano_solution_alpha_low.D_policy, params.b_grid_size, params.y_grid_size)';
B_policy_low(D_policy_low == 1) = nan;
V_r_low = reshape(arellano_solution_alpha_low.V_r, params.b_grid_size, params.y_grid_size)';
V_d_low = reshape(arellano_solution_alpha_low.V_d, params.b_grid_size, params.y_grid_size)';
V_low = reshape(arellano_solution_alpha_low.V, params.b_grid_size, params.y_grid_size)';

% For Alpha High Calibration:
Q_high = reshape(arellano_solution_alpha_high.Q, params.b_grid_size, params.y_grid_size)';
B_policy_high = reshape(arellano_solution_alpha_high.B_policy, params.b_grid_size, params.y_grid_size)'+1; %MATLAB starts indexing at 1.
D_policy_high = reshape(arellano_solution_alpha_high.D_policy, params.b_grid_size, params.y_grid_size)';
B_policy_high(D_policy_high == 1) = nan;
V_r_high = reshape(arellano_solution_alpha_high.V_r, params.b_grid_size, params.y_grid_size)';
V_d_high = reshape(arellano_solution_alpha_high.V_d, params.b_grid_size, params.y_grid_size)';
V_high = reshape(arellano_solution_alpha_high.V, params.b_grid_size, params.y_grid_size)';


%% plot: Value functions:

y_choice = 30;
figure('Position', [100, 100, 800, 600]); % Adjust the figure size as needed
hold on;
plot(B_grid, V_benchmark(y_choice, :), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle','--'); % Gray line
plot(B_grid, V_low(y_choice, :), 'r-.', 'LineWidth', 2);
plot(B_grid, V_high(y_choice, :), 'b', 'LineWidth', 2);
hold off;
xlabel('Bond Holdings', 'FontSize', 12);
ylabel('Value Functions', 'FontSize', 12);
title('Low recovery rate vs high recovery rate', 'FontSize', 14);
% Create legend with labels and numerical values
legend({['$\alpha_{base} = ' num2str(params_benchmark.alpha) '$'],  ['$\alpha_{low} = ' num2str(params_low.alpha) '$'], ['$\alpha_{high} = ' num2str(params_high.alpha) '$'],}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Southeast');
% Adjust label positions
xlabelPos = get(gca, 'Position');
xlabelPos(1) = 0.9 * xlabelPos(1);
ylabelPos = get(gca, 'Position');
ylabelPos(2) = 0.9 * ylabelPos(2);
set(gca, 'Position', xlabelPos);
set(gca, 'Position', ylabelPos);
% Set x-axis limits to start from -0.5
xlim([-0.6, max(B_grid)]);
% Set background color to white
set(gcf, 'Color', 'w');
% Save the figure as a PNG
saveas(gcf, 'Value_functions.png');


%% Plot: Prices:

% Plot for Alpha High Calibration
figure('Position', [100, 100, 800, 600]); % Adjust the figure size as needed
hold on;
plot(B_grid, Q_benchmark(y_choice, :), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle','--'); % Gray line
plot(B_grid, Q_low(y_choice, :), 'r-.', 'LineWidth', 2);
plot(B_grid, Q_high(y_choice, :), 'b', 'LineWidth', 2);
hold off;
xlabel('Bond Holdings', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$q(b^{\prime}, y)$', 'Interpreter', 'latex', 'FontSize', 12);
title('Price Schedule for Different Recovery Rates', 'FontSize', 14);
% Create legend with labels and numerical values
legend({['$\alpha_{base} = ' num2str(params_benchmark.alpha) '$'],  ['$\alpha_{low} = ' num2str(params_low.alpha) '$'], ['$\alpha_{high} = ' num2str(params_high.alpha) '$'],}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Southeast');
% Adjust label positions
xlabelPos = get(gca, 'Position');
xlabelPos(1) = 0.9 * xlabelPos(1);
ylabelPos = get(gca, 'Position');
ylabelPos(2) = 0.9 * ylabelPos(2);
set(gca, 'Position', xlabelPos);
set(gca, 'Position', ylabelPos);
% Set x-axis limits to start from -0.5
xlim([-0.6, max(B_grid)]);
% Set background color to white
set(gcf, 'Color', 'w');
% Save the figure as a PNG
saveas(gcf, 'Prices.png');


