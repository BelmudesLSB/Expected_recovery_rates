%%% This file creates an Mex file using the source code from C++ %%%

%% Clear all
clc; 
clear;

%% mex
mex CXXFLAGS="$CXXFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" main.cpp economy.cpp initialization.cpp auxiliary.cpp 

%% Common parameters:
params.b_grid_size_lowr = 75;           % Number of points in the grid for the bond price.
params.b_grid_size_highr = 121;
params.b_grid_min_lowr = -0.8;          % Minimum value of the bond price.
params.b_grid_min_highr = -1.5;
params.b_grid_max_lowr = 0.0;           % Maximum value of the bond price.
params.b_grid_max_highr = 0.0;
params.y_grid_size = 21;               % Number of points in the grid for the income.
params.y_default = 0.969;              % Maximum income under default.
params.beta = 0.953;                   % Discount factor.
params.gamma = 2;                      % Risk aversion.
params.r = 0.017;                      % Interest rate.
params.rho = 0.945;                    % Persistence of the income.
params.sigma = 0.025;                  % Standard deviation of the income.
params.theta = 0.282;                  % Probability of a re-entry.
params.max_iter = 1000;                 % Maximum number of iterations.
params.tol = 1e-7;                     % Tolerance for the convergence.
params.m = 3;                          % Number of standard deviations for the income grid.
params.alpha_lowr = 0.15;              % Low recovery on defaulted debt.
params.alpha_highr = 0.3;               % High recovery on defaulted debt.

%% Run code with both alphas;
tic;
calibrated_model_solution = main(params);
save('Solution', 'calibrated_model_solution')
save('Parameters', 'params')
toc;

%% Perform checks to see if bounds bind or not:
calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy==1) = nan;
min(calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy==0))
find(min(calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy==0)))
%% Format variables:

% Exogenous:
Y_grid = calibrated_model_solution.Y_grid;
Y_grid_default = calibrated_model_solution.Y_grid_default;
B_grid_lowr = calibrated_model_solution.B_grid_lowr;
B_grid_highr = calibrated_model_solution.B_grid_highr;
P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';
% Endogenous objects:
Q_lowr = permute(reshape(calibrated_model_solution.Q_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Q_highr = permute(reshape(calibrated_model_solution.Q_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
V = permute(reshape(calibrated_model_solution.V, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
V_r = permute(reshape(calibrated_model_solution.V_r, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
V_d = permute(reshape(calibrated_model_solution.V_d, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
B_policy_lowr = permute(reshape(calibrated_model_solution.B_policy_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
B_policy_highr = permute(reshape(calibrated_model_solution.B_policy_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
B_policy_highr(D_policy == 1) = nan;
B_policy_lowr(D_policy==1) = nan;
B_policy_lowr = B_policy_lowr + 1;
B_policy_highr = B_policy_highr + 1;

%% PLots for presentation:

% Assuming B_grid_highr contains your x-values
% V(:,params.b_grid_size_lowr, y_choice) represents the y-values for the chosen y_choice
y_choice = 12;
x_values = B_grid_highr;
y_values = V(:, params.b_grid_size_lowr, y_choice);
% Define the x-range from -0.8 to 0
x_range = [-0.8, 0];
% Find the indices that fall within the specified x-range
indices = x_values >= x_range(1) & x_values <= x_range(2);
% Plot the selected data
plot(x_values(indices), y_values(indices));
% Add labels and title if desired
xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Plot Title');

% Assuming B_grid_highr contains your x-values
% V(:,params.b_grid_size_lowr, y_choice) represents the y-values for the chosen y_choice
y_choice = 12;
x_values = B_grid_lowr;
y_values = V(params.b_grid_size_highr, :, y_choice);
% Define the x-range from -0.8 to 0
x_range = [-0.8, 0];
% Find the indices that fall within the specified x-range
indices = x_values >= x_range(1) & x_values <= x_range(2);
% Plot the selected data
plot(x_values(indices), y_values(indices));
% Add labels and title if desired
xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Plot Title');



