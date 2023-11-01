%%% This file creates an Mex file using the source code from C++ %%%

%% Clear all
clc; 
clear;

%% mex 

mex main.cpp economy.cpp initialization.cpp auxiliary.cpp 

%% Common parameters:
params.b_grid_size = 51;               % Number of points in the grid for the bond price.
params.b_grid_min = -0.6;               % Minimum value of the bond price.
params.b_grid_max = 0.00;               % Maximum value of the bond price.
params.y_grid_size = 6;                 % Number of points in the grid for the income.
params.y_default = 0.969;               % Maximum income under default.
params.beta = 0.953;                    % Discount factor.
params.gamma = 2;                       % Risk aversion.
params.r = 0.017;                       % Interest rate.
params.rho = 0.945;                     % Persistence of the income.
params.sigma = 0.025;                   % Standard deviation of the income.
params.theta = 0.282;                   % Probability of a re-entry.
params.max_iter = 2000;                 % Maximum number of iterations.
params.tol = 1e-7;                      % Tolerance for the convergence.
params.m = 3;                           % Number of standard deviations for the income grid.
params.alpha_low = 0.15;                % Low recovery on defaulted debt.
params.alpha_high = 0.3;                % High recovery on defaulted debt.

%% Run code with both alphas;

tic;
calibrated_model_solution = main(params);
save('Solution', 'calibrated_model_solution')
save('Parameters', 'params')
toc;

%% Format variables:

% Exogenous:
Y_grid = calibrated_model_solution.Y_grid;
Y_grid_default = calibrated_model_solution.Y_grid_default;
B_grid = calibrated_model_solution.B_grid;
P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';

%Endogenous: this stores everthing in 3-d matrices. Each matrix is:
%   |1, 2, 3, ..., nb  
%   |nb,nb+1,....,2nb-1,
%   |2nb,....

Q_low = permute(reshape(calibrated_model_solution.Q_low, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
Q_high = permute(reshape(calibrated_model_solution.Q_high, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V = permute(reshape(calibrated_model_solution.V, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V_r = permute(reshape(calibrated_model_solution.V_r, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V_d = permute(reshape(calibrated_model_solution.V_d, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
B_policy_low = permute(reshape(calibrated_model_solution.B_policy_low, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
B_policy_high = permute(reshape(calibrated_model_solution.B_policy_high, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);

%% Problems with the code: you tend to borrow as much as you can, close to default.

plot(B_grid, B_policy_high(50,:,4))




