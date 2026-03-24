function [lb, ub, initial_guess] = config()
% lower and upper bounds edits
% initial_guess is always computed as the midpoint.

%      r_bc  r_bp  r_cp  sigma_bc  sigma_bp  sigma_cp  sigma_p  A  sigma_A  r_p
lb = [0.01, 0.01, 0,   1, 1, 1, 1, 0,   0,    0  ];
ub = [0.25, 0.25, 0.1, 7, 7, 7, 7, 111, 0.99, 0.6];
initial_guess = 0.5 * (lb + ub);

end