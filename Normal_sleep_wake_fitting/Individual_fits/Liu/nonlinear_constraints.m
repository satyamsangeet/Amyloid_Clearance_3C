function [c, ceq] = nonlinear_constraints(params)
    r_bc     = params(1);
    r_bp     = params(2);
    sigma_bc = params(4);
    sigma_bp = params(5);

    % c <= 0 required by fmincon
    c(1) = r_bc + r_bp - 0.25;                            % wake clearance cap
    c(2) = sigma_bc * r_bc + sigma_bp * r_bp - 0.5;       % sleep clearance cap

    ceq = [];
end