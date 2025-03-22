function [betas, deflection_angles, total_deflection_angles, mach_numbers, pressure_ratios, total_pressure_ratio] = ramp_angle_calc()
    % Inputs
    M_1 = 3.2;  % Incoming Mach number
    M_n = 1.3;  % Mach number before normal shock
    n = 3;  % Number of oblique shocks
    gamma = 1.4;  % Gas constant

    % Initial guess for betas in degrees
    beta_initial = 10 * ones(1, n);

    % Lower and upper bounds for betas in degrees
    lb = 0 * ones(1, n);
    ub = 90 * ones(1, n);

    % Solve using lsqnonlin with bounds
    options = optimset('Display', 'iter');  % Display iterations
    betas = lsqnonlin(@(beta) shock_equations(beta, M_1, M_n, n, gamma), beta_initial, lb, ub, options);

    % Compute final Mach numbers and deflection angles
    [mach_numbers, deflection_angles, total_deflection_angles, pressure_ratios] = compute_final_values(M_1, betas, n, gamma);

    % Display results
    fprintf('Oblique Shock Angles (in degrees):\n');
    for i = 1:n
        fprintf('Beta %d: %.4f°\n', i, betas(i));
    end

    fprintf('\nDeflection Angles (in degrees):\n');
    for i = 1:n
        fprintf('Delta %d: %.4f°\n', i, deflection_angles(i));
    end

    fprintf('\nTotal Deflection Angles (in degrees):\n');
    for i = 1:n
        fprintf('Delta %d: %.4f°\n', i, total_deflection_angles(i));
    end

    fprintf('\nFinal Mach Numbers:\n');
    for i = 1:n+2
        fprintf('M%d: %.4f\n', i, mach_numbers(i));
    end

    fprintf('\nPressure Ratios:\n');
    for i = 1:n+1
        fprintf('π%d: %.4f\n', i, pressure_ratios(i));
    end

    % Compute total pressure recovery ratio
    total_pressure_ratio = prod(pressure_ratios);
    fprintf('\nTotal Pressure Recovery Ratio (πd): %.4f\n', total_pressure_ratio);
end

% Function for computing deflection angle for given oblique shock angle
function delta = compute_deflection(M_1, gamma, beta)
    tan_delta = (2 * cotd(beta) * (M_1^2 * sind(beta)^2 - 1)) / ((gamma + 1) * M_1^2 - 2 * (M_1^2 * sind(beta)^2 - 1));
    delta = atand(tan_delta);
end

% Function for computing new mach number after normal shock
function M_2 = compute_mach_after_normal(M_1, gamma)
    M_2_2 = (M_1^2 + 2 / (gamma - 1)) / (2 * gamma / (gamma - 1) * M_1^2 - 1);
    M_2 = sqrt(M_2_2);
end

% Function for computing new mach number after oblique shock
function M_2 = compute_mach_after_oblique(M_1, gamma, beta)
    % Compute deflection angle after oblique shock
    delta = compute_deflection(M_1, gamma, beta);

    % Compute normal component of M_1
    M_1n = compute_normal_component(M_1, beta);

    % Compute normal of M_2 by treating as normal shock with M_1n
    M_2n = compute_mach_after_normal(M_1n, gamma);

    % Compute M_2 from the normal of M_2
    M_2 = M_2n / sind(beta - delta);
end

% Function for computing the normal mach number to a shock
function M_n = compute_normal_component(M, beta)
    % Compute normal component
    M_n = M * sind(beta);
end

% Function for computing pressure ratio across a shock
function P_r = compute_pressure_ratio(M_1, beta, gamma)
    % Compute normal component of mach
    M_1n = compute_normal_component(M_1, beta);

    % Compute the pressure ratio across normal component
    num = (((gamma + 1) / 2 * M_1n^2) / (1 + (gamma - 1) / 2 * M_1n^2))^(gamma / (gamma - 1));
    den = (2 * gamma / (gamma + 1) * M_1n^2 - (gamma - 1) / (gamma + 1))^(1 / (gamma - 1));
    P_r = num / den;
end

% Function to solve for oblique shock angles
function F = shock_equations(beta, M_1, M_n, n, gamma)
    % Preallocate the array for downstream Mach numbers
    M = zeros(1, n + 1);
    M(1) = M_1;  % Initialize the first Mach number

    % Compute downstream Mach numbers using oblique shock function
    for i = 1:n
        M(i + 1) = compute_mach_after_oblique(M(i), gamma, beta(i));  % Compute the next Mach number
    end

    % Preallocate the array for the system of equations
    F = zeros(1, n);
    F(1) = M(end) - M_n;  % Flow from 3rd oblique shock should equal M_n

    % Define the system of equations
    for i = 2:n
        F(i) = M(i - 1) * sind(beta(i - 1)) - M(i) * sind(beta(i));  % Equation for oblique shock angles
    end
end

% Function to compute final Mach numbers and deflection angles
function [mach_numbers, deflection_angles, total_deflection_angles, pressure_ratios] = compute_final_values(M_1, betas, n, gamma)
    mach_numbers = zeros(1, n + 2);
    deflection_angles = zeros(1, n);
    total_deflection_angles = zeros(1, n);
    pressure_ratios = zeros(1, n + 1);
    mach_numbers(1) = M_1;

    % Compute values across oblique
    for i = 1:n
        mach_numbers(i + 1) = compute_mach_after_oblique(mach_numbers(i), gamma, betas(i));
        deflection_angles(i) = compute_deflection(mach_numbers(i), gamma, betas(i));
        total_deflection_angles(i) = sum(deflection_angles(1:i));
        pressure_ratios(i) = compute_pressure_ratio(mach_numbers(i), betas(i), gamma);
    end

    % Compute values across normal
    mach_numbers(n + 2) = compute_mach_after_normal(mach_numbers(n + 1), gamma);
    pressure_ratios(n + 1) = compute_pressure_ratio(mach_numbers(n + 1), 90, gamma);
end