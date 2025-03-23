function [betas, deflection_angles, total_deflection_angles, mach_numbers, pressure_ratios, total_pressure_ratio] = ramp_angle_calc()
    %RAMP_ANGLE_CALC Calculate oblique shock angles and properties for a multi-ramp inlet
    %   [BETAS, DEFLECTION_ANGLES, TOTAL_DEFLECTION_ANGLES, MACH_NUMBERS, PRESSURE_RATIOS, TOTAL_PRESSURE_RATIO] = RAMP_ANGLE_CALC()
    %
    %   Calculates the optimal oblique shock angles for a multi-ramp inlet 
    %   design that maximizes pressure recovery.
    %
    %   Outputs:
    %       BETAS - Oblique shock angles in degrees
    %       DEFLECTION_ANGLES - Flow deflection angles at each ramp in degrees
    %       TOTAL_DEFLECTION_ANGLES - Cumulative flow deflection angles in degrees
    %       MACH_NUMBERS - Mach numbers before and after each shock
    %       PRESSURE_RATIOS - Pressure ratios across each shock
    %       TOTAL_PRESSURE_RATIO - Overall pressure recovery ratio across all shocks
    
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

function delta = compute_deflection(M_1, gamma, beta)
    %COMPUTE_DEFLECTION Calculate flow deflection angle for an oblique shock
    %   DELTA = COMPUTE_DEFLECTION(M_1, GAMMA, BETA)
    %
    %   Calculates the flow deflection angle using the theta-beta-M relation
    %   for oblique shocks.
    %
    %   Inputs:
    %       M_1 - Upstream Mach number
    %       GAMMA - Ratio of specific heats (cp/cv)
    %       BETA - Oblique shock angle in degrees
    %
    %   Output:
    %       DELTA - Flow deflection angle in degrees
    
    tan_delta = (2 * cotd(beta) * (M_1^2 * sind(beta)^2 - 1)) / ((gamma + 1) * M_1^2 - 2 * (M_1^2 * sind(beta)^2 - 1));
    delta = atand(tan_delta);
end

function M_2 = compute_mach_after_normal(M_1, gamma)
    %COMPUTE_MACH_AFTER_NORMAL Calculate downstream Mach number after a normal shock
    %   M_2 = COMPUTE_MACH_AFTER_NORMAL(M_1, GAMMA)
    %
    %   Calculates the Mach number downstream of a normal shock using
    %   the Rankine-Hugoniot relation.
    %
    %   Inputs:
    %       M_1 - Upstream Mach number
    %       GAMMA - Ratio of specific heats (cp/cv)
    %
    %   Output:
    %       M_2 - Downstream Mach number
    
    M_2_2 = (M_1^2 + 2 / (gamma - 1)) / (2 * gamma / (gamma - 1) * M_1^2 - 1);
    M_2 = sqrt(M_2_2);
end

function M_2 = compute_mach_after_oblique(M_1, gamma, beta)
    %COMPUTE_MACH_AFTER_OBLIQUE Calculate downstream Mach number after an oblique shock
    %   M_2 = COMPUTE_MACH_AFTER_OBLIQUE(M_1, GAMMA, BETA)
    %
    %   Calculates the Mach number downstream of an oblique shock by decomposing
    %   the flow into normal and tangential components, applying normal shock
    %   relations, and then recomposing.
    %
    %   Inputs:
    %       M_1 - Upstream Mach number
    %       GAMMA - Ratio of specific heats (cp/cv)
    %       BETA - Oblique shock angle in degrees
    %
    %   Output:
    %       M_2 - Downstream Mach number
    
    % Compute deflection angle after oblique shock
    delta = compute_deflection(M_1, gamma, beta);

    % Compute normal component of M_1
    M_1n = compute_normal_component(M_1, beta);

    % Compute normal of M_2 by treating as normal shock with M_1n
    M_2n = compute_mach_after_normal(M_1n, gamma);

    % Compute M_2 from the normal of M_2
    M_2 = M_2n / sind(beta - delta);
end

function M_n = compute_normal_component(M, beta)
    %COMPUTE_NORMAL_COMPONENT Calculate normal component of Mach number
    %   M_N = COMPUTE_NORMAL_COMPONENT(M, BETA)
    %
    %   Calculates the component of Mach number normal to a shock wave.
    %
    %   Inputs:
    %       M - Mach number
    %       BETA - Shock angle in degrees
    %
    %   Output:
    %       M_N - Normal component of Mach number
    
    % Compute normal component
    M_n = M * sind(beta);
end

function P_r = compute_pressure_ratio(M_1, beta, gamma)
    %COMPUTE_PRESSURE_RATIO Calculate pressure ratio across a shock
    %   P_R = COMPUTE_PRESSURE_RATIO(M_1, BETA, GAMMA)
    %
    %   Calculates the static pressure ratio across a shock wave.
    %
    %   Inputs:
    %       M_1 - Upstream Mach number
    %       BETA - Shock angle in degrees (90 for normal shock)
    %       GAMMA - Ratio of specific heats (cp/cv)
    %
    %   Output:
    %       P_R - Static pressure ratio (p2/p1)
    
    % Compute normal component of mach
    M_1n = compute_normal_component(M_1, beta);

    % Compute the pressure ratio across normal component
    num = (((gamma + 1) / 2 * M_1n^2) / (1 + (gamma - 1) / 2 * M_1n^2))^(gamma / (gamma - 1));
    den = (2 * gamma / (gamma + 1) * M_1n^2 - (gamma - 1) / (gamma + 1))^(1 / (gamma - 1));
    P_r = num / den;
end

function F = shock_equations(beta, M_1, M_n, n, gamma)
    %SHOCK_EQUATIONS System of equations for solving oblique shock angles
    %   F = SHOCK_EQUATIONS(BETA, M_1, M_N, N, GAMMA)
    %
    %   Defines the system of equations to be minimized when solving for
    %   the optimal oblique shock angles that achieve a specified M_n.
    %   Used as an objective function for nonlinear optimization.
    %
    %   Inputs:
    %       BETA - Vector of oblique shock angles in degrees
    %       M_1 - Initial upstream Mach number
    %       M_N - Target Mach number before the terminal normal shock
    %       N - Number of oblique shocks
    %       GAMMA - Ratio of specific heats (cp/cv)
    %
    %   Output:
    %       F - Vector of equation residuals to be minimized
    
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

function [mach_numbers, deflection_angles, total_deflection_angles, pressure_ratios] = compute_final_values(M_1, betas, n, gamma)
    %COMPUTE_FINAL_VALUES Calculate all flow properties across multiple shocks
    %   [MACH_NUMBERS, DEFLECTION_ANGLES, TOTAL_DEFLECTION_ANGLES, PRESSURE_RATIOS] = COMPUTE_FINAL_VALUES(M_1, BETAS, N, GAMMA)
    %
    %   Calculates the flow properties (Mach numbers, deflection angles, etc.)
    %   across a series of oblique shocks followed by a normal shock.
    %
    %   Inputs:
    %       M_1 - Initial upstream Mach number
    %       BETAS - Vector of oblique shock angles in degrees
    %       N - Number of oblique shocks
    %       GAMMA - Ratio of specific heats (cp/cv)
    %
    %   Outputs:
    %       MACH_NUMBERS - Vector of Mach numbers before and after each shock
    %       DEFLECTION_ANGLES - Vector of flow deflection angles at each ramp
    %       TOTAL_DEFLECTION_ANGLES - Vector of cumulative flow deflection angles
    %       PRESSURE_RATIOS - Vector of pressure ratios across each shock
    
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