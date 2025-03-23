%% Non-Ideal Ramjet Analysis
% This script performs analysis of a ramjet engine considering non-ideal
% effects such as:
%   - Inlet efficiency
%   - Combustion efficiency
%   - Nozzle efficiency
%   - Pressure losses
%
% The analysis calculates performance parameters including:
%   - Thrust
%   - Specific impulse
%   - Thermal efficiency
%   - Propulsive efficiency
%   - Overall efficiency
%
% Author: Bryan Kikuta
% Date: 22/03/25

% Values
Q_R = 43500000;  % [J/kg] for JP-7 jet fuel
c_p = 1005;
R = 287;

% Pressure ratios
pi_d = 0.7926;  % Inlet
pi_b = 0.9;  % Burner
pi_n = 0.9;  % Nozzle

% Efficiencies
eta_b = 0.9;  % Burner efficiency
eta_n = 0.9;  % Nozzle efficiency

% Atmospheric conditions
rho_a = 1.225;  % [kg/m^3]
P_a = 101300;  % [kPa]
T_a = 273;  % [K]
gamma = 1.4;
M_1 = 3.2;
P_0a = P_a * (1 + (gamma - 1) / 2 * M_1^2)^(gamma / (gamma - 1));
T_0a = T_a * (1 + (gamma - 1) / 2 * M_1^2);

% Compression process
P_02 = pi_d * P_0a;
T_02 = T_0a;

% After inlet to before combustion chamber
P_03 = P_02;
T_03 = T_02;

% Combustion process
P_04 = pi_b * P_03;
T_04 = 1800;
f = ((T_04 / T_03) - 1) / ((eta_b * Q_R / (c_p * T_03)) - (T_04 / T_03));

% Expansion process
P_0e = pi_n * P_04;
T_0e = T_04;
P_e = eta_n * P_a;

M_e = sqrt((2 / (gamma - 1)) * ((1 + (gamma - 1) / 2 * M^2) * (pi_d * pi_b * pi_n * P_a / P_e)^((gamma - 1) / gamma) - 1));

T_e = T_0e / (1 + (gamma - 1) / 2 * M_e^2);
a_e = sqrt(gamma * R * T_e);
U_e = M_e * a_e;

% Determine thrust required
L_D = 7;  % Based on L/D of the XB-70 at Mach 3
n_engines = 6;
W = 2418421;  % [N] Based on XB-70 loaded weight
Th = W / L_D / n_engines;  % [N] Required thrust given weight, L/D and number of engines
