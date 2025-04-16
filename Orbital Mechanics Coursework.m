% Housekeeping
clear
clc

%% PART 1 (i)

% Define variables given in question
mu = 398600; % km^3/s^2
r_0 = [2031; 5599; 2775]; % km
v_0 = [-0.8419; -3.204; 7.081]; % km/s
a_1 = 15272; % km
e_1 = 0.3187;
i_1 = 88.39; % degrees
Omega_1 = 65.93; % degrees
omega_1 = 337.62; % degrees
theta_1 = 248.68; % degrees

% Calculating position in perifocal frame
p_1 = a_1 * (1 - e_1^2); % Semi-latus rectum, km
r_1p = p_1 / (1 + e_1 * cosd(theta_1)) .* [cosd(theta_1); sind(theta_1)]; % Position in perifocal frame, km

% Calculating velocity in perifocal frame
h_1 = sqrt(mu * p_1); % Specific angular momentum, km^2/s
v_1p = (mu / h_1) * [-cosd(theta_1); sind(theta_1) + e_1]; % Velocity in perifocal frame, km/s

% Print answer
disp('Part 1(i)')
disp(['Position vector in perifocal frame, r_1,p = ' num2str(r_1p(1)) 'i_e + ' num2str(r_1p(2)) 'i_p km'])
disp(['Velocity vector in perifocal frame, v_1,p = ' num2str(v_1p(1)) 'i_e + ' num2str(v_1p(2)) 'i_p km/s'])
disp(newline)

%% PART 1 (ii)

% Prepare rotation matrix for conversion of Cartesian unit vectors to
% perifocal unit vectors
R_omega = [cosd(omega_1) sind(omega_1) 0; -sind(omega_1) cosd(omega_1) 0; 0 0 1]; % Rotation matrix about omega
R_i = [1 0 0; 0 cosd(i_1) sind(i_1); 0 -sind(i_1) cosd(i_1)]; % Rotation matrix about i
R_Omega = [cosd(Omega_1) sind(Omega_1) 0; -sind(Omega_1) cosd(Omega_1) 0; 0 0 1]; % Rotation matrix about Omega
R = R_omega * R_i * R_Omega;

% Apply rotation matrices and find perifocal unit vectors

i_e_hat = R' * [0; 1; 0]; % i_e axis in Cartesian coordinates
i_w_hat = R' * [0; 0; 1]; % i_p axis in Cartesian coordinates
i_p_hat = R' * [1; 0; 0]; % i_w axis in Cartesian coordinates

% Print answer
disp('Part 1(ii)')
disp('Perifocal reference frame unit vectors:')
disp(['i_e = [' num2str(i_e_hat(1)) ', ' num2str(i_e_hat(2)) ', ' num2str(i_e_hat(3)) ']'])
disp(['i_w = [' num2str(i_w_hat(1)) ', ' num2str(i_w_hat(2)) ', ' num2str(i_w_hat(3)) ']'])
disp(['i_p = [' num2str(i_p_hat(1)) ', ' num2str(i_p_hat(2)) ', ' num2str(i_p_hat(3)) ']'])
disp(newline)

%% PART 2 (i)

% Define values given in question
r_1 = [-4085.1; -9918.1; -11215.3]; % km
v_1 = [1.911; 4.131; -2.135]; % km/s
Delta_t = 6.753 * 60^2; % seconds

% Calculate change in true anomaly
Delta_theta = acosd(dot(r_0,r_1) / (norm(r_0) * norm(r_1))); % degrees

% Calculate chord length
c = sqrt(norm(r_1)^2 + norm(r_0)^2 -2 * norm(r_1) * norm(r_0) * cosd(Delta_theta)); % km

% Calculate minimum-energy semi-major axis
a_m = (norm(r_0) + norm(r_1) + c) / 4; % km

% Print answer
disp('Part 2(i)')
disp(['Minimum-energy semi-major axis, a_m = ' num2str(a_m) ' km'])
disp(newline)

%% PART 2 (ii)

s = 2 * a_m; % Minimum-energy major axis, km

% Calculate transfer time for parabolic orbit
Delta_t_p = (sqrt(2 / mu) / 3) * (s^(3/2) - sind(Delta_theta) * (s - c)^(3/2)); % seconds

% Delta_t_p < Delta_t provided, so this corresponds to an elliptical
% transfer orbit

% Print answer
disp('Part 2(ii)')
disp(['Transfer time for parabolic orbit, Delta_t_p = ' num2str(Delta_t_p) ' secs'])
disp(['Transfer time provided, Delta_t = ' num2str(Delta_t) ' secs'])
disp('Delta_t > Delta_t_p so this represents an elliptical transfer orbit')
disp(newline)

%% PART 2 (iii)

beta_m = 2 * asin(sqrt((s - c) / (2 * a_m))); % rads

% Calculate Delta_t for minimum energy semi-major axis
% Use positive value for short way
Delta_t_m = sqrt(s^3 / (8 * mu)) * (pi - beta_m + sin(beta_m)); % seconds

% Delta_t > Delta_t_m and Delta_theta < pi
% So sin(alpha) < 0, sin(beta) >= 0
% and alpha = 2*pi - alpha_0
% beta = beta_0

% Print answer
disp('Part 2(iii)')
disp(['Transfer time for minimum energy orbit, Delta_t_m = ' num2str(Delta_t_m) ' secs'])
disp('Delta_t > Delta_t_m and Delta_theta < pi')
disp('So sin(alpha) < 0, sin(beta) >= 0 and alpha = 2*pi - alpha_0, beta = beta_0 ')
disp(newline)

%% PART 2 (iv)

% Initialise functions to solve Lambert's equation

% Auxiliary angles
alpha = @(a) (2 * pi - 2 * asin(sqrt(s / (2 * a))));
beta = @(a) (2 * asin(sqrt((s - c) / (2 * a))));

% Function for fsolve (rearranged Lambert's equation)
func = @(a) ((a^3 / mu)^0.5 * (alpha(a) - beta(a) - (sin(alpha(a)) - sin(beta(a)))) - Delta_t);

% Calculate semi-major axis for transfer orbit, a_01
a_init = 2*a_m; % Initial guess for a_01, km
a_01 = fzero(func,a_init); % Find a_01 using fzero function, km

% Print answer
disp('Part 2(iv)')
disp(['Semi-major axis of transfer orbit, a_01 = ' num2str(a_01) ' km'])
disp(newline)

%% PART 2 (v)

% Calculate alpha and beta for transfer orbit
alpha_01 = 2 * pi - 2 * asin(sqrt(s / (2 * a_01))); % rads
beta_01 =  2 * asin(sqrt((s - c) / (2 * a_01))); % rads

% Calculate A and B
A = sqrt(mu / (4 * a_01)) * cot(alpha_01 / 2);
B = sqrt(mu / (4 * a_01)) * cot(beta_01 / 2);

% Find normalised position vectors at start and end of transfer orbit
u_0 = r_0 / norm(r_0);
u_1 = r_1 / norm(r_1);

% Find normalised chord vector
u_c = (r_1 - r_0) / c;

% Calculate terminal velocity vectors
v_0_transfer = (B + A) * u_c + (B - A) * u_0; % km/s
v_1_transfer = (B + A) * u_c - (B - A) * u_1; % km/s

% Calculate velocity impulse vectors
Delta_v_0 = v_0_transfer - v_0; % km/s
Delta_v_1 = v_1_transfer - v_1; % km/s

% Print answer
disp('Part 2(v)')
disp(['Departure velocity impulse vector, Delta_v_0 = [' num2str(Delta_v_0(1)) ', ' num2str(Delta_v_0(2)) ', ' num2str(Delta_v_0(3)) '] km/s'])
disp(['Arrival velocity impulse vector, Delta_v_1 = [' num2str(Delta_v_1(1)) ', ' num2str(Delta_v_1(2)) ', ' num2str(Delta_v_1(3)) '] km/s'])