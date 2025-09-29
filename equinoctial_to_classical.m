function elements1 = equinoctial_to_classical(sigma_point)
% EQUINOCTIAL_TO_CLASSICAL Converts equinoctial coordinates to classical orbital elements
%
% Inputs (sigma_point vector):
% a - Semi-major axis [m]
% h - Eccentricity component h = e * sin(omega + Omega) [-]
% k - Eccentricity component k = e * cos(omega + Omega) [-]
% p - Inclination component p = tan(i/2) * sin(Omega) [-]
% q - Inclination component q = tan(i/2) * cos(Omega) [-]
% l - Mean longitude [deg]

% Extract equinoctial elements
a = sigma_point(1);
h = sigma_point(2);
k = sigma_point(3);
p = sigma_point(4);
q = sigma_point(5);
l = sigma_point(6);

% Standard gravitational parameter for Earth (m^3/s^2)
mu = 3.986004418e14;

% Convert mean longitude to radians if it's in degrees
if abs(l) > 2*pi
    l_rad = deg2rad(l);
else
    l_rad = l;  % Already in radians
end

% Calculate eccentricity
elements1.eccentricity = sqrt(h^2 + k^2);
e = elements1.eccentricity;

% Calculate inclination
tan_half_i = sqrt(p^2 + q^2);
i_rad = 2 * atan(tan_half_i);
elements1.inclination = rad2deg(i_rad);

% Calculate RAAN (Right Ascension of Ascending Node)
if tan_half_i == 0
    % Equatorial orbit - RAAN is undefined, set to 0
    Omega_rad = 0;
else
    Omega_rad = atan2(p, q);
end
elements1.raan = rad2deg(Omega_rad);

% Calculate argument of periapsis + RAAN (longitude of periapsis)
if e == 0
    % Circular orbit - argument of periapsis is undefined
    varpi_rad = 0;
    omega_rad = 0;
else
    varpi_rad = atan2(h, k);  % This is omega + Omega
    omega_rad = varpi_rad - Omega_rad;  % Extract omega
end
elements1.argp = rad2deg(omega_rad);

% Calculate mean anomaly
% For equinoctial elements: M = l - varpi = l - (omega + Omega)
M_rad = l_rad - varpi_rad;
elements1.meanAnomaly = rad2deg(M_rad);

% Calculate mean motion from semi-major axis using Kepler's third law
n_rad = sqrt(mu / a^3); % rad/s
elements1.meanmotion = rad2deg(n_rad); % deg/s

% Normalize angles to [0, 360) degrees
elements1.raan = mod(elements1.raan, 360);
elements1.argp = mod(elements1.argp, 360);
elements1.meanAnomaly = mod(elements1.meanAnomaly, 360);

% Ensure positive angles
if elements1.raan < 0
    elements1.raan = elements1.raan + 360;
end
if elements1.argp < 0
    elements1.argp = elements1.argp + 360;
end
if elements1.meanAnomaly < 0
    elements1.meanAnomaly = elements1.meanAnomaly + 360;
end

end