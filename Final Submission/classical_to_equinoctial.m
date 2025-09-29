function [equinoctial, P] = classical_to_equinoctial(elements1)
    % CLASSICAL_TO_EQUINOCTIAL Converts classical orbital elements to equinoctial coordinates
    %
    % Inputs:
    %   elements1 - Structure with fields:
    %     .meanmotion or .MeanMotion - Mean motion [deg/s]
    %     .eccentricity - Eccentricity [-]
    %     .inclination - Inclination [deg]
    %     .raan - Right Ascension of Ascending Node [deg]
    %     .argp - Argument of Periapsis [deg]
    %     .meanAnomaly - Mean Anomaly [deg]
    %
    % Outputs:
    %   a - Semi-major axis [m]
    %   h - Eccentricity component h = e * sin(argp + raan) [-]
    %   k - Eccentricity component k = e * cos(argp + raan) [-]
    %   p - Inclination component p = tan(i/2) * sin(raan) [-]
    %   q - Inclination component q = tan(i/2) * cos(raan) [-]
    %   l - Mean longitude l = raan + argp + meanAnomaly [deg]
    
    % Standard gravitational parameter for Earth (m^3/s^2)

    equinoctial = struct();
    uncertainty = struct();

    mu = 3.986004418e14;
    
    % Extract classical elements and handle both naming conventions
    if isfield(elements1, 'meanmotion')
        n = elements1.meanmotion; % deg/s
    elseif isfield(elements1, 'MeanMotion')
        n = elements1.MeanMotion; % deg/s
    else
        error('Mean motion field not found. Expected "meanmotion" or "MeanMotion"');
    end
    
    e = elements1.Eccentricity;       % eccentricity
    i = elements1.Inclination;        % deg
    Omega = elements1.RightAscensionOfAscendingNode;           % deg
    omega = elements1.ArgumentOfPeriapsis;           % deg
    M = elements1.MeanAnomaly;        % deg
    
    % Convert angular quantities to radians for calculations
    n_rad = deg2rad(n);               % rad/s
    i_rad = deg2rad(i);               % rad
    Omega_rad = deg2rad(Omega);       % rad
    omega_rad = deg2rad(omega);       % rad
    M_rad = deg2rad(M);               % rad
    
    % Calculate semi-major axis from mean motion using Kepler's third law
    % n = sqrt(mu/a^3) -> a = (mu/n^2)^(1/3)
    a = (mu / n_rad^2)^(1/3);         % m
   

    % Calculate equinoctial elements
    h = e * sin(omega_rad + Omega_rad);
    k = e * cos(omega_rad + Omega_rad);
    p = tan(i_rad/2) * sin(Omega_rad);
    q = tan(i_rad/2) * cos(Omega_rad);
    
    % Mean longitude (convert back to degrees)
    l_rad = Omega_rad + omega_rad + M_rad;
    l = rad2deg(l_rad);
    
    % Ensure mean longitude is in [0, 360) range
    l = mod(l, 360);

    equinoctial.a = a;
    equinoctial.h = h;
    equinoctial.k = k;
    equinoctial.p = p;
    equinoctial.q = q;
    equinoctial.l = l;

    %{
    perc.a = 0.0005;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.005;  % 0.5% of l

    % Calculate uncertainties as percentages of nominal values
    uncertainty.a = perc.a * equinoctial.a;
    uncertainty.h = perc.h * equinoctial.h;
    uncertainty.k = perc.k * equinoctial.k;
    uncertainty.p = perc.p * equinoctial.p;
    uncertainty.q = perc.q * equinoctial.q;
    uncertainty.l = perc.l * equinoctial.l;
    %}

    uncertainty.a = 100;
    uncertainty.h = 0.0001;
    uncertainty.k = 0.0001;
    uncertainty.p = 0.0001;
    uncertainty.q = 0.0001;
    uncertainty.l = 0.0001;
    
    P = diag([uncertainty.a^2, ...
              uncertainty.h^2, ...
              uncertainty.k^2, ...
              uncertainty.p^2, ...
              uncertainty.q^2, ...
              uncertainty.l^2]);
end
