function [next_Particle_ECI, Weights] = Generate_Monte_Carlo(mu, P, startTime, stopTime, sampleTime, numSamples, File)
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

fprintf('\n=== Processing Multiple N Values ===\n');

% Time the particle generation
fprintf('Generating particles...\n');
[Weights, Equinoctial_Particles] = GenerateGSF(mu, P, numSamples);

numParticles = numSamples;

% Pre-allocate the output array
next_Particle_ECI = zeros(6, numTimeSteps, numParticles);

% Use parfor for the main particle propagation loop
parfor i = 1:numParticles
    classical_elements = equinoctial_to_classical(Equinoctial_Particles(:, i));
    
    % Create a local satellite scenario for this worker
    sc_local = satelliteScenario(startTime, stopTime, sampleTime);
    
    Keplerian_Particles_i = [
        classical_elements.meanmotion; % mean motion
        classical_elements.eccentricity; % eccentricity
        classical_elements.inclination; % inclination
        classical_elements.raan; % right ascension of ascending node
        classical_elements.argp; % argument of periapsis
        classical_elements.meanAnomaly % mean anomaly
    ];
    
    elems = struct( ...
        'meanMotion', Keplerian_Particles_i(1) * 86400 / 360, ... % rev/day (primary input)
        'eccentricity', Keplerian_Particles_i(2), ...
        'inclination', Keplerian_Particles_i(3), ...
        'raan', Keplerian_Particles_i(4), ...
        'argPeriapsis', Keplerian_Particles_i(5), ...
        'meanAnomaly', Keplerian_Particles_i(6) ...
    );
    
    % Create unique filename for each sigma point and component
    outputFile = 'updated.xml';
    updateOMMFile(File, outputFile, elems);
    sat = satellite(sc_local, outputFile);
    
    el = orbitalElements(sat);
    
    % Get states for the new scenario
    [positions, velocities, times] = states(sat);
    
    % Store the propagated states
    next_Particle_ECI(:, :, i) = [positions; velocities];
    delete(outputFile)
end

end
function [Weights, Equinoctial_Particles] = GenerateGSF(mu, P, numSamples)
    data = mvnrnd(mu', P, numSamples); % Note: mu' to make it row vector
    Equinoctial_Particles = data';
    Weights = mvnpdf(data, mu', P);
    Weights = Weights/sum(Weights);
end