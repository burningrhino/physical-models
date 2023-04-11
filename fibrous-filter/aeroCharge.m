function boltzmannStruct = aeroCharge(particleSizes)
%AEROCHARGE Return charge and neutral fraction of aerosols by radii.
% particleSizes - array of particle diameters in meters
% temperature   - ambient temperature in Kelvin
% Calculate the Boltzmann equilibium charge average and neutral fraction
% for aerosols of a certaian radius. The shape is assumed to be a sphere.

% elementary charge
e = 1.602176634*10^(-19);
% permittivity free space | farads per meter
eps0 = 8.8541878128 * 10^(-12);
% Boltzmann constant | Joules per mole - Kelvin
kB = 1.380649 * 10^(-23);
% temperature | Kelvins
T = 298.15;

chargeMultiple = (1:1000)';

% convert particle size vector to row if needed
if iscolumn(particleSizes)
    particleSizes = particleSizes';
end

aerosolRadius = particleSizes * 0.5;

%% electrostatic energy of a sphere of charge
electroPotent = (3/5) * 1/(4*pi*eps0) * (chargeMultiple * e).^2 ./ aerosolRadius;

exponent1 = exp(-electroPotent/(kB*T));

exponent2 = 2 * e * chargeMultiple .* exponent1;

% partition function
partitionFunction = 2 * sum(exponent1)';

%% magnitude and average charge for charged aerosols by radii
chargeMagnitude = sum(exponent2)'./partitionFunction;

boltzmannStruct.numberCharge = chargeMagnitude ./ e;
boltzmannStruct.neutralFraction = 1 ./ (1 + partitionFunction);
boltzmannStruct.chargedFraction = 1 - boltzmannStruct.neutralFraction;

end

