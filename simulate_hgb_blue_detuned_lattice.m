function result = simulate_hgb_blue_detuned_lattice()
%SIMULATE_HGB_BLUE_DETUNED_LATTICE Counter-propagating hollow-Gaussian beams.
% Two HGB waists are at z=-5 cm and z=+5 cm and propagate toward each other.
% SI units are used internally. Intensity is W/m^2 and potential is joules.

%% Physical and numerical configuration
cfg.laserWavelength = 760e-9;          % blue-detuned laser wavelength [m]
cfg.transitionWavelength = 780.241e-9; % Rb-87 D2 transition wavelength [m]
cfg.linewidthHz = 6.065e6;             % natural linewidth Gamma/(2*pi) [Hz]
cfg.totalPower = 0.20;                 % fixed power of both beams together [W]
cfg.sourceZ = [-0.05, 0.05];           % HGB waist positions [m]
cfg.relativePhase = 0;                 % phase of right beam relative to left [rad]

% Bounded grid search. These bounds make "maximum trap depth" well posed.
cfg.orderCandidates = 1:6;
cfg.waistCandidates = (0.20:0.05:1.00)*1e-3; % waist radius [m]
cfg.powerFractionCandidates = 0.20:0.05:0.80; % fraction in left beam

cfg.transverseRange = 2.5e-3;          % x/y half-width [m]
cfg.transverseSamples = 501;
cfg.macroZSamples = 801;               % envelope samples over source-to-source range
cfg.latticePeriodsShown = 8;           % central standing-wave periods displayed
cfg.latticeZSamples = 1601;
cfg.outputDirectory = 'results_hgb';
cfg.makeFigures = true;

validateConfiguration(cfg);
if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end

constants = atomicConstants(cfg);
x = linspace(-cfg.transverseRange, cfg.transverseRange, cfg.transverseSamples);

%% Optimize HGB order, waist and power division at fixed total power
radialOptimizationAxis = linspace(0, cfg.transverseRange, 2001);
nTrials = numel(cfg.orderCandidates)*numel(cfg.waistCandidates)* ...
    numel(cfg.powerFractionCandidates);
history = zeros(nTrials, 5); % order, waist, left fraction, barrier J, barrier uK
trial = 0;
best.barrier = -Inf;

for order = cfg.orderCandidates
    for waist = cfg.waistCandidates
        for leftFraction = cfg.powerFractionCandidates
            trial = trial + 1;
            intensity = interferenceIntensity(radialOptimizationAxis, 0, order, ...
                waist, leftFraction, cfg, constants.k);
            potential = constants.dipoleCoefficient*intensity;
            barrier = max(potential) - potential(1);
            history(trial, :) = [order, waist, leftFraction, barrier, ...
                barrier/constants.kB*1e6];
            if barrier > best.barrier
                best.order = order;
                best.waist = waist;
                best.leftPowerFraction = leftFraction;
                best.barrier = barrier;
            end
        end
    end
end

best.leftPower = cfg.totalPower*best.leftPowerFraction;
best.rightPower = cfg.totalPower-best.leftPower;
best.barrierMicroK = best.barrier/constants.kB*1e6;

%% z=0 transverse intensity and blue-detuned optical potential
[X, Y] = meshgrid(x, x);
R = hypot(X, Y);
intensityXY = interferenceIntensity(R, 0, best.order, best.waist, ...
    best.leftPowerFraction, cfg, constants.k);
potentialXY = constants.dipoleCoefficient*intensityXY;
potentialXYMicroK = potentialXY/constants.kB*1e6;

%% Propagation: macroscopic envelope plus a wavelength-resolved lattice view
zMacro = linspace(cfg.sourceZ(1), cfg.sourceZ(2), cfg.macroZSamples);
[XMacro, ZMacro] = meshgrid(x, zMacro);
[leftMacro, rightMacro] = hgbPair(XMacro, ZMacro, best.order, best.waist, ...
    best.leftPowerFraction, cfg, constants.k);
% A 10-cm plot cannot numerically resolve a sub-micron standing-wave period.
% Plot its phase-averaged envelope here and show the exact interference below.
intensityEnvelopeXZ = abs(leftMacro).^2 + abs(rightMacro).^2;

latticePeriod = cfg.laserWavelength/2;
zZoom = linspace(-cfg.latticePeriodsShown*latticePeriod/2, ...
    cfg.latticePeriodsShown*latticePeriod/2, cfg.latticeZSamples);
[XZoom, ZZoom] = meshgrid(x, zZoom);
[leftZoom, rightZoom] = hgbPair(XZoom, ZZoom, best.order, best.waist, ...
    best.leftPowerFraction, cfg, constants.k);
fieldLatticeXZ = leftZoom + exp(1i*cfg.relativePhase)*rightZoom;
intensityLatticeXZ = abs(fieldLatticeXZ).^2;
potentialLatticeXZMicroK = constants.dipoleCoefficient*intensityLatticeXZ/ ...
    constants.kB*1e6;

%% Save machine-readable results
historyTable = array2table(history, 'VariableNames', ...
    {'order', 'waist_m', 'left_power_fraction', 'barrier_J', 'barrier_microK'});
writetable(historyTable, fullfile(cfg.outputDirectory, 'optimization_history.csv'));
bestTable = table(best.order, best.waist, best.leftPowerFraction, best.leftPower, ...
    best.rightPower, best.barrier, best.barrierMicroK, ...
    'VariableNames', {'order', 'waist_m', 'left_power_fraction', 'left_power_W', ...
    'right_power_W', 'barrier_J', 'barrier_microK'});
writetable(bestTable, fullfile(cfg.outputDirectory, 'optimized_parameters.csv'));

result.config = cfg;
result.constants = constants;
result.optimized = best;
result.optimizationHistory = historyTable;
result.x = x;
result.zMacro = zMacro;
result.intensityEnvelopeXZ = intensityEnvelopeXZ;
result.zZoom = zZoom;
result.fieldLatticeXZ = fieldLatticeXZ;
result.intensityLatticeXZ = intensityLatticeXZ;
result.potentialLatticeXZMicroK = potentialLatticeXZMicroK;
result.intensityXY = intensityXY;
result.potentialXY = potentialXY;
result.potentialXYMicroK = potentialXYMicroK;
save(fullfile(cfg.outputDirectory, 'hgb_lattice_result.mat'), 'result', '-v7.3');

if cfg.makeFigures
    saveResultFigures(result);
end

fprintf('Optimized HGB: m=%d, w0=%.3f mm, Pleft/Ptotal=%.2f\n', ...
    best.order, best.waist*1e3, best.leftPowerFraction);
fprintf('Blue-detuned transverse barrier at z=0: %.3f microkelvin\n', ...
    best.barrierMicroK);
fprintf('Results saved in %s\n', cfg.outputDirectory);
end

function validateConfiguration(cfg)
assert(cfg.sourceZ(1) < 0 && cfg.sourceZ(2) > 0, ...
    'The two source planes must lie on opposite sides of z=0.');
assert(abs(cfg.sourceZ(1)+cfg.sourceZ(2)) < 10*eps(max(abs(cfg.sourceZ))), ...
    'This implementation expects source positions symmetric about z=0.');
assert(cfg.laserWavelength < cfg.transitionWavelength, ...
    'laserWavelength must be blue-detuned (shorter than the transition wavelength).');
assert(cfg.totalPower > 0, 'totalPower must be positive.');
assert(all(cfg.orderCandidates >= 1 & mod(cfg.orderCandidates, 1) == 0), ...
    'HGB orders must be positive integers.');
assert(all(cfg.waistCandidates > 0), 'All waist candidates must be positive.');
assert(all(cfg.powerFractionCandidates > 0 & cfg.powerFractionCandidates < 1), ...
    'Power fractions must lie strictly between zero and one.');
end

function constants = atomicConstants(cfg)
constants.c = 299792458;
constants.kB = 1.380649e-23;
constants.k = 2*pi/cfg.laserWavelength;
constants.omega0 = 2*pi*constants.c/cfg.transitionWavelength;
constants.omegaLaser = 2*pi*constants.c/cfg.laserWavelength;
constants.detuning = constants.omegaLaser-constants.omega0;
constants.Gamma = 2*pi*cfg.linewidthHz;
constants.dipoleCoefficient = 3*pi*constants.c^2*constants.Gamma/ ...
    (2*constants.omega0^3*constants.detuning);
assert(constants.detuning/constants.Gamma > 100, ...
    'The selected laser is not sufficiently far detuned for this model.');
end

function intensity = interferenceIntensity(radius, z, order, waist, leftFraction, cfg, k)
[leftField, rightField] = hgbPair(radius, z, order, waist, leftFraction, cfg, k);
totalField = leftField + exp(1i*cfg.relativePhase)*rightField;
intensity = abs(totalField).^2;
end

function [leftField, rightField] = hgbPair(radius, z, order, waist, leftFraction, cfg, k)
leftDistance = z-cfg.sourceZ(1);
rightDistance = cfg.sourceZ(2)-z;
leftField = hgbField(radius, leftDistance, order, waist, ...
    cfg.totalPower*leftFraction, cfg.laserWavelength, k);
rightField = hgbField(radius, rightDistance, order, waist, ...
    cfg.totalPower*(1-leftFraction), cfg.laserWavelength, k);
end

function field = hgbField(radius, distance, order, waist0, power, wavelength, k)
% Paraxial HGB model, normalized so integral(|E|^2 dA) equals beam power.
rayleighRange = pi*waist0^2/wavelength;
waist = waist0*sqrt(1+(distance/rayleighRange).^2);
gouy = atan(distance/rayleighRange);
inverseCurvature = distance./(distance.^2+rayleighRange^2);
normalization = sqrt(power*2^(2*order+1)./ ...
    (pi*waist.^2*gamma(2*order+1)));
radialFactor = (radius.^2./waist.^2).^order;
amplitude = normalization.*radialFactor.*exp(-radius.^2./waist.^2);
phase = k*distance + 0.5*k*radius.^2.*inverseCurvature - (2*order+1)*gouy;
field = amplitude.*exp(1i*phase);
end

function saveResultFigures(result)
cfg = result.config;
xMm = result.x*1e3;

figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
subplot(1, 2, 1);
imagesc(xMm, result.zMacro*1e2, result.intensityEnvelopeXZ);
axis xy; colorbar; xlabel('x [mm]'); ylabel('z [cm]');
title('Phase-averaged HGB propagation envelope [W/m^2]');
subplot(1, 2, 2);
imagesc(xMm, result.zZoom*1e6, result.intensityLatticeXZ);
axis xy; colorbar; xlabel('x [mm]'); ylabel('z [um]');
title('Resolved counter-propagating interference lattice [W/m^2]');
saveas(gcf, fullfile(cfg.outputDirectory, 'hgb_propagation.png'));

figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
subplot(1, 2, 1);
imagesc(xMm, xMm, result.intensityXY); axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('z=0 intensity [W/m^2]');
subplot(1, 2, 2);
imagesc(xMm, xMm, result.potentialXYMicroK); axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('z=0 blue-detuned potential [microkelvin]');
saveas(gcf, fullfile(cfg.outputDirectory, 'xy_intensity_and_potential.png'));

figure('Color', 'w');
scatter3(result.optimizationHistory.waist_m*1e3, ...
    result.optimizationHistory.order, result.optimizationHistory.barrier_microK, ...
    28, result.optimizationHistory.left_power_fraction, 'filled');
xlabel('w0 [mm]'); ylabel('HGB order m'); zlabel('Barrier [microkelvin]');
cb = colorbar; cb.Label.String = 'Left-beam power fraction'; grid on;
title('Bounded HGB parameter optimization');
saveas(gcf, fullfile(cfg.outputDirectory, 'optimization.png'));
close all;
end
