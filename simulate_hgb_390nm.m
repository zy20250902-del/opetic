function result = simulate_hgb_390nm(z1, z2)
%SIMULATE_HGB_390NM Counter-propagating 390-nm hollow-Gaussian beams.
% RESULT = SIMULATE_HGB_390NM(Z1,Z2) places the +z beam waist at Z1 and the
% -z beam waist at Z2 (metres). The observation plane is their midpoint.

if nargin < 1, z1 = -0.05; end
if nargin < 2, z2 = 0.05; end

%% User-adjustable optical parameters (SI units)
cfg.wavelength = 390e-9;
cfg.z1 = z1;
cfg.z2 = z2;
cfg.observationZ = (z1+z2)/2;
cfg.hgbOrder = 1;
cfg.waist1 = 0.40e-3;
cfg.waist2 = 0.40e-3;
cfg.power1 = 0.10;
cfg.power2 = 0.10;
cfg.relativePhase = 0;
% Replace this illustrative value with the real dynamic polarizability of
% the particle at 390 nm. Negative polarizability describes blue detuning.
cfg.dynamicPolarizabilitySI = -1.0e-39; % [C*m^2/V]

cfg.transverseHalfWidth = 2.0e-3;
cfg.transverseSamples = 501;
cfg.macroZSamples = 801;
cfg.latticePeriodsShown = 10;
cfg.latticeZSamples = 2001;
cfg.outputDirectory = 'results_hgb_390nm';
cfg.makeFigures = true;

validateConfiguration(cfg);
if ~isfolder(cfg.outputDirectory), mkdir(cfg.outputDirectory); end

c = 299792458;
epsilon0 = 8.8541878128e-12;
k = 2*pi/cfg.wavelength;
% U = C_U*I and F_grad = -gradient(U).
potentialCoefficient = -cfg.dynamicPolarizabilitySI/(2*epsilon0*c);
x = linspace(-cfg.transverseHalfWidth, cfg.transverseHalfWidth, ...
    cfg.transverseSamples);
dx = x(2)-x(1);

%% Macroscopic xz propagation envelope
zMacro = linspace(cfg.z1, cfg.z2, cfg.macroZSamples);
[XMacro, ZMacro] = meshgrid(x, zMacro);
[field1Macro, field2Macro] = hgbPair(XMacro, ZMacro, cfg, k);
% A centimetre-scale grid cannot resolve a 195-nm lattice period.
intensityEnvelopeXZ = abs(field1Macro).^2+abs(field2Macro).^2;

%% Wavelength-resolved xz interference and z-directed gradient force
latticePeriod = cfg.wavelength/2;
zHalfWidth = cfg.latticePeriodsShown*latticePeriod/2;
zFine = linspace(cfg.observationZ-zHalfWidth, ...
    cfg.observationZ+zHalfWidth, cfg.latticeZSamples);
dz = zFine(2)-zFine(1);
[XFine, ZFine] = meshgrid(x, zFine);
[field1Fine, field2Fine] = hgbPair(XFine, ZFine, cfg, k);
fieldXZ = field1Fine+exp(1i*cfg.relativePhase)*field2Fine;
intensityXZ = abs(fieldXZ).^2;
potentialXZ = potentialCoefficient*intensityXZ;
[dPotentialDx, dPotentialDz] = gradient(potentialXZ, dx, dz);
forceXZ_X = -dPotentialDx;
forceXZ_Z = -dPotentialDz;

%% xy intensity and transverse gradient force at the observation plane
[X, Y] = meshgrid(x, x);
radius = hypot(X, Y);
[field1XY, field2XY] = hgbPair(radius, cfg.observationZ, cfg, k);
fieldXY = field1XY+exp(1i*cfg.relativePhase)*field2XY;
intensityXY = abs(fieldXY).^2;
potentialXY = potentialCoefficient*intensityXY;
[dPotentialDx, dPotentialDy] = gradient(potentialXY, dx, dx);
forceXY_X = -dPotentialDx;
forceXY_Y = -dPotentialDy;
forceXYMagnitude = hypot(forceXY_X, forceXY_Y);

% An ideal HGB is dark on its geometrical axis. Also extract the useful
% axial profile through the bright annulus, in addition to the true axis.
centerIndex = (cfg.transverseSamples+1)/2;
[~, positivePeakIndex] = max(intensityXY(centerIndex, centerIndex:end));
brightRingIndex = centerIndex+positivePeakIndex-1;
brightRingRadius = x(brightRingIndex);

%% Return and save all numerical data
result.config = cfg;
result.x = x;
result.zMacro = zMacro;
result.intensityEnvelopeXZ = intensityEnvelopeXZ;
result.zFine = zFine;
result.fieldXZ = fieldXZ;
result.intensityXZ = intensityXZ;
result.potentialXZ = potentialXZ;
result.forceXZ_X = forceXZ_X;
result.forceXZ_Z = forceXZ_Z;
result.intensityXY = intensityXY;
result.potentialXY = potentialXY;
result.forceXY_X = forceXY_X;
result.forceXY_Y = forceXY_Y;
result.forceXYMagnitude = forceXYMagnitude;
result.onAxisIntensity = intensityXZ(:, centerIndex);
result.onAxisForceZ = forceXZ_Z(:, centerIndex);
result.brightRingRadius = brightRingRadius;
result.brightRingIntensityZ = intensityXZ(:, brightRingIndex);
result.brightRingForceZ = forceXZ_Z(:, brightRingIndex);
result.potentialCoefficient = potentialCoefficient;
save(fullfile(cfg.outputDirectory, 'hgb_390nm_result.mat'), 'result', '-v7.3');

summaryTable = table(cfg.z1, cfg.z2, cfg.observationZ, cfg.wavelength, ...
    cfg.hgbOrder, cfg.waist1, cfg.waist2, cfg.power1, cfg.power2, ...
    cfg.dynamicPolarizabilitySI, max(intensityXY(:)), ...
    max(abs(forceXZ_Z(:))), max(forceXYMagnitude(:)), ...
    'VariableNames', {'z1_m', 'z2_m', 'observation_z_m', 'wavelength_m', ...
    'hgb_order', 'waist1_m', 'waist2_m', 'power1_W', 'power2_W', ...
    'polarizability_SI', 'max_xy_intensity_W_m2', 'max_abs_Fz_N', ...
    'max_transverse_force_N'});
writetable(summaryTable, fullfile(cfg.outputDirectory, 'simulation_summary.csv'));

if cfg.makeFigures, saveFigures(result); end
fprintf('390-nm HGB: z1=%.4f m, z2=%.4f m, xy plane z=%.4f m\n', ...
    cfg.z1, cfg.z2, cfg.observationZ);
fprintf('Maximum |Fz| = %.3e N; maximum transverse |F| = %.3e N\n', ...
    max(abs(forceXZ_Z(:))), max(forceXYMagnitude(:)));
end

function validateConfiguration(cfg)
assert(isfinite(cfg.z1) && isfinite(cfg.z2) && cfg.z1 < cfg.z2, ...
    'Require finite waist positions with z1 < z2.');
assert(cfg.wavelength > 0, 'wavelength must be positive.');
assert(cfg.hgbOrder >= 1 && mod(cfg.hgbOrder, 1) == 0, ...
    'hgbOrder must be a positive integer.');
assert(cfg.waist1 > 0 && cfg.waist2 > 0, 'Beam waists must be positive.');
assert(cfg.power1 > 0 && cfg.power2 > 0, 'Beam powers must be positive.');
assert(cfg.transverseSamples >= 3 && mod(cfg.transverseSamples, 2) == 1, ...
    'transverseSamples must be an odd integer of at least 3.');
assert(cfg.latticeZSamples >= 3, 'latticeZSamples must be at least 3.');
end

function [field1, field2] = hgbPair(radius, z, cfg, k)
distance1 = z-cfg.z1; % beam 1 propagates in +z
distance2 = cfg.z2-z; % beam 2 propagates in -z
field1 = hgbField(radius, distance1, cfg.hgbOrder, cfg.waist1, ...
    cfg.power1, cfg.wavelength, k);
field2 = hgbField(radius, distance2, cfg.hgbOrder, cfg.waist2, ...
    cfg.power2, cfg.wavelength, k);
end

function field = hgbField(radius, distance, order, waist0, power, wavelength, k)
% Power-normalized paraxial HGB: integral |field|^2 dA = power.
rayleighRange = pi*waist0^2/wavelength;
waist = waist0*sqrt(1+(distance/rayleighRange).^2);
gouyPhase = atan(distance/rayleighRange);
inverseCurvature = distance./(distance.^2+rayleighRange^2);
normalization = sqrt(power*2^(2*order+1)./ ...
    (pi*waist.^2*gamma(2*order+1)));
radialFactor = (radius.^2./waist.^2).^order;
amplitude = normalization.*radialFactor.*exp(-radius.^2./waist.^2);
phase = k*distance+0.5*k*radius.^2.*inverseCurvature- ...
    (2*order+1)*gouyPhase;
field = amplitude.*exp(1i*phase);
end

function saveFigures(result)
cfg = result.config;
xMm = result.x*1e3;
figure('Color', 'w', 'Position', [100, 100, 1500, 450]);
subplot(1, 3, 1);
imagesc(xMm, result.zMacro*1e2, result.intensityEnvelopeXZ);
axis xy; colorbar; xlabel('x [mm]'); ylabel('z [cm]');
title('Macroscopic propagation envelope [W/m^2]');
subplot(1, 3, 2);
imagesc(xMm, (result.zFine-cfg.observationZ)*1e6, result.intensityXZ);
axis xy; colorbar; xlabel('x [mm]'); ylabel('z-z_{obs} [um]');
title('Resolved z-direction interference [W/m^2]');
subplot(1, 3, 3);
imagesc(xMm, xMm, result.intensityXY); axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('xy-plane intensity [W/m^2]');
saveas(gcf, fullfile(cfg.outputDirectory, 'intensity_maps.png'));

figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
subplot(1, 2, 1);
imagesc(xMm, (result.zFine-cfg.observationZ)*1e6, result.forceXZ_Z*1e12);
axis xy; colorbar; xlabel('x [mm]'); ylabel('z-z_{obs} [um]');
title('z-directed gradient force F_z [pN]');
subplot(1, 2, 2);
imagesc(xMm, xMm, result.forceXYMagnitude*1e12); axis image xy; colorbar; hold on;
step = max(1, floor(numel(xMm)/25)); indices = 1:step:numel(xMm);
quiver(xMm(indices), xMm(indices), result.forceXY_X(indices, indices), ...
    result.forceXY_Y(indices, indices), 'w');
xlabel('x [mm]'); ylabel('y [mm]');
title('xy gradient-force magnitude and direction [pN]');
saveas(gcf, fullfile(cfg.outputDirectory, 'gradient_force_maps.png'));

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
subplot(2, 1, 1);
plot((result.zFine-cfg.observationZ)*1e6, result.brightRingIntensityZ);
ylabel('Intensity [W/m^2]'); grid on;
title(sprintf('Axial profile at bright ring r = %.3f mm', ...
    result.brightRingRadius*1e3));
subplot(2, 1, 2);
plot((result.zFine-cfg.observationZ)*1e6, result.brightRingForceZ*1e12);
xlabel('z-z_{obs} [um]'); ylabel('F_z [pN]'); grid on;
saveas(gcf, fullfile(cfg.outputDirectory, 'axial_profiles.png'));
close all;
end
