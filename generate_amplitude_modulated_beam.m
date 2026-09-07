%% Pure-amplitude generation of an order-n hollow Gaussian beam
% The incident Gaussian has a flat spatial phase.  A real, non-negative
% transmission mask therefore changes only its amplitude, not its phase.

clear; close all; clc;

%% Parameters
lambda = 632.8e-9;       % wavelength [m] (used only for documentation here)
wG     = 1.5e-3;         % incident Gaussian waist [m]
w0     = 1.0e-3;         % target beam waist [m]
n      = 2;               % non-negative integer beam order
E0     = 1;               % incident on-axis complex amplitude
L      = 4.0e-3;          % half-width of the simulated/masked aperture [m]
N      = 1001;            % number of samples along each transverse axis
alpha  = -1.0e-39;        % real dynamic polarizability Re[alpha(omega)] [C*m^2/V]

assert(n >= 0 && n == floor(n), 'n must be a non-negative integer.');
assert(wG > 0 && w0 > 0 && L > 0 && N >= 2, ...
    'Beam waists, aperture size, and grid size must be positive.');

%% Coordinate grid and input Gaussian
x = linspace(-L, L, N);
[X, Y] = meshgrid(x, x);
r2 = X.^2 + Y.^2;

% A constant phase could be included by replacing E0 with E0*exp(1i*phi0).
EG = E0 .* exp(-r2 ./ wG^2);

%% Real, non-negative amplitude transmission
% The unscaled mask is the ratio of the desired radial amplitude and EG.
% C is chosen on the finite aperture so that 0 <= Tn <= 1, as required for
% a passive amplitude mask.  This scaling changes power but not beam shape.
Traw = (r2 ./ w0^2).^n .* exp(-r2 ./ w0^2 + r2 ./ wG^2);
Tmax = max(Traw(:));
if Tmax == 0
    C = 1;                % defensive fallback for a degenerate grid
else
    C = 1 / Tmax;
end
Tn = C .* Traw;

% Pure-amplitude modulation: no extra exp(1i*phi) factor is applied.
Eout = Tn .* EG;

%% Numerical verification
Etarget = C .* E0 .* (r2 ./ w0^2).^n .* exp(-r2 ./ w0^2);
relativeError = norm(Eout(:) - Etarget(:)) / max(norm(Etarget(:)), eps);

fprintf('lambda = %.1f nm, order n = %d\n', lambda * 1e9, n);
fprintf('Mask range: [%.6g, %.6g]\n', min(Tn(:)), max(Tn(:)));
fprintf('Relative field error: %.3e\n', relativeError);

%% Optical dipole-trap potential
% Here Eout is the peak complex electric-field amplitude, defined through
% E(t) = Re[Eout*exp(-1i*omega*t)].  In vacuum,
%
%        I(x,y) = (epsilon0*c/2)*|Eout(x,y)|^2,
%        U(x,y) = -(1/4)*Re[alpha(omega)]*|Eout(x,y)|^2
%               = -Re[alpha(omega)]*I(x,y)/(2*epsilon0*c).
%
% alpha > 0 (usually red detuning) attracts particles to high intensity;
% alpha < 0 (usually blue detuning) attracts them to the dark region.
epsilon0 = 8.8541878128e-12;    % vacuum permittivity [F/m]
c0       = 299792458;           % speed of light in vacuum [m/s]
kB       = 1.380649e-23;        % Boltzmann constant [J/K]

Iout = 0.5 .* epsilon0 .* c0 .* abs(Eout).^2;  % output intensity [W/m^2]
Udip = -0.25 .* alpha .* abs(Eout).^2;         % optical potential [J]
UuK  = Udip ./ kB .* 1e6;                      % potential in microkelvin

fprintf('Optical-potential range: [%.6g, %.6g] J\n', ...
    min(Udip(:)), max(Udip(:)));
fprintf('Equivalent U/k_B range: [%.6g, %.6g] microkelvin\n', ...
    min(UuK(:)), max(UuK(:)));

%% Visualize input, mask, output intensity, and phase difference
extentMM = x * 1e3;
figure('Color', 'w', 'Name', 'Pure-amplitude beam generation');

subplot(2, 3, 1);
imagesc(extentMM, extentMM, abs(EG).^2);
axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('Input Gaussian intensity');

subplot(2, 3, 2);
imagesc(extentMM, extentMM, Tn);
axis image xy; colorbar; caxis([0, 1]);
xlabel('x [mm]'); ylabel('y [mm]'); title('Amplitude transmission T_n');

subplot(2, 3, 3);
imagesc(extentMM, extentMM, abs(Eout).^2);
axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('Output intensity');

% Compare phase only where the output amplitude is nonzero; phase at a
% field zero (including the center for n > 0) is mathematically undefined.
phaseDifference = nan(size(Eout));
valid = abs(Eout) > max(abs(Eout(:))) * 1e-12;
phaseDifference(valid) = angle(Eout(valid) .* conj(EG(valid)));

subplot(2, 3, 4);
imagesc(extentMM, extentMM, phaseDifference);
axis image xy; colorbar; caxis([-pi, pi]);
xlabel('x [mm]'); ylabel('y [mm]');
title('\phi_{out} - \phi_{Gaussian} [rad]');

subplot(2, 3, 5);
imagesc(extentMM, extentMM, UuK);
axis image xy; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title('U_{dip}/k_B [\muK]');

% A center-line cut makes the radial trapping structure easier to inspect.
subplot(2, 3, 6);
centerIndex = ceil(N / 2);
plot(extentMM, UuK(centerIndex, :), 'LineWidth', 1.5);
grid on;
xlabel('x [mm]'); ylabel('U_{dip}/k_B [\muK]');
title('Potential along y = 0');

colormap(parula);
