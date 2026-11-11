function gaussian_beam_obstruction(varargin)
%GAUSSIAN_BEAM_OBSTRUCTION Simulate a centrally blocked Gaussian beam.
%   GAUSSIAN_BEAM_OBSTRUCTION() propagates a 695 nm TEM00 Gaussian beam
%   along +z, places an opaque disk at z = 15 cm, and calculates the field
%   through z = 25 cm. Figures and numerical data are written to output/.
%
%   Name-value pairs can override Wavelength, Waist, ObstacleZ,
%   ObstacleRadius, FinalZ, Window, Samples, ZSamples, and OutputDirectory.

    parameters = parse_inputs(varargin{:});
    if parameters.FinalZ <= parameters.ObstacleZ
        error('FinalZ must be greater than ObstacleZ.');
    end
    if ~exist(parameters.OutputDirectory, 'dir')
        mkdir(parameters.OutputDirectory);
    end

    n = parameters.Samples;
    dx = parameters.Window / n;
    coordinate = (-n/2:n/2-1) * dx;
    [xGrid, yGrid] = meshgrid(coordinate, coordinate);

    % At z=0 the beam waist has a planar wavefront and unit peak field.
    initialField = exp(-(xGrid.^2 + yGrid.^2) / parameters.Waist^2);
    referenceIntensity = max(abs(initialField(:)).^2);

    fieldAtObstacle = fresnel_propagate(initialField, parameters.ObstacleZ, ...
        parameters.Wavelength, dx);
    opaqueDisk = (xGrid.^2 + yGrid.^2) < parameters.ObstacleRadius^2;
    fieldAfterObstacle = fieldAtObstacle;
    fieldAfterObstacle(opaqueDisk) = 0;
    finalField = fresnel_propagate(fieldAfterObstacle, ...
        parameters.FinalZ - parameters.ObstacleZ, parameters.Wavelength, dx);

    z = linspace(0, parameters.FinalZ, parameters.ZSamples);
    intensity = zeros(n, n, parameters.ZSamples, 'single');
    for index = 1:parameters.ZSamples
        if z(index) < parameters.ObstacleZ
            field = fresnel_propagate(initialField, z(index), ...
                parameters.Wavelength, dx);
        else
            field = fresnel_propagate(fieldAfterObstacle, ...
                z(index) - parameters.ObstacleZ, parameters.Wavelength, dx);
        end
        intensity(:, :, index) = single(abs(field).^2 / referenceIntensity);
    end

    save_xy_planes(parameters.OutputDirectory, coordinate, initialField, ...
        fieldAfterObstacle, finalField, referenceIntensity, parameters);
    save_longitudinal(parameters.OutputDirectory, coordinate, z, intensity, parameters);
    save_3d(parameters.OutputDirectory, coordinate, z, intensity, parameters);

    wavelength = parameters.Wavelength; %#ok<NASGU>
    waist = parameters.Waist; %#ok<NASGU>
    obstacleZ = parameters.ObstacleZ; %#ok<NASGU>
    obstacleRadius = parameters.ObstacleRadius; %#ok<NASGU>
    x = coordinate; %#ok<NASGU>
    y = coordinate; %#ok<NASGU>
    save(fullfile(parameters.OutputDirectory, 'intensity_data.mat'), ...
        'x', 'y', 'z', 'intensity', 'wavelength', 'waist', ...
        'obstacleZ', 'obstacleRadius', '-v7.3');
end


function parameters = parse_inputs(varargin)
    parser = inputParser;
    positiveScalar = @(value) isnumeric(value) && isscalar(value) && ...
        isfinite(value) && value > 0;
    integerScalar = @(value) positiveScalar(value) && value == round(value);
    addParameter(parser, 'Wavelength', 695e-9, positiveScalar);
    addParameter(parser, 'Waist', 1e-3, positiveScalar);
    addParameter(parser, 'ObstacleZ', 0.15, positiveScalar);
    addParameter(parser, 'ObstacleRadius', 0.4e-3, positiveScalar);
    addParameter(parser, 'FinalZ', 0.25, positiveScalar);
    addParameter(parser, 'Window', 8e-3, positiveScalar);
    addParameter(parser, 'Samples', 512, integerScalar);
    addParameter(parser, 'ZSamples', 101, integerScalar);
    addParameter(parser, 'OutputDirectory', 'output', ...
        @(value) ischar(value) || (isstring(value) && isscalar(value)));
    parse(parser, varargin{:});
    parameters = parser.Results;
    parameters.OutputDirectory = char(parameters.OutputDirectory);
end


function propagated = fresnel_propagate(field, distance, wavelength, spacing)
%FRESNEL_PROPAGATE Apply the paraxial transfer function in frequency space.
    n = size(field, 1);
    frequencies = ifftshift((-floor(n/2):ceil(n/2)-1) / (n * spacing));
    [fx, fy] = meshgrid(frequencies, frequencies);
    transferFunction = exp(-1i * pi * wavelength * distance .* (fx.^2 + fy.^2));
    propagated = ifft2(fft2(field) .* transferFunction);
end


function save_xy_planes(outputDirectory, coordinate, initialField, ...
        obstacleField, finalField, referenceIntensity, parameters)
    figureHandle = figure('Color', 'white', 'Position', [100, 100, 1400, 430]);
    fields = {initialField, obstacleField, finalField};
    positions = [0, parameters.ObstacleZ, parameters.FinalZ];
    for index = 1:3
        subplot(1, 3, index);
        imagesc(1e3 * coordinate, 1e3 * coordinate, ...
            abs(fields{index}).^2 / referenceIntensity, [0, 1]);
        axis image xy;
        xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('z = %.0f cm', 100 * positions(index)));
        colormap(gca, inferno_colormap(256));
        colorbar;
    end
    print(figureHandle, fullfile(outputDirectory, 'xy_intensity_planes.png'), ...
        '-dpng', '-r220');
    close(figureHandle);
end


function save_longitudinal(outputDirectory, coordinate, z, intensity, parameters)
    centerIndex = floor(numel(coordinate) / 2) + 1;
    xzIntensity = squeeze(intensity(centerIndex, :, :));
    figureHandle = figure('Color', 'white', 'Position', [100, 100, 900, 480]);
    imagesc(100 * z, 1e3 * coordinate, xzIntensity, [0, 1]);
    axis xy;
    xlabel('z (cm)'); ylabel('x (mm)');
    title('Longitudinal intensity (y = 0)');
    colormap(inferno_colormap(256)); colorbar;
    hold on;
    plot([100 100] * parameters.ObstacleZ, 1e3 * [coordinate(1), coordinate(end)], ...
        '--c', 'LineWidth', 1.2);
    hold off;
    print(figureHandle, fullfile(outputDirectory, 'longitudinal_xz_intensity.png'), ...
        '-dpng', '-r220');
    close(figureHandle);
end


function save_3d(outputDirectory, coordinate, z, intensity, parameters)
    xyStride = max(1, floor(numel(coordinate) / 64));
    zStride = max(1, floor(numel(z) / 60));
    reducedCoordinate = coordinate(1:xyStride:end);
    reducedZ = z(1:zStride:end);
    reducedIntensity = intensity(1:xyStride:end, 1:xyStride:end, 1:zStride:end);
    [xGrid, yGrid, zGrid] = meshgrid(reducedCoordinate, reducedCoordinate, reducedZ);
    visible = reducedIntensity >= 0.08;

    figureHandle = figure('Color', 'white', 'Position', [100, 100, 1000, 700]);
    scatter3(1e3 * xGrid(visible), 1e3 * yGrid(visible), 100 * zGrid(visible), ...
        5, double(reducedIntensity(visible)), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.1);
    colormap(inferno_colormap(256)); colorbar; caxis([0, 1]);
    hold on;
    angle = linspace(0, 2*pi, 150);
    plot3(1e3 * parameters.ObstacleRadius * cos(angle), ...
        1e3 * parameters.ObstacleRadius * sin(angle), ...
        100 * parameters.ObstacleZ * ones(size(angle)), 'c-', 'LineWidth', 2);
    hold off;
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (cm)');
    title('3D normalized intensity (I/I_0 >= 0.08)');
    axis tight; grid on; view(38, 24);
    print(figureHandle, fullfile(outputDirectory, 'xyz_intensity_3d.png'), ...
        '-dpng', '-r220');
    close(figureHandle);
end


function colors = inferno_colormap(numberOfColors)
%INFERNO_COLORMAP Toolbox-free dark-to-bright perceptual approximation.
    anchors = [0.0015 0.0005 0.0139; 0.2201 0.0160 0.3180; ...
        0.5783 0.1480 0.4044; 0.8650 0.3168 0.2261; ...
        0.9876 0.6453 0.0399; 0.9884 0.9984 0.6449];
    anchorPositions = linspace(0, 1, size(anchors, 1));
    colors = interp1(anchorPositions, anchors, linspace(0, 1, numberOfColors), 'pchip');
    colors = min(max(colors, 0), 1);
end
