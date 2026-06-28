function [rateMap, spacingVec, orientationVec, phaseVec, gridField_locs] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, ...
    gridSpacing_variability, arenaSize, gridOrientation, gridOrientation_variability, nNeurons, arenaResolution, poissonFlag, oneD_flag)
%-------------------------------------------------------------------------%
%   This function generate the synthetic (idealized) rate maps via the sum
%   of three two-dimensional sinusoid gratings (as Solstad et al., 2006).
%   Grid orientations and spacings are sampled randomly, according to the
%   variability defined by gridSpacing_variability and
%   gridOrientation_variability. Grid phase is uniformly sampled. 
%
%   Written by WTR 06/05/2023 // Last updated by WTR 06/03/2026
%-------------------------------------------------------------------------%
%% Globals
arenaX = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
arenaY = (-arenaSize(2)/2):arenaResolution:(arenaSize(2)/2);

%% Ratemaps
rateMap = zeros(length(arenaX), length(arenaY), nNeurons); 
spacingVec = zeros(1, nNeurons);
orientationVec = zeros(1, nNeurons);
phaseVec = zeros(2, nNeurons);
n_modules = length(gridOrientation);

for nn = 1:nNeurons

    module_id = ceil(nn * n_modules / nNeurons);

    orientation = normrnd(gridOrientation(module_id), gridOrientation_variability) * pi / 180;
    spacing = normrnd(gridSpacing(module_id), gridSpacing_variability); 
    if oneD_flag ~= 1
        if nn > 1
            phase = rand(1, 2) .* arenaSize; 
        else
            phase = zeros(1, 2);
        end
    else
       phase = zeros(1, 2);
    end

    orientationVec(nn) = orientation;
    spacingVec(nn) = spacing;
    phaseVec(:, nn) = phase;
    
    for xx = 1:length(arenaX)
        for yy = 1:length(arenaY)
            k = 4 * pi / (sqrt(3) * spacing);
            k1 = k / sqrt(2) * [cos(orientation + pi / 12) + sin(orientation + pi / 12), cos(orientation + pi / 12) - sin(orientation + pi / 12)]; 
            k2 = k / sqrt(2) * [cos(orientation + 5 * pi / 12) + sin(orientation + 5 * pi / 12), cos(orientation + 5 * pi / 12) - sin(orientation + 5 * pi / 12)];
            k3 = k / sqrt(2) * [cos(orientation + 3 * pi / 4) + sin(orientation + 3 * pi / 4), cos(orientation + 3 * pi / 4) - sin(orientation + 3 * pi / 4)]; 
            waveVectors = [k1; k2; k3];
            rateMap(xx, yy, nn) = gridFiring_max * 2/3 * (1/3 * sum(cos(waveVectors * [arenaX(xx) + phase(1); arenaY(yy) + phase(2)])) + 1/2);
        end
    end

    if nn == 1
        [gridField_locs_x, gridField_locs_y] = find(abs(rateMap(:, :, nn) - gridFiring_max) < 0.10); % set a little arbitrarily but seems to pick up all the grid fields for the 2D case
        gridField_locs = [gridField_locs_x, gridField_locs_y];
    end

    if poissonFlag == 1
        rateMap(:, :, nn) =  poissrnd(rateMap(:, :, nn));
    end

end

