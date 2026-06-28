%% GRID CELL DISTANCE PIPELINE
%-------------------------------------------------------------------------%
%   This script generates a population of grid cells, with varying
%   properties, and computes the correlation as a function of distance.
%   Based on the simulation code from Redman et al. eLife (2024). 
%
%   Written by WTR 01/13/2025 // Last updated by WTR 06/02/2026
%-------------------------------------------------------------------------%
%close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals      
gridSpacing = [0.85, sqrt(2) * 0.85, 0.85 / sqrt(2)];          % mean spacing of grid cells - in meters    
gridOrientation = [0, 0, 0];                     % mean orientation of grid cells - in degrees
arenaSize = [20, 20];                             % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = [128];                                % size of populations of grid cells to be tested
gridFiring_max = 15;                             % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = 1;                                      % standard deviation used for sampling grid orientation
spacingStd = 0.05;                               % standard deviation used for sampling grid spacing
arenaResolution = 0.05;                           % size of each bin used for discretizing space
nSamples = 25;                                   % number of grid populations to sample from
poissonFlag = 1;                                 % flag to sample from a Poisson distribution 

saveFlag = 0;                                    % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of size of grid population on decoding accuracy
corr0Variable = zeros(length(nNeurons), ceil(arenaSize(1)/arenaResolution + 1), nSamples);
corr90Variable = zeros(length(nNeurons), ceil(arenaSize(2)/arenaResolution + 1), nSamples);
corr45Variable = zeros(length(nNeurons), ceil(arenaSize(2)/arenaResolution + 1), nSamples);
corr0Fixed = zeros(length(nNeurons), ceil(arenaSize(1)/arenaResolution + 1), nSamples);
corr90Fixed = zeros(length(nNeurons), ceil(arenaSize(2)/arenaResolution + 1), nSamples);
corr0MultiFixed = zeros(length(nNeurons), ceil(arenaSize(1)/arenaResolution + 1), nSamples);
corr90MultiFixed = zeros(length(nNeurons), ceil(arenaSize(2)/arenaResolution + 1), nSamples);
corr_allVariable = cell(length(nNeurons), nSamples);
dist_allVariable = cell(length(nNeurons), nSamples);

for ss = 1:nSamples
    ss
    for nn = 1:length(nNeurons)
        [rateMap_Fixed, spacing_Fixed, orientation_Fixed, phase_Fixed, gridField_locsFixed] = ...
            RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing(1), 0, arenaSize, gridOrientation(1), 0, nNeurons(nn), arenaResolution, poissonFlag, oneD_flag);
        [corr0Fixed(nn, :, ss), corr90Fixed(nn, :, ss), ~, ~] = DISTANCE_CORR(rateMap_Fixed, gridField_locsFixed);

        [rateMap_Variable, spacing_Variable, orientation_Variable, phase_Variable, gridField_locsVariable] = ...
            RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing(1), spacingStd, arenaSize, gridOrientation(1), oriStd, nNeurons(nn), arenaResolution, poissonFlag, oneD_flag);
%         if ss == 1 && nn == 1
%              corr_allVariable = zeros(length(nNeurons), length(gridField_locsVariable), nSamples);
%              dist_allVariable = zeros(length(nNeurons), length(gridField_locsVariable), nSamples);
%         end
        [corr0Variable(nn, :, ss), corr90Variable(nn, :, ss), corr45Variable(nn, :, ss), corr_allVariable{nn, ss}, dist_allVariable{nn, ss}] = DISTANCE_CORR(rateMap_Variable, gridField_locsVariable);

        [rateMap_MultiFixed, spacing_MultiFixed, orientation_MultiFixed, phase_MultiFixed, gridField_locsMultiFixed] = ...
            RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, 0, arenaSize, gridOrientation, 0, nNeurons(nn), arenaResolution, poissonFlag, oneD_flag);
        [corr0MultiFixed(nn, :, ss), corr90MultiFixed(nn, :, ss), ~] = DISTANCE_CORR(rateMap_MultiFixed, gridField_locsMultiFixed);

    end
end

corr0Variable = median(corr0Variable, 3);
corr90Variable = median(corr90Variable, 3);
corr45Variable = median(corr45Variable, 3);
corr0Fixed = median(corr0Fixed, 3);
corr90Fixed = median(corr90Fixed, 3);
corr0MultiFixed = median(corr0MultiFixed, 3);
corr90MultiFixed = median(corr90MultiFixed, 3);

%% Saving data
if saveFlag == 1
    save(strcat(savePath, '/fixed_correlation_0_vs_distance_', num2str(gridSpacing), '.mat'), "corr0Fixed");
    save(strcat(savePath, '/fixed_correlation_90_vs_distance_', num2str(gridSpacing), '.mat'), "corr90Fixed");
    save(strcat(savePath, '/variable_correlation_0_vs_distance_', num2str(gridSpacing), '.mat'), "corr0Variable");
    save(strcat(savePath, '/variable_correlation_90_vs_distance_', num2str(gridSpacing), '.mat'), "corr90Variable");
    save(strcat(savePath, '/variable_correlation_all_vs_distance_', num2str(gridSpacing), '.mat'), "corr_allVariable");
end

%% Plotting
if oneD_flag
    figure
    imagesc(squeeze(rateMap_Fixed((arenaSize(1) / arenaResolution / 2):end, :, :))');
    title('One module: No variability')

    figure
    imagesc(squeeze(rateMap_MultiFixed((arenaSize(1) / arenaResolution / 2):end, :, :))');
    title('Three modules: No variability')

    figure
    imagesc(squeeze(rateMap_Variable((arenaSize(1) / arenaResolution / 2):end, :, :))');
    title('One module: Variability')
end

distance0 = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
distance90 = (-arenaSize(2)/2):arenaResolution:(arenaSize(2)/2);
bins_per_field0 = round(gridSpacing(1) / arenaResolution);
bins_per_field90 = sqrt(3) * gridSpacing(1) / arenaResolution;
ids0 = 1:size(rateMap_Fixed, 1);
ids90 = 1:size(rateMap_Fixed, 2);
id_fields0 = [flip(ids0(ceil(size(rateMap_Fixed, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Fixed, 1)/2 + bins_per_field0):bins_per_field0:end)];
id_fields90 =  [flip(ids90(ceil(size(rateMap_Fixed, 2)/2):(-bins_per_field90):1)), ids90(ceil(size(rateMap_Fixed, 2)/2 + bins_per_field90):bins_per_field90:end)];

% One module fixed properties
figure; 
legend_string = [];
if oneD_flag ~= 1
    subplot(1, 2, 1)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0Fixed(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k')    
        else
            plot(distance0(id_fields0), corr0Fixed(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('No variability 0 degree correlation')

    subplot(1, 2, 2)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance90(id_fields90), corr90Fixed(nn, id_fields90), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance90(id_fields90), corr90Fixed(nn, id_fields90), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(2)/2 + 0.2, 0.5, 1.1])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('No variability 90 degree correlation')

else
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0Fixed(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance0(id_fields0), corr0Fixed(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('No variability 0 degree correlation')
end

if saveFlag == 1
    savefig(strcat(savePath, 'One_module_fixed_properties.fig'));
end

% Multiple modules fixed properties
figure;
legend_string = [];
if oneD_flag ~= 1
    subplot(1, 2, 1)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0MultiFixed(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance0(id_fields0), corr0MultiFixed(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('Multi module no variability 0 degree correlation')

    subplot(1, 2, 2)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance90(id_fields90), corr90MultiFixed(nn, id_fields90), 'ko','MarkerFaceColor', 'k');
        else
            plot(distance90(id_fields90), corr90MultiFixed(nn, id_fields90), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(2)/2 + 0.2, 0.5, 1.05])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('Multi module no variability 90 degree correlation')

else
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0MultiFixed(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance0(id_fields0), corr0MultiFixed(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    legend(legend_string)
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    xlabel('Distance')
    ylabel('Correlation')
    title('Multi module no variability 0 degree correlation')
end

if saveFlag == 1
    savefig(strcat(savePath, 'Three_modules_fixed_properties.fig'));
end

% One module variable properties
figure; 
legend_string = [];
if oneD_flag ~= 1
    subplot(1, 2, 1)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0Variable(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance0(id_fields0), corr0Variable(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    legend(legend_string)
    xlabel('Distance')
    ylabel('Correlation')
    title('Variability 0 degree correlation')

    subplot(1, 2, 2)
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance90(id_fields90), corr90Variable(nn, id_fields90), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance90(id_fields90), corr90Variable(nn, id_fields90), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    axis([-0.2, arenaSize(2)/2 + 0.2, 0.5, 1.05])
    % axis square
    legend(legend_string)
    xlabel('Distance')
    ylabel('Correlation')
    title('Variability 90 degree correlation')

else
    for nn = 1:length(nNeurons)
        legend_string = [legend_string, {strcat('N = ', num2str(nNeurons(nn)))}];
        if length(nNeurons) == 1
            plot(distance0(id_fields0), corr0Variable(nn, id_fields0), 'ko', 'MarkerFaceColor', 'k');
        else
            plot(distance0(id_fields0), corr0Variable(nn, id_fields0), 'o', 'MarkerFaceColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)], 'MarkerEdgeColor', [(nn - 1)/(length(nNeurons) - 1), 0, (length(nNeurons) - nn)/(length(nNeurons) - 1)]); hold on
        end
    end
    axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
    % axis square
    legend(legend_string)
    xlabel('Distance')
    ylabel('Correlation')
    title('Variability 0 degree correlation')
end

if saveFlag == 1
    savefig(strcat(savePath, 'One_module_variable_properties.fig'));
end

if oneD_flag ~= 1
    dist = cell2mat(dist_allVariable) * arenaResolution;
    corr = cell2mat(corr_allVariable); 

    dist_binned = 0:gridSpacing(1):(sqrt(2) * (arenaSize(1)/2)); 
    corr_binned = zeros(1, length(dist_binned));

    [~, bins] = min(abs(dist_binned - dist'), [], 2); 
    
    for ii = 1:length(dist_binned)
        corr_binned(ii) = median(corr(bins == ii));
    end
end

figure 
plot(dist_binned, corr_binned, 'ko', 'MarkerFaceColor', 'k')
axis([-0.2, arenaSize(1)/2 + 0.2, 0.5, 1.05])
xlabel('Distance (m)')
ylabel('Correlation')


















