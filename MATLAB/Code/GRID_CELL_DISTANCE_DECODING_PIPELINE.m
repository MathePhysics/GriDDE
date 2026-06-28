%% GRID CELL DISTANCE DECODING PIPELINE
%-------------------------------------------------------------------------%
%   This script generates a population of grid cells, with varying
%   properties, and uses the population correlation vector (with respect to
%   the starting location) to decode distance. 
%
%   Based on the simulation code from Redman et al. eLife (2024). Decoding 
%   analysis is inspired by the experiments of Tennant et al. (2018). 
%
%   Written by WTR 10/09/2025 // Last updated by WTR 05/12/2026
%-------------------------------------------------------------------------%
close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals      
gridSpacing = 0.8; % 0.8                             % mean spacing of grid cells - in meters    
gridOrientation = 0; %0                          % mean orientation of grid cells - in degrees
arenaSize = [20, 0]; %[20, 0]                            % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = [64]; % [64]                               % size of populations of grid cells to be tested
gridFiring_max = 15;  % 15                           % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = 0; %0                                     % standard deviation used for sampling grid orientation
spacingStd = 0.05; % 0.05                              % standard deviation used for sampling grid spacing
arenaResolution = 0.01; % 0.01                       % size of each bin used for discretizing space
nSamples = 7;  % 7                                  % number of grid populations to sample from
nDecode = 10; % 10                                   % number of noisy ratemaps to sample from
poissonFlag = 1;                                     % flag to sample from a Poisson distribution 
decodingThresh = 0.05; % 0.05                       % threshold for determining when there is sufficient similarity for decoding

saveFlag = 0;                                    % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of size of grid population on decoding accuracy
for ss = 1:nSamples
    ss
    [rateMap, ~, ~, ~] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, spacingStd, arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, poissonFlag, oneD_flag);
    [error, field_distances] = DISTANCE_DECODING(rateMap, nDecode, arenaSize, arenaResolution, gridSpacing, decodingThresh);

    if ss == 1
        nFields = length(field_distances);
        decoding_error_vs_distance = nan(nSamples, nFields, nDecode);
    end

    decoding_error_vs_distance(ss, :, :) = error;

end

% decoding_error_vs_distance = abs(decoding_error_vs_distance);

%% Saving data
if saveFlag == 1
    save(strcat(savePath, '/decoding_accuracy_vs_distance_lambda_', num2str(gridSpacing), '_sigma_', num2str(spacingStd), '_decoding_thresh_', num2str(decodingThresh), '.mat'), "decoding_error_vs_distance");
end

%% Plotting
field_distances_plot = [0.88, 1.18, 1.63, 2.3, 3.31, 4.81];
plot_ids = zeros(1, length(field_distances_plot));
for ii = 1:length(field_distances_plot)
    [~, plot_ids(ii)] = min(abs(field_distances - field_distances_plot(ii)));
end

figure
plot(field_distances(plot_ids), nanmean(nanmean(decoding_error_vs_distance(:, plot_ids, :), [3]), [1]), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on 
for ii = 1:length(plot_ids)
    plot([field_distances(plot_ids(ii)), field_distances(plot_ids(ii))], [nanmean(nanmean(decoding_error_vs_distance(:, plot_ids(ii), :), [3]), [1]) - nanstd(nanmean(decoding_error_vs_distance(:, plot_ids(ii), :), [3]), [1]), ...
        nanmean(nanmean(decoding_error_vs_distance(:, plot_ids(ii), :), [3]), [1]) + nanstd(nanmean(decoding_error_vs_distance(:, plot_ids(ii), :), [3]), [1])], 'k-')
end
ylabel('Error (m)')
xlabel('Location')
axis([0, 5.1, -4.0, 0.5])
title('Error vs distance: Model')
if saveFlag == 1
    savefig(strcat(savePath, '/decoding_accuracy_vs_distance_lambda_', num2str(gridSpacing), '_sigma_', num2str(spacingStd), '_decoding_thresh_', num2str(decodingThresh), '.fig'));
end



