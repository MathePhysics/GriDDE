%% GRID CELL DISTANCE SCALING PIPELINE
%-------------------------------------------------------------------------%
%   This script generates a population of grid cells, with varying
%   properties, and computes the correlation as a function of distance. How
%   different properties of the grid code interact with the distance coding
%   is investigated.
%
%   Based on the simulation code from Redman et al. eLife (2025). 
%
%   Written by WTR 02/20/2025 // Last updated by WTR 10/28/2025
%-------------------------------------------------------------------------%
close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals      
gridSpacing = [0.4:0.05:0.8];         % mean spacing of grid cells - in meters    
gridOrientation = [0];                % mean orientation of grid cells - in degrees
arenaSize = [20, 0];                 % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = 64;                % size of populations of grid cells to be tested
gridFiring_max = 15;                % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = 0;                         % standard deviation used for sampling grid orientation
spacingStd = [0.01:0.005:0.05];                  % standard deviation used for sampling grid spacing
arenaResolution = 0.05;               % size of each bin used for discretizing space
nSamples = 10;                      % number of grid populations to sample from
poissonFlag = 1;                    % flag to sample from a Poisson distribution 
baseline_percentile = 1;

saveFlag = 0;                       % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of size of grid population on distance coding
FWHM = zeros(length(gridSpacing), length(spacingStd));
if oneD_flag == 1
    corr0 = zeros(length(gridSpacing), length(spacingStd), nSamples, arenaSize(1) / arenaResolution + 1);
end

for ll = 1:length(gridSpacing)
    ll
    for dd = 1:length(spacingStd)
        for ss = 1:nSamples
            [rateMap_Variable, spacing_Variable, orientation_Variable, phase_Variable] = ...
                RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing(ll), spacingStd(dd), arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, poissonFlag, oneD_flag);
            [corr0Variable, ~] = DISTANCE_CORR(rateMap_Variable);
            corr0(ll, dd, ss, :) = corr0Variable;
        end

        distance0 = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
        bins_per_field0 = round(gridSpacing(ll) / arenaResolution);
        ids0 = 1:size(rateMap_Variable, 1);
        id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];

        C = squeeze(median(corr0(ll, dd, :, id_fields0), 3));
        
        sorted_C = sort(C);
        baseline_corr = prctile(C, baseline_percentile);
        peak_height = sorted_C(end - 1) -  baseline_corr;
        half_max = peak_height / 2 + baseline_corr;
        [~, half_max_id] = min(abs(half_max - C));
        FWHM(ll, dd) =  abs(distance0(id_fields0(half_max_id)));

    end
end

%% Saving data
if saveFlag 
    if oneD_flag == 1
        save(strcat(savePath, 'HWHM_vs_grid_properties_1d.mat'), 'FWHM');
        save(strcat(savePath, 'corr0_vs_grid_properties_1d.mat'), 'corr0');
    else
        save(strcat(savePath, 'HWHM_vs_grid_properties_2d.mat'), 'FWHM');
    end
end

%% Plotting
% FWHM vs grid properties
figure
imagesc(FWHM); hold on
title('FWHM vs grid properties')
yt = get(gca, 'YTick');
yt_label = gridSpacing; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('Grid spacing (m)')
xt = get(gca, 'XTick');
xt_label = spacingStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\Delta_\lambda (m)')

% Plotting example correlation curves for different grid spacings 
C = mean(corr0, 3);
distance0 = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
grid_spacings_plot = [1:length(gridSpacing)];

figure 
legend_key = [];
for ii = 1:length(grid_spacings_plot)
    bins_per_field0 = round(gridSpacing(ii) / arenaResolution);
    ids0 = 1:size(rateMap_Variable, 1);
    id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];
    plot(distance0(id_fields0), squeeze(C(ii, 3, id_fields0)), 'o-', 'Color', [1 - ii * 0.1, 0, ii * 0.1], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\lambda = ', num2str(gridSpacing(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.7, 1.1])


% Plotting example difference in correlation curves for different grid
% spacings 
figure 
legend_key = [];
for ii = 1:length(grid_spacings_plot)
    bins_per_field0 = round(gridSpacing(ii) / arenaResolution);
    ids0 = 1:size(rateMap_Variable, 1);
    id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];    
    backwards_DoC = [abs(squeeze(C(ii, 3, id_fields0(2:(end))) - C(ii, 3, id_fields0(1:(end - 1)))))', nan];
    forwards_DoC = [nan, abs(squeeze(C(ii, 3, id_fields0(2:(end))) - C(ii, 3, id_fields0(1:(end - 1)))))'];
    DoC = nanmean([backwards_DoC; forwards_DoC], 1);
    plot(distance0(id_fields0), DoC, '-', 'Color', [1 - ii * 0.1, 0, ii * 0.1], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\lambda = ', num2str(gridSpacing(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.0, 0.1])

% Plotting example correlation curves for different grid spacing
% variabilities
grid_spacings_std_plot = [1:length(spacingStd)];

figure 
legend_key = [];
for ii = 1:length(grid_spacings_std_plot)
    bins_per_field0 = round(gridSpacing(3) / arenaResolution);
    ids0 = 1:size(rateMap_Variable, 1);
    id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];
    plot(distance0(id_fields0), squeeze(C(3, ii, id_fields0)), '-', 'Color', [1 - ii * 0.1, 0, ii * 0.1], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\sigma_\lambda = ', num2str(spacingStd(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.7, 1.1])

% Plotting example difference in correlation curves for different grid
% spacing variabilities
figure 
legend_key = [];
for ii = 1:length(grid_spacings_std_plot)
    bins_per_field0 = round(gridSpacing(3) / arenaResolution);
    ids0 = 1:size(rateMap_Variable, 1);
    id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];
    backwards_DoC = [abs(squeeze(C(3, ii, id_fields0(2:(end))) - C(3, ii, id_fields0(1:(end - 1)))))', nan];
    forwards_DoC = [nan, abs(squeeze(C(3, ii, id_fields0(2:(end))) - C(3, ii, id_fields0(1:(end - 1)))))'];
    DoC = nanmean([backwards_DoC; forwards_DoC], 1);
    plot(distance0(id_fields0), DoC, '-', 'Color', [1 - ii * 0.1, 0, ii * 0.1], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\sigma_\lambda = ', num2str(spacingStd(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.0, 0.1])

% Computing maximum difference in correlation (and its location)
max_DoC = zeros(length(gridSpacing), length(spacingStd));
max_DoC_loc = zeros(length(gridSpacing), length(spacingStd));

for ii = 1:length(gridSpacing)
    for jj = 1:length(spacingStd)
        bins_per_field0 = round(gridSpacing(ii) / arenaResolution);
        ids0 = 1:size(rateMap_Variable, 1);
        id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];    
        backwards_DoC = [abs(squeeze(C(ii, jj, id_fields0(2:(end))) - C(ii, jj, id_fields0(1:(end - 1)))))', nan];
        forwards_DoC = [nan, abs(squeeze(C(ii, jj, id_fields0(2:(end))) - C(ii, jj, id_fields0(1:(end - 1)))))'];
        DoC = nanmean([backwards_DoC; forwards_DoC], 1);
        DoC = DoC((length(DoC) / 2):end);
        id_fields0 = id_fields0((length(id_fields0) / 2):end);
        [max_DoC(ii, jj), max_DoC_loc_id] = max(DoC(3:end));
        max_DoC_loc(ii, jj) = distance0(id_fields0(max_DoC_loc_id + 2));
    end
end

% Plotting maximum difference in correlation 
figure
imagesc(max_DoC)
title('Maximum DoC vs grid properties')
yt = get(gca, 'YTick');
yt_label = gridSpacing; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('Grid spacing (m)')
xt = get(gca, 'XTick');
xt_label = spacingStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\Delta_\lambda (m)')
clim([0, 0.05])


% Plotting maximum difference in correlation location
figure
imagesc(max_DoC_loc)
title('Maximum DoC location vs grid properties')
yt = get(gca, 'YTick');
yt_label = gridSpacing; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('Grid spacing (m)')
xt = get(gca, 'XTick');
xt_label = spacingStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\Delta_\lambda (m)')
clim([0, 12])



