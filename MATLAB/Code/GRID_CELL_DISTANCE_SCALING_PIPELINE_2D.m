%% GRID CELL DISTANCE SCALING PIPELINE 2D
%-------------------------------------------------------------------------%
%   This script generates a population of grid cells, with varying
%   properties, and computes the correlation as a function of distance. How
%   different properties of the grid code interact with the distance coding
%   is investigated.
%
%   Based on the simulation code from Redman et al. eLife (2025). 
%
%   Written by WTR 10/28/2025 // Last updated by WTR 06/03/2025
%-------------------------------------------------------------------------%
close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals      
gridSpacing = [0.5];         % mean spacing of grid cells - in meters    
gridOrientation = [0];                % mean orientation of grid cells - in degrees
arenaSize = [30, 30];                 % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = 64;                % size of populations of grid cells to be tested
gridFiring_max = 15;                % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = [0.5:0.5:2];                         % standard deviation used for sampling grid orientation
spacingStd = [0.01:0.01:0.05];                  % standard deviation used for sampling grid spacing
arenaResolution = 0.05;               % size of each bin used for discretizing space
nSamples = 25;                      % number of grid populations to sample from
poissonFlag = 1;                    % flag to sample from a Poisson distribution 
baseline_percentile = 1;

saveFlag = 1;                       % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of size of grid population on distance coding
FWHM = zeros(length(spacingStd), length(oriStd));
corr0 = zeros(length(spacingStd), length(oriStd), nSamples, arenaSize(1) / arenaResolution + 1);
corr_all = cell(length(spacingStd), length(oriStd), nSamples);
dist_all = cell(length(spacingStd), length(oriStd), nSamples);

for dd = 1:length(spacingStd)
    dd
    for oo = 1:length(oriStd)
        for ss = 1:nSamples
            [rateMap_Variable, spacing_Variable, orientation_Variable, phase_Variable, gridField_locsVariable] = ...
                RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, spacingStd(dd), arenaSize, gridOrientation, oriStd(oo), nNeurons, arenaResolution, poissonFlag, oneD_flag);
            [corr0Variable, corr90Variable, corr45Variable, corr_allVariable, dist_allVariable] = DISTANCE_CORR(rateMap_Variable, gridField_locsVariable);
            corr0(dd, oo, ss, :) = corr0Variable;
            corr_all{dd, oo, ss} = corr_allVariable;
            dist_all{dd, oo, ss} = dist_allVariable;
        end

        distance0 = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
        bins_per_field0 = round(gridSpacing / arenaResolution);
        ids0 = 1:size(rateMap_Variable, 1);
        id_fields0 = [flip(ids0(ceil(size(rateMap_Variable, 1)/2):(-bins_per_field0):1)), ids0(ceil(size(rateMap_Variable, 1)/2 + bins_per_field0):bins_per_field0:end)];

        C = squeeze(median(corr0(dd, oo, :, id_fields0), 3));
        
        sorted_C = sort(C);
        baseline_corr = prctile(C, baseline_percentile);
        peak_height = sorted_C(end - 1) -  baseline_corr;
        half_max = peak_height / 2 + baseline_corr;
        [~, half_max_id] = min(abs(half_max - C));
        FWHM(dd, oo) =  abs(distance0(id_fields0(half_max_id)));

        

    end
end

%% Saving data
if saveFlag == 1
    if oneD_flag == 1
        save(strcat(savePath, 'HWHM_vs_grid_properties_1d.mat'), 'FWHM');
        save(strcat(savePath, 'corr0_vs_grid_properties_1d.mat'), 'corr0');
    else
        save(strcat(savePath, 'HWHM_vs_grid_properties_2d.mat'), 'FWHM');
        save(strcat(savePath, 'corr0_vs_grid_properties_2d.mat'), 'corr0');
        save(strcat(savePath, 'corr_all_vs_grid_properties_2d.mat'), 'corr_all');
        save(strcat(savePath, 'dist_all_vs_grid_properties_2d.mat'), 'dist_all');
    end
end

%% Plotting
% FWHM vs grid properties
figure
imagesc(FWHM); hold on
title('FWHM vs grid properties')
yt = get(gca, 'YTick');
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick');
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')

% Plotting example correlation curves for different grid spacing
% variabilities
C = mean(corr0, 3);

figure 
legend_key = [];
for ii = 1:length(spacingStd)
    bins_per_field0 = round(gridSpacing / arenaResolution);
    ids0 = 1:size(corr0, 4);
    id_fields0 = [flip(ids0(ceil(size(corr0, 4)/2):(-bins_per_field0):1)), ids0(ceil(size(corr0, 4)/2 + bins_per_field0):bins_per_field0:end)];
    plot(distance0(id_fields0), squeeze(C(ii, 1, id_fields0)), '-', 'Color', [1 - ii * 0.1, 0, ii * 0.1], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\sigma_\lambda = ', num2str(spacingStd(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.7, 1.1])

% Plotting example correlation curves for different grid orientation
% variabilities
C = mean(corr0, 3);

figure 
legend_key = [];
for ii = 1:length(oriStd)
    bins_per_field0 = round(gridSpacing / arenaResolution);
    ids0 = 1:size(corr0, 4);
    id_fields0 = [flip(ids0(ceil(size(corr0, 4)/2):(-bins_per_field0):1)), ids0(ceil(size(corr0, 4)/2 + bins_per_field0):bins_per_field0:end)];
    plot(distance0(id_fields0), squeeze(C(1, ii, id_fields0)), '-', 'Color', [1 - ii * 0.2, 0, ii * 0.2], 'LineWidth', 1.5); hold on 
    legend_key = [legend_key, {strcat('\sigma_\theta = ', num2str(oriStd(ii)))}];
end
legend(legend_key)
axis([0, 10, 0.5, 1.1])

% Computing maximum difference in correlation (and its location)
max_DoC = zeros(length(spacingStd), length(oriStd));
max_DoC_loc = zeros(length(spacingStd), length(oriStd));

for ii = 1:length(spacingStd)
    for jj = 1:length(oriStd)
        bins_per_field0 = round(gridSpacing / arenaResolution);
        ids0 = 1:size(corr0, 4);
        id_fields0 = [flip(ids0(ceil(size(corr0, 4)/2):(-bins_per_field0):1)), ids0(ceil(size(corr0, 4)/2 + bins_per_field0):bins_per_field0:end)];    
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
yt = get(gca, 'YTick'); yt = yt(2:2:end);
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick'); xt = xt(2:2:end);
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')
clim([0, 0.10])

% Plotting maximum difference in correlation location
figure
imagesc(max_DoC_loc)
title('Maximum DoC location vs grid properties')
yt = get(gca, 'YTick'); yt = yt(2:2:end);
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick'); xt = xt(2:2:end);
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')
clim([0, 5])

%% All distance correlation properties
% Computing all distance correlation FWHM 
FWHM_all = zeros(length(spacingStd), length(oriStd));
max_DoC_all = zeros(length(spacingStd), length(oriStd));
max_DoC_loc_all = zeros(length(spacingStd), length(oriStd));

if oneD_flag ~= 1
    for dd = 1:length(spacingStd)
        for oo = 1:length(oriStd)
            dist = [];
            corr = [];
            for nn = 1:nSamples
                dist = [dist, dist_all{dd, oo, nn}];
                corr = [corr, corr_all{dd, oo, nn}];
            end
            dist = dist * arenaResolution;

            dist_binned = 0:gridSpacing(1):(sqrt(2) * (arenaSize(1)/2)); 
            corr_binned = zeros(1, length(dist_binned));

            [~, bins] = min(abs(dist_binned - dist'), [], 2); 
    
            for ii = 1:length(dist_binned)
                corr_binned(ii) = median(corr(bins == ii));
            end

            sorted_C = sort(corr_binned);
            baseline_corr = prctile(corr_binned, baseline_percentile);
            peak_height = sorted_C(end - 1) -  baseline_corr;
            half_max = peak_height / 2 + baseline_corr;
            [~, half_max_id] = min(abs(half_max - corr_binned));
            FWHM_all(dd, oo) =  abs(dist_binned(half_max_id));

            backwards_DoC = [abs(corr_binned(2:end) - corr_binned(1:(end - 1))), nan];
            forwards_DoC = [nan, abs(corr_binned(2:end) - corr_binned(1:(end - 1)))];
            DoC = nanmean([backwards_DoC; forwards_DoC], 1);
            [max_DoC_all(dd, oo), max_DoC_loc_id] = max(DoC(3:end));
            max_DoC_loc_all(dd, oo) = dist_binned(max_DoC_loc_id + 2);

        end
    end
end

figure
imagesc(FWHM_all); hold on
title('FWHM vs grid properties')
yt = get(gca, 'YTick');
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick');
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')

figure
imagesc(max_DoC_all)
title('Maximum DoC vs grid properties')
yt = get(gca, 'YTick'); yt = yt(2:2:end);
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick'); xt = xt(2:2:end);
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')
clim([0, 0.10])

figure
imagesc(max_DoC_loc_all)
title('Maximum DoC location vs grid properties')
yt = get(gca, 'YTick'); yt = yt(2:2:end);
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
ylabel('\sigma_\lambda (m)')
xt = get(gca, 'XTick'); xt = xt(2:2:end);
xt_label = oriStd; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
xlabel('\sigma_\theta (\circ)')
clim([0, 5])





