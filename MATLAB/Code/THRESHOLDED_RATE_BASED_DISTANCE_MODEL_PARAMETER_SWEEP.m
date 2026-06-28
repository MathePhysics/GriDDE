%% THRESHOLDED RATE BASED DISTANCE MODEL PIPELINE PARAMETER SWEEP
%-------------------------------------------------------------------------%
%   This script uses a basic thresholded spiking rate model to encode
%   distance in a semi-biologically plausible manner and performs a
%   parameter sweep to determine how grid spacing, grid spacing 
%   variability, and threshold determine distance specificity. 
%
%   Written by WTR 05/13/2026 // Last updated by WTR 05/13/2026
%-------------------------------------------------------------------------%
close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals 
gridOrientation = 0; %0                          % mean orientation of grid cells - in degrees
arenaSize = [20, 0]; %[20, 0]                            % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = [128]; % [64]                               % size of populations of grid cells to be tested
gridFiring_max = 15;  % 15                           % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = 0; %0                                     % standard deviation used for sampling grid orientation
arenaResolution = 0.05; % 0.01                       % size of each bin used for discretizing space
nSamples = 25;  % 7                                  % number of grid populations to sample from
nDecode = 10; % 10                                   % number of noisy ratemaps to sample from
poissonFlag = 1;                                     % flag to sample from a Poisson distribution 
saveFlag = 1;                                    % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of grid spacing and spiking threshold on distance specificity 
gridSpacing = [0.4:0.05:0.8];   
iThresh = [0.5:0.05:0.90];   
spacingStd = [0.02];

dist = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
median_peakWidth = zeros(length(gridSpacing), length(iThresh));
maximum_peakLocation = zeros(length(gridSpacing), length(iThresh));

for ii = 1:length(gridSpacing)
    ii
    for jj = 1:length(iThresh)
        spikingRate = zeros(nSamples, ceil(arenaSize(1)/arenaResolution) + 1);
        for ss = 1:nSamples
            [rateMap, ~, ~, ~] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing(ii), spacingStd, arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, poissonFlag, oneD_flag);
            
            synapticInput = sum(rateMap, [3]);
            thresholded_synapticInput = synapticInput; 
            thresholded_synapticInput(thresholded_synapticInput < (nNeurons * gridFiring_max * iThresh(jj))) = 0;
            spikingRate(ss, :) = thresholded_synapticInput;
        end

        norm_spikingRate = spikingRate ./ max(spikingRate, [], 2);       
        [peaks, peak_locs, widths, proms] = findpeaks(mean(norm_spikingRate, [1]), dist, MinPeakProminence = 0.3);
        median_peakWidth(ii, jj) = median(widths(peak_locs > 0));
        maximum_peakLocation(ii, jj) = max(peak_locs);

    end 
end

%% Saving data
if saveFlag == 1
    save(strcat(savePath, '/median_peak_width_vs_grid_spacing.mat'), "median_peakWidth");
    save(strcat(savePath, '/maximum_peakLocation_vs_grid_spacing.mat'), "maximum_peakLocation");
end

%% Plotting 
figure
imagesc(median_peakWidth)
ylabel('Grid spacing')
yt = get(gca, 'YTick');
yt_label = gridSpacing; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
xlabel('I_{thresh}')
xt = get(gca, 'XTick');
xt_label = iThresh; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
title('Grid spacing vs. FWHM')
colorbar()
if saveFlag == 1
    savefig(strcat(savePath, '/grid_spacing_vs_FWHM_spiking_model_sigma_lambda_', num2str(spacingStd(1)), '.fig'));
end

figure
imagesc(maximum_peakLocation)
ylabel('Grid spacing')
yt = get(gca, 'YTick');
yt_label = gridSpacing; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
xlabel('I_{thresh}')
xt = get(gca, 'XTick');
xt_label = iThresh; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
title('Grid spacing vs. max. peak location')
colorbar()
if saveFlag == 1
    savefig(strcat(savePath, '/grid_spacing_vs_max_peak_location_spiking_model_sigma_lambda_', num2str(spacingStd(1)), '.fig'));
end

%% Grid cell width vs. rate based distance model peaks
gridSpacing = [0.4:0.05:0.8];  
spacingStd = 0.02;
dist = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
median_gridWidth = zeros(1, length(gridSpacing)); 

for ii = 1:length(gridSpacing)
    ii
    [rateMap, ~, ~, ~] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing(ii), spacingStd, arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, 0, oneD_flag);
    
    gridWidths = zeros(1, nNeurons);
    for jj = 1:nNeurons
        [peaks, peak_locs, widths, proms] = findpeaks(rateMap(:, :, jj), dist, MinPeakProminence = 0.3);
        gridWidths(jj) = median(widths(peak_locs > 0));
    end
    median_gridWidth(ii) = median(gridWidths);
end

%% Plotting 
figure 
mdl1 = fitlm(median_gridWidth', median_peakWidth(:, iThresh == 0.5));
mdl2 = fitlm(median_gridWidth', median_peakWidth(:, iThresh == 0.7));
mdl3 = fitlm(median_gridWidth', median_peakWidth(:, iThresh == 0.9));
plot(mdl1); hold on 
plot(mdl2);
plot(mdl3)
plot([0.15, 0.45], [0.15, 0.45], 'k-')
axis([0.15, 0.45, 0, 0.45])
xlabel('\lambda')
ylabel('FWHM')
if saveFlag == 1
    savefig(strcat(savePath, '/grid_width_vs_spiking_model_width_spacing_std_', num2str(spacingStd), '.fig'));
end

%% Role of grid spacing variability and spiking threshold on distance specificity 
gridSpacing = [0.5];   
iThresh = [0.5:0.05:0.90];   
spacingStd = [0.01:0.005:0.05];

dist = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
median_peakWidth = zeros(length(gridSpacing), length(iThresh));
maximum_peakLocation = zeros(length(gridSpacing), length(iThresh));

for ii = 1:length(spacingStd)
    ii
    for jj = 1:length(iThresh)
        spikingRate = zeros(nSamples, ceil(arenaSize(1)/arenaResolution) + 1);
        for ss = 1:nSamples
            [rateMap, ~, ~, ~] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, spacingStd(ii), arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, poissonFlag, oneD_flag);
            synapticInput = sum(rateMap, [3]);
            thresholded_synapticInput = synapticInput; 
            thresholded_synapticInput(thresholded_synapticInput < (nNeurons * gridFiring_max * iThresh(jj))) = 0;
            spikingRate(ss, :) = thresholded_synapticInput;
        end

        norm_spikingRate = spikingRate ./ max(spikingRate, [], 2);       
        [peaks, peak_locs, widths, proms] = findpeaks(mean(norm_spikingRate, [1]), dist, MinPeakProminence = 0.3);
        median_peakWidth(ii, jj) = median(widths(peak_locs > 0));
        maximum_peakLocation(ii, jj) = max(peak_locs);

    end 
end

%% Saving data
if saveFlag == 1
    save(strcat(savePath, '/median_peak_width_vs_grid_spacing_std.mat'), "median_peakWidth");
    save(strcat(savePath, '/maximum_peakLocation_vs_grid_spacing_std.mat'), "maximum_peakLocation");
end

%% Plotting 
figure
imagesc(median_peakWidth)
ylabel('Grid spacing variability')
yt = get(gca, 'YTick');
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
xlabel('I_{thresh}')
xt = get(gca, 'XTick');
xt_label = iThresh; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
title('Grid spacing variability vs. FWHM')
colorbar()
if saveFlag == 1
    savefig(strcat(savePath, '/grid_spacing_std_vs_FWHM_spiking_model_lambda_', num2str(gridSpacing(1)), '.fig'));
end

figure
imagesc(maximum_peakLocation)
ylabel('Grid spacing variability')
yt = get(gca, 'YTick');
yt_label = spacingStd; 
set(gca, 'YTick', yt, 'YTickLabel', yt_label);
xlabel('I_{thresh}')
xt = get(gca, 'XTick');
xt_label = iThresh; 
set(gca, 'XTick', xt, 'XTickLabel', xt_label);
title('Grid spacing variability vs. max. peak location')
colorbar()
if saveFlag == 1
    savefig(strcat(savePath, '/grid_spacing_std_vs_max_peak_location_spiking_model_lambda_', num2str(gridSpacing), '.fig'));
end


