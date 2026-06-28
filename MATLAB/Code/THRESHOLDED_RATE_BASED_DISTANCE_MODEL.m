%% THRESHOLDED RATE BASED DISTANCE MODEL PIPELINE
%-------------------------------------------------------------------------%
%   This script uses a basic thresholded spiking rate model to encode
%   distance in a semi-biologically plausible manner. 
%
%   Written by WTR 05/12/2026 // Last updated by WTR 05/13/2026
%-------------------------------------------------------------------------%
close all 
clear all
 
%% Set random seed 
rng(1)

%% Globals      
gridSpacing = 0.5; % 0.8                             % mean spacing of grid cells - in meters    
gridOrientation = 0; %0                          % mean orientation of grid cells - in degrees
arenaSize = [20, 0]; %[20, 0]                            % size of arena - in meters       
oneD_flag = arenaSize(2) == 0;                   % whether or not the arena is 1D
nNeurons = [128]; % [64]                               % size of populations of grid cells to be tested
gridFiring_max = 15;  % 15                           % maximum firing rate for idealized model (Poisson neurons can fire more)
oriStd = 0; %0                                     % standard deviation used for sampling grid orientation
spacingStd = 0.02; % 0.05                              % standard deviation used for sampling grid spacing
arenaResolution = 0.05; % 0.01                       % size of each bin used for discretizing space
nSamples = 25;  % 7                                  % number of grid populations to sample from
nDecode = 10; % 10                                   % number of noisy ratemaps to sample from
poissonFlag = 1;                                     % flag to sample from a Poisson distribution 
iMin = 0.5;                                 % minimum amount of synaptic input needed before neuron spikes
saveFlag = 0;                                    % whether or not to save all plots and matrices
savePath = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';

%% Role of spiking threshold on distance specificity 
spikingRate = zeros(nSamples, ceil(arenaSize(1)/arenaResolution) + 1);

for ss = 1:nSamples
    ss
    [rateMap, ~, ~, ~] = RATE_MAP_MULTI_FOURIER(gridFiring_max, gridSpacing, spacingStd, arenaSize, gridOrientation, oriStd, nNeurons, arenaResolution, poissonFlag, oneD_flag);
    synapticInput = sum(rateMap, [3]);
    thresholded_synapticInput = synapticInput; 
    thresholded_synapticInput(thresholded_synapticInput < (nNeurons * gridFiring_max * iMin)) = 0;
    spikingRate(ss, :) = thresholded_synapticInput;
end

norm_spikingRate = spikingRate ./ max(spikingRate, [], 2);
dist = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
[peaks, peak_locs, widths, proms] = findpeaks(mean(norm_spikingRate, [1]), dist, MinPeakProminence = 0.3);
median_peak_width = median(widths(peak_locs > 0))

%% Plotting 
figure 
plot(dist, mean(norm_spikingRate, [1]), 'k-'); hold on
plot(peak_locs, peaks, 'ro')
xlim([0, arenaSize(1)/2])
xlabel('Distance (m)')
ylabel('Mean normalized spiking rate')
title('Spike rate vs. distance')
if saveFlag == 1
    savefig(strcat(savePath, '/spiking_rate_vs_distance_lambda_', num2str(gridSpacing), '_sigma_', num2str(spacingStd), '_iMin_', num2str(iMin), '.fig'));
end


