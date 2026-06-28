%% Tennant et al. 2018 behavioral data 
%-------------------------------------------------------------------------%
%   This script loads in the publicly available data associated with
%   Tennant et al. 2018 Cell Reports and generates a plot analogous to
%   their Figure 3D. 
%
%   Written by WTR 09/24/2025 // Last updated by WTR 04/30/2026
%-------------------------------------------------------------------------%
%% Loading data
data_path = '/Users/willredman/Documents/Grid cell variability distance coding/Data/';
file_name = 'Figure3_D_0100.csv';

X = readmatrix(strcat(data_path, file_name));

%% Processing data 
n_mice = max(X(:, 1));
reward_zone_lengths = unique(X(:, 2));
n_reward_zones = length(reward_zone_lengths);

D = reshape(X(:, 3:end), [n_reward_zones, n_mice, 3]);

%% Plotting 
results_path = '/Users/willredman/Documents/Grid cell variability distance coding/Results/';
save_flag = 1;

figure; hold on
fill([reward_zone_lengths(1), reward_zone_lengths(end), reward_zone_lengths(end), reward_zone_lengths(1)], [reward_zone_lengths(1) - 20, reward_zone_lengths(end) - 20, reward_zone_lengths(end), reward_zone_lengths(1)], 'g-', 'FaceAlpha', 0.3)
for ii = 1:n_mice
    for jj = 1:n_reward_zones
        plot(reward_zone_lengths(jj), D(jj, ii, 1), 'b^');
        plot(reward_zone_lengths(jj), D(jj, ii, 2), 'co');
        plot(reward_zone_lengths(jj), D(jj, ii, 3), 'rs');
    end
end
plot(reward_zone_lengths, nanmean(D(:, :, 1), 2), 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
plot(reward_zone_lengths, nanmean(D(:, :, 2), 2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8)
plot(reward_zone_lengths, nanmean(D(:, :, 3), 2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
xlabel('Location of reward zone (cm)')
ylabel('Location (cm)')
if save_flag == 1
    savefig(strcat(results_path, 'Tennant_behavioral_data.fig'));
end

figure; hold on 
plot(reward_zone_lengths / 100, nanmean(D(:, :, 2) - reward_zone_lengths, 2) / 100, 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8); hold on 
for ii = 1:n_reward_zones
    plot([reward_zone_lengths(ii), reward_zone_lengths(ii)] / 100, ...
        [nanmean(D(ii, :, 2) - reward_zone_lengths(ii), 2) - nanstd(D(ii, :, 2) - reward_zone_lengths(ii)), ...
        nanmean(D(ii, :, 2) - reward_zone_lengths(ii), 2) + nanstd(D(ii, :, 2) - reward_zone_lengths(ii))] / 100, 'c-');
end
xlabel('Location of reward zone (m)')
ylabel('Error (m)')
title('Error vs distance (Non-beaconed): Mouse')
axis([0, 5.0, -2.5, 0.5])
if save_flag == 1
    savefig(strcat(results_path, 'Error_vs_distance_non_beaconed.fig'));
end

figure; hold on 
plot(reward_zone_lengths / 100, nanmean(D(:, :, 3) - reward_zone_lengths, 2) / 100, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on 
for ii = 1:n_reward_zones
    plot([reward_zone_lengths(ii), reward_zone_lengths(ii)] / 100, ...
        [nanmean(D(ii, :, 3) - reward_zone_lengths(ii), 2) - nanstd(D(ii, :, 3) - reward_zone_lengths(ii)), ...
        nanmean(D(ii, :, 3) - reward_zone_lengths(ii), 2) + nanstd(D(ii, :, 3) - reward_zone_lengths(ii))] / 100, 'r-');
end
xlabel('Location of reward zone (m)')
ylabel('Error (m)')
title('Error vs distance (Probe): Mouse')
axis([0, 5.0, -2.5, 0.5])
if save_flag == 1
    savefig(strcat(results_path, 'Error_vs_distance_probe.fig'));
end

%% Doing statistics
n_bootstrap = 50; 

% Looking at probe trials
probe_error = (D(:, :, 3) - reward_zone_lengths); 
probe_error4 = probe_error(4, :); 
probe_error4(isnan(probe_error4)) = [];
probe_error5 = probe_error(5, :); 
probe_error5(isnan(probe_error5)) = [];
mean_error_bootstrapped = zeros(n_bootstrap, 2);

[~,p_probe] = kstest2(probe_error4, probe_error5, 'Tail', 'Larger')

for ii = 1:n_bootstrap
    mean_error_bootstrapped(ii, 1) = mean(probe_error4(randi(length(probe_error4), 1, length(probe_error4))));
    mean_error_bootstrapped(ii, 2) = mean(probe_error5(randi(length(probe_error5), 1, length(probe_error5))));
end

[~,p_probe_bootstrap] = kstest2(mean_error_bootstrapped(:, 1), mean_error_bootstrapped(:, 2), 'Tail', 'Larger')

figure()
histogram(mean_error_bootstrapped(:, 1), binEdges = -200:20:200); hold on 
histogram(mean_error_bootstrapped(:, 2), binEdges = -200:20:200)
legend([{'4'}, {'5'}])
if save_flag == 1
    savefig(strcat(results_path, 'Bootstrap_probe.fig'));
end

% Looking at non-beaconed trials
probe_abs_error = abs(D(:, :, 2) - reward_zone_lengths); 
probe_error4 = probe_abs_error(4, :); 
probe_error4(isnan(probe_error4)) = [];
probe_error5 = probe_abs_error(5, :); 
probe_error5(isnan(probe_error5)) = [];
mean_error_bootstrapped = zeros(n_bootstrap, 2);

[~,p_non_beacon] = kstest2(probe_error4, probe_error5, 'Tail', 'Smaller')

for ii = 1:n_bootstrap
    mean_error_bootstrapped(ii, 1) = mean(probe_error4(randi(length(probe_error4), 1, length(probe_error4))));
    mean_error_bootstrapped(ii, 2) = mean(probe_error5(randi(length(probe_error5), 1, length(probe_error5))));
end

[~,p_non_beacon_bootstrap] = kstest2(mean_error_bootstrapped(:, 1), mean_error_bootstrapped(:, 2), 'Tail', 'Smaller')

figure()
histogram(mean_error_bootstrapped(:, 1), binEdges = 0:20:200); hold on 
histogram(mean_error_bootstrapped(:, 2), binEdges = 0:20:200)
legend([{'4'}, {'5'}])
if save_flag == 1
    savefig(strcat(results_path, 'Bootstrap_non_beacon.fig'));
end




