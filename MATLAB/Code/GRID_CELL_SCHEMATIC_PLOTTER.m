%% Grid cell schematic plotter
%-------------------------------------------------------------------------%
%   This script plots the activity of three grid cells in the red, green,
%   and blue channels (repsectively), as well as plots their joint
%   activity. Used to illustrate how small differences in grid parameters
%   can lead to complex interference patterns that contain local spatial
%   information. 
%
%   Written by WTR 11/18/2023 // Last updated by WTR 01/17/2025
%-------------------------------------------------------------------------%
%% 
function GRID_CELL_SCHEMATIC_PLOTTER(rateMap, neurons2plot, gridSpacing, arenaResolution, arenaSize, saveFlag, fig_id)

% Plotting rate maps overlaid
empty_channel = zeros(size(rateMap, 1), size(rateMap, 2));
red_channel = rateMap(:, :, neurons2plot(1)) / max(rateMap(:, :, neurons2plot(1)), [], 'all');
green_channel = rateMap(:, :, neurons2plot(2)) / max(rateMap(:, :, neurons2plot(2)), [], 'all');

figure
imagesc(cat(3, red_channel, green_channel, empty_channel))
if saveFlag
    savefig(strcat('Figures/Grid cell distance schematic figure/', fig_id, '_joint_grid_cell_activity.fig'))
end
pause

% Plotting activity of each grid cell at one of the grid's peaks 
[~, nY, ~] = size(rateMap);
centerY = ceil(nY / 2);
distance = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
bins_per_field = gridSpacing / arenaResolution;
ids = 1:size(rateMap, 1);
id_fields = [flip(ids(ceil(size(rateMap, 1)/2):(-bins_per_field):1)), ids(ceil(size(rateMap, 1)/2 + bins_per_field):bins_per_field:end)];

red_activity = red_channel(id_fields, centerY);
green_activity = green_channel(id_fields, centerY);

figure
plot(distance(id_fields), red_activity, 'ro-', 'MarkerFaceColor', 'r'); hold on 
plot(distance(id_fields), green_activity, 'go-', 'MarkerFaceColor', 'g')
xlabel('Distance in physical space')
ylabel('Activity')
axis([-2, 2, 0, 1.2])
if saveFlag
    savefig(strcat('Figures/Grid cell distance schematic figure/', fig_id, '_activity_at_grid_fields.fig'))
end
pause

% Plotting distance in activity space as a function of physical distance
figure
plot(distance(id_fields), sqrt((red_activity - green_activity).^2), 'ko-', 'MarkerFaceColor', 'k'); hold on 
xlabel('Distance in physical space')
ylabel('Distance in neural activity space')
axis([-2, 2, 0, 1.5])
if saveFlag
    savefig(strcat('Figures/Grid cell distance schematic figure/', fig_id, '_neural_distance_at_grid_fields.fig'))
end
pause


end