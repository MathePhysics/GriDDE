function [corr0, corr90, corr45, corr_all, dist_all] = DISTANCE_CORR(rateMap, gridField_locs)
%-------------------------------------------------------------------------%
%   This function computes the population correlation, as a function of
%   distance from the center of the environment, 
%
%   Written by WTR 01/13/2025 // Last updated by WTR 06/02/2026
%-------------------------------------------------------------------------%
%%
[nX, nY, ~] = size(rateMap);
centerX = ceil(nX / 2);
centerY = ceil(nY / 2);

corr0 = zeros(1, nX);
corr90 = zeros(1, nY);
corr45 = zeros(1, min([nX, nY]));
corr_all = zeros(1, length(gridField_locs));
dist_all = zeros(1, length(gridField_locs));


norm_center = sqrt(squeeze(rateMap(centerX, centerY, :))' * squeeze(rateMap(centerX, centerY, :)));

for xx = 1:nX
    norm_xx = sqrt(squeeze(rateMap(xx, centerY, :))' * squeeze(rateMap(xx, centerY, :)));
    corr0(xx) = (squeeze(rateMap(xx, centerY, :))' * squeeze(rateMap(centerX, centerY, :))) / (norm_xx * norm_center);
end

for yy = 1:nY
    norm_yy = sqrt(squeeze(rateMap(centerX, yy, :))' * squeeze(rateMap(centerX, yy, :)));
    corr90(yy) = (squeeze(rateMap(centerX, yy, :))' * squeeze(rateMap(centerX, centerY, :))) / (norm_yy * norm_center);
end

for zz = 1:min([nX, nY])
    norm_zz = sqrt(squeeze(rateMap(zz, zz, :))' * squeeze(rateMap(zz, zz, :)));
    corr45(zz) = (squeeze(rateMap(zz, zz, :))' * squeeze(rateMap(centerX, centerY, :))) / (norm_zz * norm_center);
end

for ii = 1:length(gridField_locs)
    norm_ii = sqrt(squeeze(rateMap(gridField_locs(ii, 1), gridField_locs(ii, 2), :))' * squeeze(rateMap(gridField_locs(ii, 1), gridField_locs(ii, 2), :)));
    corr_all(ii) = (squeeze(rateMap(gridField_locs(ii, 1), gridField_locs(ii, 2), :))' * squeeze(rateMap(centerX, centerY, :))) / (norm_ii * norm_center);
    dist_all(ii) = sqrt((gridField_locs(ii, 1) - centerX)^2 + (gridField_locs(ii, 2) - centerY)^2);
end
