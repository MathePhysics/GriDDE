function [corr0, corr90] = DISTANCE_CORR_DECODING(rateMap1, rateMap2)
%-------------------------------------------------------------------------%
%   This function computes the population correlation, as a function of
%   distance from the center of the environment, using one ratemap as a
%   template and the other as a set of observations.
%
%   Written by WTR 01/13/2025 // Last updated by WTR 01/13/2025
%-------------------------------------------------------------------------%
%%
[nX, nY, ~] = size(rateMap1);
centerX = ceil(nX / 2);
centerY = ceil(nY / 2);

corr0 = zeros(1, nX);
corr90 = zeros(1, nY);

norm_center = sqrt(squeeze(rateMap1(centerX, centerY, :))' * squeeze(rateMap1(centerX, centerY, :)));

for xx = 1:nX
    norm_xx = sqrt(squeeze(rateMap2(xx, centerY, :))' * squeeze(rateMap2(xx, centerY, :)));
    corr0(xx) = (squeeze(rateMap2(xx, centerY, :))' * squeeze(rateMap1(centerX, centerY, :))) / (norm_xx * norm_center);
end

for yy = 1:nY
    norm_yy = sqrt(squeeze(rateMap2(centerX, yy, :))' * squeeze(rateMap2(centerX, yy, :)));
    corr90(yy) = (squeeze(rateMap2(centerX, yy, :))' * squeeze(rateMap1(centerX, centerY, :))) / (norm_yy * norm_center);
end