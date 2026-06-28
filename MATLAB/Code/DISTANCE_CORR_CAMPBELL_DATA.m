function [corr0] = DISTANCE_CORR_CAMPBELL_DATA(rateMap)
%-------------------------------------------------------------------------%
%   This function computes the population correlation, as a function of
%   distance from the start of the environment, on real data recorded by
%   Campbell et al. (2021)
%
%   Written by WTR 05/27/2026 // Last updated by WTR 05/27/2026
%-------------------------------------------------------------------------%
%%
[nX, ~] = size(rateMap);
start = 1;

corr0 = zeros(1, nX);

norm_start = sqrt(rateMap(start, :) * rateMap(start, :)');

for xx = 1:nX
    norm_xx = sqrt(rateMap(xx, :) * rateMap(xx, :)');
    corr0(xx) = (rateMap(xx, :) * rateMap(start, :)') / (norm_start * norm_xx);
end

