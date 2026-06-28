function [ decoding_error_vs_distance, field_distances ] = DISTANCE_DECODING(rateMap, n_decode, arenaSize, arenaResolution, gridSpacing, thresh)
%-------------------------------------------------------------------------%
%   This function samples noisy ratemaps as Poisson processes from the
%   idealized ratemap. From these noisy ratemaps, a cross-validated
%   average distance vs. correlation template is constructed. The remaining
%   ratemap is treated then as a "test". Inspired by Tennant et al.'s 
%   (2018) experiment, we look for the first location where the observed
%   correlation is within thresh % of the template's target distance. 
%
%   Written by WTR 10/09/2025 // Last updated by WTR 10/09/2025
%-------------------------------------------------------------------------%
%% Generate nDecode samples 
[n_x, n_y, n_N] = size(rateMap); 

X = zeros(n_x, n_y, n_N, n_decode); 
for ii = 1:n_N
    for jj = 1:n_decode
        X(:, :, ii, jj) = poissrnd(rateMap(:, :, ii)); 
    end
end

%% Identify location of grid peaks
distances = (-arenaSize(1)/2):arenaResolution:(arenaSize(1)/2);
bins_per_field = round(gridSpacing / arenaResolution);
ids = 1:n_x;
id_fields = [ids(ceil(n_x/2 + bins_per_field):bins_per_field:end)];
field_distances = distances(id_fields);
n_fields = length(id_fields);

%% Performing leave one out decoding 
decoding_error_vs_distance = nan(n_fields, n_decode);

for ii = 1:n_decode
    tr_ids = [1:(ii - 1), (ii + 1):n_decode];
    te_ids = ii;

    tr_dist_corr = zeros(n_x, length(tr_ids));
    for jj = 1:length(tr_ids)
        tr_dist_corr(:, jj) = DISTANCE_CORR(X(:, :, :, tr_ids(jj)));
    end
    tr_dist_corr = median(tr_dist_corr, [2]); tr_dist_corr = tr_dist_corr(id_fields);
    
    te_X = X(:, :, :, te_ids);
    te_dist_corr = DISTANCE_CORR(te_X); te_dist_corr = te_dist_corr(id_fields);

    for jj = 1:n_fields
        corr_rel_diff = abs(te_dist_corr(jj) - tr_dist_corr) / te_dist_corr(jj); 
        sim_ids = find(corr_rel_diff < thresh);
        if ~isempty(sim_ids)
            decoded_dist = field_distances(min(sim_ids));
            decoding_error_vs_distance(jj, ii) =  decoded_dist - field_distances(jj);
        end
    end


end

end