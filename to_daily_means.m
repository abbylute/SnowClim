% convert data on hours_in_ts time step to daily means
function [data_daily_mean] = to_daily_means(data_on_timestep, hours_in_ts)
    data_reshape = reshape(data_on_timestep,...
        24/hours_in_ts,...
        size(data_on_timestep,1)/(24/hours_in_ts),...
        size(data_on_timestep,2));
    data_daily_mean = squeeze(nanmean(data_reshape,1));
end
