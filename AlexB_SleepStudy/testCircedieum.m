        my_favourite_colour     = [0.8500 0.3250 0.0980];
        subjectID = subID;
        hemisphere = hemi;
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');

        time_stamps = tmData.actTime(:);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.2999e+09) = nan;

        values = mSunful;

        figure
        set(gcf,'Units','Normalized','Position',[.2 .4 .6 .25])
        plot_zscored_timeseries(time_stamps, values, my_favourite_colour)

        figure
        set(gcf,'Units','Normalized','Position',[.3 .3 .4 .3])

        time_res        = 1; % time resolution (in hours) for periodogram
        max_period      = 72; % Maximum period of 1 week = 168 hours
        do_normalise    = true; % Whether to normalise the periodogram

        % Calculate periodogram
        [psd_estimate, time_periods] = circadian_periodogram(time_stamps, values, time_res, max_period);

        % Plot the periodogram
        plot_periodogram(psd_estimate, time_periods, do_normalise, my_favourite_colour)

        time_res        = 1;
        n_shuffles      = 200;
        shuffle_type    = 'circshift';

        % Get the proportion of variance explained by time of day
        var_explained = variance_explained_by_timeofday(time_stamps, values, time_res);

        % Get the variance explained for shuffled data n_shuffles times to see whether var_explained is significant
        [~, var_explained_p] = get_shuffled_var_explained(time_stamps, values, time_res, n_shuffles, shuffle_type);

%         % Plot a fit based on time of day to the data across days
%         figure
%         plot_timeofday_fit(time_stamps, values, time_res, my_favourite_colour)
% 
%         % Add the variance explained by time of day & p-val to the figure title
%         title(['Var explained by TOD: ' num2str(var_explained) ', p =' num2str(var_explained_p)])
% 
%         figure
%         plot_timeofday_fit(time_stamps, values, time_res, my_favourite_colour,'polar')
% 
%         time_res    = 1;
%         stat        = 'mean';
% 
%         figure
%         circadian_rose(time_stamps, values, time_res, stat, my_favourite_colour)
% 
%         time_res            = 1; % Temporal resolution of time bins (= heatmap pixels)
%         percentile_cutoff   = 2; % For the colour scale of the heatmap, ignore the top and bottom x% of data
% 
%         [circadian_matrix, ~] = make_circadian_matrix(time_stamps, values, time_res);
% 
%         figure
%         plot_circadian_matrix(circadian_matrix, percentile_cutoff, my_favourite_colour);
% 
%         close all
% 
%         circadian_summary_figure(time_stamps, values,time_res)

        %         % Get shuffled data points
        %         shuffled_data_points    = within_day_shuffle(time_stamps, values, 'circshift');
        %
        %         [circadian_matrix, time_edges] = make_circadian_matrix(time_stamps, shuffled_data_points, time_res);
        %
        %
        %         figure
        %         plot_circadian_matrix(circadian_matrix, percentile_cutoff, my_favourite_colour);