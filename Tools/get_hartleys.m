function get_hartleys(data_hartleys, target_SUs, lag, save_vars, data)
% the function takes data_hartleys from the previous step in
% MasterFile_LSRProcessing. target_SUs are the single units we want to plot
% lag tells us at what lag we will plot

% Find Single unit response
Robs = double(data_hartleys.Robs);
% Find number of single units = rows in Robs
% nSU = size(data_hartleys.Robs, 1);

for i = target_SUs
    % for each single unit, make a new figure
    fig = figure;
    % change the size of the figure [xposition yposition xsize ysize]
    fig.Position = [0 0 900 700];
        % find plots for j=1 (lum), j=2 (L-M), j=3 (S?)
    % for lag = 1:6
        for j = 1:3
            % Pull ith single unit
            single_unit_response = Robs(i, :);
            % Shift the response earlier in time to align with the stimuli
            % by a certain lag
            single_unit_response = circshift(single_unit_response, -lag);
            % set the first values to zero corresponding to how many lags
            % there are
            single_unit_response(end-lag+1:end) = 0;
            % Find the hartleys that correspond to luminance only
            % Isolate the data_hartleys.hartley_metas which = j in 4th row 
            hartley_meta = data_hartleys.hartley_metas;
            % Find columns (timepoints) that do not correspond to j
            non_j_cols = hartley_meta(4, :) ~=j;
            % Create a new matrix with only the columns which correspond to j
            hartley_meta_j = hartley_meta;
            hartley_meta_j(:, non_j_cols) = 0;
                
            % find plots for k=1 (orientation), k=2 (spatial frequency)
            for k = 1:2
                % Find orientation data - second row of hartley_metas
                % Find spatial frequency - first row of hartley_metas
                orientation_or_freq = hartley_meta_j(k, :);
                
                % Calculate average response for each orientation 
                % Example using grpstats to compute group means
                averaged_response = grpstats(single_unit_response', orientation_or_freq', {'mean'});
                
                % find each possible orientation
                unique_orientations_or_freq = unique(orientation_or_freq)';

                subplot(4, 3, 2*(j-1)+k)
                % legend('Lag 1','Lag 2','Lag 3','Lag 4','Lag 5','Lag 6');

                hold on

                plot(unique_orientations_or_freq, averaged_response)

                if k==1
                    orientation_or_freq_name = 'frequency';
                elseif k==2
                    orientation_or_freq_name = 'orientation';
                end
                if j==1
                    color = ' Luminance';
                elseif j==2
                    color = ' L-M';
                elseif j==3
                    color = ' S?';
                end
                xlabel(orientation_or_freq_name)
                ylabel('response')
                [t,s] = title(['Averaged ', orientation_or_freq_name, ' response'], [color,' SU # ', num2str(i),' Lag = ', num2str(lag),]);
            end
        end
    % end
    % Pull ith single unit
    single_unit_response = Robs(i, :);
    % Shift the response earlier in time to align with the stimuli
    % by a certain lag
    % lag = 5;
    single_unit_response = circshift(single_unit_response, -lag);
    % set the first values to zero corresponding to how many lags
    % there are
    single_unit_response(end-lag+1:end) = 0;
    % Generate a heat map with orientation on y and response on x
    % and the color is the strength of the activation
    % Pull each orientation and frequency
    orientation = hartley_meta_j(1, :);
    frequency = hartley_meta_j(2, :);
    % Find each unique orientation and frequency so we can iterate
    % over them
    unique_orientations = unique(orientation);
    unique_frequencies = unique(frequency);
    % Iterate each unique orientation, frequency, and time point. Add each activation to the variable activation
    for l = 1:length(unique_orientations)
        unique_orientation = unique_orientations(l);
        for m = 1:length(unique_frequencies)
            unique_frequency = unique_frequencies(m);
            orientation_indicies = find(orientation == unique_orientation);
            frequency_indicies = find(frequency == unique_frequency);
            common_indicies = intersect(orientation_indicies, frequency_indicies);
            activations = single_unit_response(common_indicies);
            % Average each activation
            activation = mean(activations, 'all');
            if unique_orientation == 0 
                % disp(activation)
                % disp(activations)
                % check if activations has NaN values
                % disp(any(isnan(activations)))
                % disp(isempty(activations))
                % disp(common_indicies)
                orientation_indicies_0orientation = orientation_indicies;
                frequency_indicies_0orientation = frequency_indicies;
                common_indicies_0orientation = common_indicies;
                activation_0orientation = activation;
                freq_at_0orientation = frequency(orientation_indicies_0orientation);
            end
            % Add each activation to a matrix - orientations on x
            % axis (l) and frequencies on y axis (m)
            activation_matrix(l, m) = activation;
        end
    end
    subplot(4, 3, 9)
    imagesc(activation_matrix)
    ylabel('Orientation')
    xlabel('Frequency')
    x_tick_labels = num2cell(unique_frequencies);
    y_tick_labels = num2cell(unique_orientations);
    xticks(1:length(x_tick_labels));
    xticklabels(x_tick_labels);
    yticks(1:length(y_tick_labels));
    yticklabels(y_tick_labels);
    % Set 'Orientation vs Frequency' as title and rest of info as subtitle
    [t,s] = title('Orientation vs Frequency ', ['SU # ', num2str(i), ' Lag = ', num2str(lag)]);
    if save_vars.to_save==1
        % save the STA as pdf saveas(figure, filename)
    	saveas(fig,[save_vars.outputdir save_vars.titlestr '_lag_' num2str(lag) '_SU_' num2str(i) '_Hartley.pdf'])
        % save the STA as a jpg
        saveas(fig,[save_vars.outputdir save_vars.titlestr '_lag_' num2str(lag) '_SU_' num2str(i) '_Hartley.jpg'])
    end
end

if ~isempty(data)
    % Now plot the relevant STA for comparison - taken from get_sta.m
    binned_SU = [single(data.Robs'), single(data.RobsMU')];
    use_inds = data.valid_data;
    use_inds(end-10:end) = []; %cut last few indices to avoid artifacts
    stim_shift=permute(data.stim,[4 1 2 3]);
    NT=size(data.ETtrace,2);
    stim_deltas = zeros(2,NT); 
    stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
    stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;
    cur_STA1(lag,:) = binned_SU(use_inds+lag,target_SUs)' * stim2(use_inds,:);
    curlim=max(abs(cur_STA1(:)')); 
        if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots
    cur_STA2=reshape(cur_STA1(lag,:),60,180);
    cur_STA2(:,[60 120])=curlim;
    subplot(4,3,7)
    
    %Now, we deviate from get_sta.m to plot only the lum portion
    cur_STA2 = cur_STA2(:,1:60);
    imagesc(cur_STA2); clim([-curlim curlim]); pbaspect([3 1 1])
    title('Lum STA - fix ratio___')
    
    % Now do the 2D fast Fourier transform
    fft_STA = fft2(cur_STA2);
    fft_STA = abs(fftshift(fft_STA));
    subplot(4,3,8)
    imagesc(fft_STA)
    colorbar
    title('2D FFT of STA')
    
    
    % divide the graph into quadrants
    % search for the pair of quadrants that have the largest peak
    % Find the size of the matrix
    [m,n] = size(fft_STA);
    % Divide the matrix into four quadrants
    TopLeft = fft_STA(1:m/2, 1:n/2);
    TopRight = fft_STA(1:m/2, n/2+1:end);
    BottomLeft = fft_STA(m/2+1:end, 1:n/2);
    BottomRight = fft_STA(m/2+1:end, n/2+1:end);
    % Determine the maximum in each Quadrant
    % Find row and column of Top Left max (TL)
    [max_TL, linearIndex] = max(TopLeft(:));
    [row_TL, col_TL] = ind2sub(size(TopLeft), linearIndex);
    % Find row and column of Top Right max (TR)
    [max_TR, linearIndex] = max(TopRight(:));
    [row_TR, col_TR] = ind2sub(size(TopRight), linearIndex);
    % Adjust top right column so that the additional columns are added in
    col_TR = col_TR + n/2;
    % Find row and column of Bottom right
    [max_BR, linearIndex] = max(BottomRight(:));
    [row_BR, col_BR] = ind2sub(size(BottomRight), linearIndex);
    % Adjust bottom right matrix rows and columns
    row_BR = row_BR + m/2;
    col_BR = col_BR + n/2;
    % Find row and column of Bottom Left
    [max_BL, linearIndex] = max(BottomLeft(:));
    [row_BL, col_BL] = ind2sub(size(BottomRight), linearIndex);
    % Adjust row of bottom left
    row_BL = row_BL + m/2;
    % pull the corresponding maximum value from fourier
    % Main diagonal = top left to bottom right
    mainDiagonalAvg = (max_TL + max_BR)/2;
    counterDiagonalAvg = (max_TR + max_BL)/2;
    subplot(4,3,10)
    imagesc(fft_STA)
    title('Slope from quadrants')
    hold on
    if mainDiagonalAvg > counterDiagonalAvg
        % slope of main diagnoal
        slope = -(col_BR-col_TL)/(row_BR-row_TL);
        % disp(['Main diagonal slope:' num2str(slope)])
        plot([col_BR, col_TL], [row_BR, row_TL], 'r', 'LineWidth', 2);
        legend(['Slope: ', num2str(slope)])
    else
        % slope of counter diagonal
        slope = -(col_TR-col_BL)/(row_TR-row_BL);
        % disp(['Counter diagonal slope:' num2str(slope)])
        % plot a line between the fourier peaks and extend it using the slope
        plot([col_TR, col_BL], [row_TR, row_BL], 'r', 'LineWidth', 2);
        legend(['Slope: ', num2str(slope)])
    end
    
    subplot(4,3,11)
    imagesc(cur_STA2)
    title('Lum STA with Slope Line')
    hold on;
    x1=30; y1=30; % middle of plot
    x2=x1+30; y2=y1+30*slope;
    line([x1,x2], [y1,y2], 'Color', 'r');
    x3=0; y3=y1-30*slope;
    line([x1,x3],[y1,y3], 'Color','r');

end