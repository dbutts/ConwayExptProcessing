%% read in plexon filepath
clear;
pl2path = '/home/bizon/Data/'; %'plexon .pl2 path
filenameP = '250110_114356_Jacomo'; % '250218_173539_Jacomo 241210_154655_Jacomo' experiment name
addpath(genpath('/home/bizon/Git/ConwayExptProcessing')) % add folder to path
thisSessionFile = [pl2path filenameP '.pl2'];
%% read from plexon file
[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI05'); % XR
[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI06'); % YR
[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07'); % XL
[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08'); % YL

[~, ~, ~, ~, PlexET_ad(5,:)] = plx_ad_v(thisSessionFile, 'AI01'); % sync strobes

eyetrace = PlexET_ad(1:5,:)*(90/1000); % convert to arcmin *gains/(analog scale)


%% plot eyetraces only during fixation - Cameron
winSize = 1:length(eyetrace); % 1:length(eyetrace) or 1000;  % the size of the window
num_windows = floor(length(eyetrace)/range(winSize)); % how many windows
start_window =  1; % First window
x_analysis = true; % true = run analysis in x direction, false = run analysis in y direction
quick_plot = false; % whether or not to not plot saccades and correction - long time in larger datasets
speedThres = 100000; % 3 % assumed speed of each saccade - set to very high number to ignore saccades 
saccadeLen = 10; % assumed length of each saccade
x_fixation = [-25,25]; % need to set x and y during each experiment
y_fixation = [-25,25]; 

for ii = start_window:num_windows % iterate over windows
    figure(1)
    time2plot = winSize + max(winSize)*(ii-1);

    subplot(8,1,1);
        YL = eyetrace(4, time2plot);%-mean(eyetrace(4, time2plot));
        plot(YL); hold on;
        YR = eyetrace(2, time2plot);%-mean(eyetrace(2, time2plot));
        plot(YR); hold off;
        legend({'Left Eye','Right Eye'})
        xlim([0, max(time2plot)-min(time2plot)])
        ylabel('Y trace (arcmin)')
        ylim([-1000 1000]) % scale to be a constant ylim to compare between windows easier

    subplot(8,1,2); % plot eyetraces
        XL = eyetrace(3, time2plot);%-mean(eyetrace(3, time2plot));
        plot(XL); hold on;
        XR = eyetrace(1, time2plot);%-mean(eyetrace(1, time2plot));
        plot(XR); hold off;
        legend({'Left Eye','Right Eye'})
        ylabel('X trace (arcmin)')
        xlim([0, max(time2plot)-min(time2plot)])
        ylimEye = ylim;
        ylim([-1000 1000]) % comment out for saccade detection
    
    disp('Subplot 3')
    subplot(8,1,3); % plot differences
        if x_analysis
            eyediff = XL - XR;
            plot(eyediff)
            ylabel('X difference (arcmin)')
        else
            eyediff = YL - YR;
            plot(eyediff)
            ylabel('Y difference (arcmin)')
        end
        xSqrt = sqrt(mean(eyediff.^2));
        xlim([0, max(time2plot)-min(time2plot)])
        ylimDiff = ylim;
    
    disp('Subplot 4')
    subplot(8,1,4) % plot speed
        if x_analysis
            XL_vel = abs(diff(eyetrace(3, time2plot))); % take derivative of eyetrace to find speed
            XR_vel = abs(diff(eyetrace(1, time2plot)));
            plot(XL_vel); hold on
            plot(XR_vel); hold off
        else
            YL_vel = abs(diff(eyetrace(4, time2plot))); % take derivative of eyetrace to find speed
            YR_vel = abs(diff(eyetrace(2, time2plot)));
            plot(YL_vel); hold on
            plot(YR_vel); hold off
        end
        legend({'Left Eye','Right Eye'})
        ylabel('Speed')
        xlabel('Time (ms)')
        xlim([0, max(time2plot)-min(time2plot)])
    
    
    subplot(8,1,5) % plot saccades
        if x_analysis
            saccades = (XL_vel > speedThres) | (XR_vel > speedThres);
        else
            saccades = (YL_vel > speedThres) | (YR_vel > speedThres);
        end
        if ~quick_plot
            plot(saccades)
            ylabel('Saccades')
            xlim([0, max(time2plot)-min(time2plot)])
        end
    
    
    transitions = diff([1 saccades 1]); % pad the start and end and take derivative
    start_indices = find(transitions == -1);  % Transition from 1 to 0
    end_indices = find(transitions == 1);  % Transition from 0 to 1
    startends = [start_indices+saccadeLen; end_indices-saccadeLen]' + min(time2plot); % concatenate, add saccade length, and add the start window
    
    subplot(8,1,6) % plot saccade indicies
        if ~quick_plot
            scatter(start_indices, zeros(length(start_indices)), 'g'); hold on
            scatter(end_indices, zeros(length(end_indices)),'r', 'MarkerEdgeAlpha', 0.5); hold off
            ylabel('Start (green) end (red)')
            xlim([0, max(time2plot)-min(time2plot)])
        end
    
    disp('Subplot 7')
    subplot(8,1,7) % plot eyetraces divided by saccades
    L_cor = NaN(1,length(time2plot)); 
    R_cor = NaN(1,length(time2plot));
        for i=1:size(startends,1)
            indx = startends(i,1):startends(i,2);
            if x_analysis
                L_cor(indx-min(time2plot)) = eyetrace(3, indx); % Left_corrected eyetraces after removing the saccades
                R_cor(indx-min(time2plot)) = eyetrace(1, indx);
                % L_cor(indx-min(time2plot)) = eyetrace(3, indx) - mean(eyetrace(3,indx));
                % R_cor(indx-min(time2plot)) = eyetrace(1, indx) - mean(eyetrace(1,indx));
            else
                % L_cor(indx-min(time2plot)) = eyetrace(4, indx) - mean(eyetrace(4,indx));
                % R_cor(indx-min(time2plot)) = eyetrace(2, indx) - mean(eyetrace(2,indx));
                L_cor(indx-min(time2plot)) = eyetrace(4, indx);
                R_cor(indx-min(time2plot)) = eyetrace(2, indx);
            end
        end
        % Set values to NaN that are outside of our window
        if x_analysis
            L_cor(L_cor < x_fixation(1)) = NaN;
            L_cor(L_cor > x_fixation(2)) = NaN;
            R_cor(R_cor < x_fixation(1)) = NaN;
            R_cor(R_cor > x_fixation(2)) = NaN;
        else
            L_cor(L_cor < y_fixation(1)) = NaN;
            L_cor(L_cor > y_fixation(2)) = NaN;
            R_cor(R_cor < y_fixation(1)) = NaN;
            R_cor(R_cor > y_fixation(2)) = NaN;
        end
        if ~quick_plot
            plot(L_cor); hold on
            plot(R_cor); hold off
            legend({'Left Eye','Right Eye'})
            xlabel('Divided by saccades and restricted to (0,0)')
            if x_analysis
                ylabel('X trace (Arcmin)')
            else 
                ylabel('Y trace (Arcmin)')
        end    
        xlim([0, max(time2plot)-min(time2plot)])
        ylim(ylimEye); % plot with same limits as original traces
        end
    disp('Subplot 8')
    subplot(8,1,8) % plot square differences
        diff_cor = L_cor - R_cor;
        if ~quick_plot
            plot(diff_cor)
            if x_analysis
                ylabel('X difference (arcmin)')
            else 
                ylabel('Y difference (arcmin)')
            end
        end
        % only include valid (not NaN) values in root means squared
        % calculation
        validL = ~isnan(L_cor);
        validR = ~isnan(R_cor);
        validIndicies = validL & validR;
        xdiff_valid = L_cor(validIndicies) - R_cor(validIndicies);
        xSqrt_cor = sqrt(mean(xdiff_valid.^2)); % find root mean squared of difference
        xlim([0, max(time2plot)-min(time2plot)])
        ylim(ylimDiff)

    titleText = sprintf(['start = %d end = %d ' ...
                         'sqr diff = %.4f ' ...
                         'corrected sqr diff = %.4f\n ' ...
                         'saccade length = %d ' ...
                         'speed Thresh = %d isX = %d'], ...
                         min(time2plot), max(time2plot), xSqrt, xSqrt_cor, saccadeLen, speedThres, x_analysis);
    sgtitle(titleText)
    if num_windows ~= 1
        pause
    end
end

%% plot 2D histogram - Cameron
l_analysis = true; % true = run analysis in x direction, false = run analysis in y direction

figure(2)
if l_analysis
    % XL = eyetrace(3, :);
    % YL = eyetrace(4, :);
    [counts, xedges, yedges] = histcounts2(XL, YL, 300);  % 300 bins for each axis
else
    % XR = eyetrace(1, :);
    % YR = eyetrace(2, :);
    [counts, xedges, yedges] = histcounts2(XR, YR, 300);  % 300 bins for each axis
end
% Shift xedges and yedges to the center of each bin
xcenters = (xedges(1:end-1) + xedges(2:end)) / 2;
ycenters = (yedges(1:end-1) + yedges(2:end)) / 2;

% Create a 2D histogram plot using pcolor
pcolor(xcenters, ycenters, counts'); 
shading flat;  % Remove grid lines and make it look smooth
xlabel('X-axis');
ylabel('Y-axis');
title(['2D Histogram using histcounts2 left analysis =' num2str(l_analysis)]);
colorbar;  % Display color scale


%% plot Fixation points in fivedot

fixaton = zeros(length(ExptTrials),2);
result = cell(length(ExptTrials), 1);
figure(2)
for i = 1:length(ExptTrials)
    currentfixation = ExptTrials{i,1}.m_pt2iFixationSpot; % current fixation spot in pixels
    fixation(i,1:2) = currentfixation;
    % result{i} = ExptTrials{i,1}.m_strctTrialOutcome.m_strResult;
    scatter(currentfixation(1),currentfixation(2))
    xlim([840 1080]) % found in ExptTrials{1, 1}.apt2iFixationSpots  
    ylim([420 660])
    disp(['current fixation x y: ' num2str(currentfixation)]);
    trialTime = ExptTrials{i, 7}.ET_times(1)-ExptTrials{i, 7}.ET_times(end);
    disp(['Current trial time: ' num2str(trialTime)]);
    disp(ExptTrials{i, 7}.ET_times(1))
    pause
end
scatter(fixation(:,1),fixation(:,2))

%% plot strobes, pupil, eyetraces

figure(2)
for i = 1:length(ExptTrials)
    subplot(4,1,1)
    plot(ExptTrials{i,7}.ET_trace(1,:)); 
    legend({'Sync strobe'})

    subplot(4,1,2)
    plot(ExptTrials{i,7}.ET_trace(3,:)); hold on
    plot(ExptTrials{i,7}.ET_trace(4,:)); hold off
    legend({'x pupil','y pupil'})

    subplot(4,1,3)
    plot(ExptTrials{i,7}.ET_trace(5,:)); hold on
    plot(ExptTrials{i,7}.ET_trace(7,:)); hold off
    legend({'R x','L x'})

    subplot(4,1,4)
    plot(ExptTrials{i,7}.ET_trace(6,:)); hold on
    plot(ExptTrials{i,7}.ET_trace(8,:)); hold off
    legend({'R y','L y'})

    pause
end