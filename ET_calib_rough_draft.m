% Shift according to Kofiko eye traces

calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), {trial.m_strTrialType});

upSampleFactor = 2;
Kofiko_Xpix_frameRate_cellArray = ...
    cellfun(@(t,x, start, stop, fixspot, n) interp1(t, x - fixspot, linspace(start, stop, n), 'linear'),...
    KofikoET_cellArrays.Kofiko_ET_TS_PlexonTime_cellArray(2*find(ccloudTrialIdx)),...
    KofikoET_cellArrays.Kofiko_Xpix_cellArray(2*find(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStartTimes(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStopTimes(ccloudTrialIdx)),...
    X_fixationSpot(ccloudTrialIdx),...
    num2cell(upSampleFactor * stimTiming.numFrames(ccloudTrialIdx)), ...
    'UniformOutput', false);

Kofiko_Ypix_frameRate_cellArray = ...
    cellfun(@(t,y, start, stop, fixspot, n) interp1(t, y - fixspot, linspace(start, stop, n), 'linear'),...
    KofikoET_cellArrays.Kofiko_ET_TS_PlexonTime_cellArray(2*find(ccloudTrialIdx)),...
    KofikoET_cellArrays.Kofiko_Ypix_cellArray(2*find(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStartTimes(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStopTimes(ccloudTrialIdx)),...
    Y_fixationSpot(ccloudTrialIdx),...
    num2cell(upSampleFactor * stimTiming.numFrames(ccloudTrialIdx)), ...
    'UniformOutput', false);

% ddpi

dpi_fname = '/Volumes/lsr-conway/DATA/monkey_ephys/Sprout/260622/RawDDPI-2026Jun22-155430/RawDDPI-2026Jun22-155430.txt';
ET = readtable(dpi_fname);
% ddpi timestamps
t_dpi = ET.RightSeconds;
t_dpi = t_dpi - t_dpi(1); % time stamps in s
% ddpi sync signal
sync_dpi = ET.Int0;
sync_dpi = sync_dpi - min(sync_dpi);
sync_dpi = sync_dpi & 1;
sync_dpi_diff = diff(sync_dpi);
t_rising_dpi = t_dpi(sync_dpi_diff > 0);
t_falling_dpi = t_dpi(sync_dpi_diff < 0);
% get rid of unpaired falling edge at beginning or rising edge at end
if size(t_rising_dpi,1) > size(t_falling_dpi,1) % unpaired rising edge at the end
    t_rising_dpi = t_rising_dpi(1:size(t_falling_dpi,1));
elseif size(t_rising_dpi,1) < size(t_falling_dpi,1)
    t_falling_dpi = t_falling_dpi(2:size(t_rising_dpi,1));
end
dt_dpi = median(t_falling_dpi - t_rising_dpi);
delays_dpi = (t_rising_dpi(2:end) - t_falling_dpi(1:end-1))./dt_dpi;
delays_dpi_int = int8(delays_dpi);
sync_ad = PlexET_ad_calib(:,1);
thresh = mean([min(sync_ad), max(sync_ad)]);
sync_ad_bin = sync_ad > thresh;
sync_ad_diff = diff(sync_ad_bin);
t_rising_plexon = t_plexon(sync_ad_diff > 0);
t_falling_plexon = t_plexon(sync_ad_diff < 0);
if size(t_rising_plexon,1) > size(t_falling_plexon,1) % unpaired rising edge at the end
    t_rising_plexon = t_rising_plexon(1:size(t_falling_plexon,1));
elseif size(t_rising_plexon,1) < size(t_falling_plexon,1)
    t_falling_plexon = t_falling_plexon(2:size(t_rising_plexon,1));
end
dt_plexon = median(t_falling_plexon - t_rising_plexon);
% delays
delays_plexon = (t_rising_plexon(2:end) - t_falling_plexon(1:end-1))./dt_plexon;
delays_plexon_int = int8(delays_plexon);
% alignment
% cross-correlate delays
[c,lags] = xcorr(delays_plexon_int, delays_dpi_int);
best_lag = lags(c == max(c));
x = t_rising_plexon(best_lag+1:end);
t_rising_plexon_matched = x(1:numel(t_rising_dpi));
b_dpi_plexon = [ones(size(t_rising_dpi)) t_rising_dpi]\t_rising_plexon_matched';
% put dpi signal in plexon time
t_dpi_plexon = [ones(size(t_dpi)) t_dpi]*b_dpi_plexon;
dpi_cellArrays.dpi_ts_PlexonTime_cellArray = binByStimIntervals(t_dpi_plexon, t_dpi_plexon, stimTiming.stimIntervals);
cr_x_Right = ET.RightCR1X;
cr_y_Right = ET.RightCR1Y;
cr_x_Left = ET.LeftCR1X;
cr_y_Left = ET.LeftCR1Y;
p4_x_Right = ET.RightCR4X;
p4_y_Right = ET.RightCR4Y;
p4_x_Left = ET.LeftCR4X;
p4_y_Left = ET.LeftCR4Y;
dpi_raw_Right = [cr_x_Right - p4_x_Right, cr_y_Right - p4_y_Right];
dpi_raw_Left = [cr_x_Left - p4_x_Left, cr_y_Left - p4_y_Left];

dpi_pupilArea_Right = ET.RightPupilWidth .* ET.RightPupilHeight;
dpi_pupilArea_Left =  ET.LeftPupilWidth .* ET.LeftPupilHeight;

dpi_cellArrays.RightX = binByStimIntervals(t_dpi_plexon, dpi_raw_Right(:,1), stimTiming.stimIntervals);
dpi_cellArrays.RightY = binByStimIntervals(t_dpi_plexon, dpi_raw_Right(:,2), stimTiming.stimIntervals);
dpi_cellArrays.LeftX = binByStimIntervals(t_dpi_plexon, dpi_raw_Left(:,1), stimTiming.stimIntervals);
dpi_cellArrays.LeftY = binByStimIntervals(t_dpi_plexon, dpi_raw_Left(:,2), stimTiming.stimIntervals);

dpi_cellArrays.RightPupilArea =  binByStimIntervals(t_dpi_plexon,dpi_pupilArea_Right, stimTiming.stimIntervals);
dpi_cellArrays.LeftPupilArea =  binByStimIntervals(t_dpi_plexon, dpi_pupilArea_Left, stimTiming.stimIntervals);

RPA = vertcat(dpi_cellArrays.RightPupilArea{2*find(isCalibrationTrial)});
LPA = vertcat(dpi_cellArrays.LeftPupilArea{2*find(isCalibrationTrial)});

blinkIdx = RPA < 0.5.*max(RPA) & LPA < 0.5*max(LPA);

fixation_spot = vertcat(trial(isCalibrationTrial).m_pt2iFixationSpot);

%% X

inlierThresh = 95;
saccThresh = 75;

X_fix = fixation_spot(:,1) - 960;

x_samplesPerTrial = cellfun(@numel,dpi_cellArrays.RightX(2*find(isCalibrationTrial)));
X_fix_rep = repelem(X_fix, x_samplesPerTrial);
X_fix_rep_transition = [0; diff(X_fix_rep)] ~= 0;
h = ones(301,1);
X_transitionIdx = circshift(conv(X_fix_rep_transition, h, 'same'), 150);

xr = vertcat(dpi_cellArrays.RightX{2*find(isCalibrationTrial)});
xl = vertcat(dpi_cellArrays.LeftX{2*find(isCalibrationTrial)});

% xr = sgolayfilt(xr, 2, 21);
% xl = sgolayfilt(xl, 2, 21);

xr_diff = [0; diff(xr)];
xl_diff = [0; diff(xl)];

xr_saccIdx = abs(xr_diff) > prctile(abs(xr_diff(~blinkIdx)), saccThresh);
xl_saccIdx = abs(xl_diff) > prctile(abs(xl_diff(~blinkIdx)), saccThresh);

xr_outlierIdx =  abs(xr - median(xr)) > prctile(abs(xr - median(xr)), inlierThresh);
xl_outlierIdx =  abs(xl - median(xl)) > prctile(abs(xl - median(xl)), inlierThresh);

xr_idx = ~xr_outlierIdx & ~blinkIdx & ~xr_saccIdx & ~X_transitionIdx;
xl_idx = ~xl_outlierIdx & ~blinkIdx & ~ xl_saccIdx & ~X_transitionIdx;


%% Y

Y_fix = fixation_spot(:,2) - 540;

y_samplesPerTrial = cellfun(@numel,dpi_cellArrays.RightY(2*find(isCalibrationTrial)));
Y_fix_rep = repelem(Y_fix, y_samplesPerTrial);
Y_fix_rep_transition = [0; diff(Y_fix_rep)] ~= 0;
h = ones(301,1);
Y_transitionIdx = circshift(conv(Y_fix_rep_transition, h, 'same'), 150);

yr = vertcat(dpi_cellArrays.RightY{2*find(isCalibrationTrial)});
yl = vertcat(dpi_cellArrays.LeftY{2*find(isCalibrationTrial)});

% yr = sgolayfilt(yr, 2, 21);
% yl = sgolayfilt(yl, 2, 21);

yr_diff = [0; diff(yr)];
yl_diff = [0; diff(yl)];

yr_saccIdx = abs(yr_diff) > prctile(abs(yr_diff(~blinkIdx)), saccThresh);
yl_saccIdx = abs(yl_diff) > prctile(abs(yl_diff(~blinkIdx)), saccThresh);

yr_outlierIdx =  abs(yr - median(yr)) > prctile(abs(yr - median(yr)), inlierThresh);
yl_outlierIdx =  abs(yl - median(yl)) > prctile(abs(yl - median(yl)), inlierThresh);

yr_idx = ~yr_outlierIdx & ~blinkIdx & ~yr_saccIdx & ~Y_transitionIdx;
yl_idx = ~yl_outlierIdx & ~blinkIdx & ~ yl_saccIdx & ~Y_transitionIdx;

%%% FIT
xr_included = xr(xr_idx & yr_idx);
yr_included = yr(xr_idx & yr_idx);

B_xr = [ones(size(xr_included)), xr_included, xr_included.*yr_included, xr_included.^2] \ X_fix_rep(xr_idx & yr_idx);
B_yr = [ones(size(yr_included)), yr_included, xr_included.*yr_included, yr_included.^2] \ Y_fix_rep(xr_idx & yr_idx);

xr_calib = B_xr(1) + B_xr(2)*xr + B_xr(3).*xr.*yr+ B_xr(4)*xr.^2;
yr_calib = B_yr(1) + B_yr(2)*yr + B_yr(3).*xr.*yr + B_yr(4)*yr.^2;

yl_included = yl(xl_idx & yl_idx);
xl_included = xl(xl_idx & yl_idx);

B_xl = [ones(size(xl_included)), xl_included, xl_included.*yl_included, xl_included.^2] \ X_fix_rep(xl_idx & yl_idx);
B_yl = [ones(size(yl_included)), yl_included, xl_included.*yl_included, yl_included.^2] \ Y_fix_rep(xl_idx &yl_idx);

xl_calib = B_xl(1) + B_xl(2)*xl + B_xl(3).*xl.*yr + B_xl(4)*xl.^2;
yl_calib = B_yl(1) + B_yl(2)*yl + B_yl(3).*xl.*yl + B_yl(4)*yl.^2;

%% PLOT
figure, scatter(xr_calib, yr_calib, 'filled', 'o', 'MarkerFaceColor', 'b','MarkerFaceAlpha',0.05);
hold on
scatter(X_fix, Y_fix, 'ro', 'lineWidth', 2);
axis equal
xlim([-200 200])
ylim([-200 200])

figure, scatter(xl_calib, yl_calib, 'filled', 'o', 'MarkerFaceColor', 'b','MarkerFaceAlpha',0.05);
hold on
scatter(X_fix, Y_fix, 'ro', 'lineWidth', 2);
axis equal
xlim([-200 200])
ylim([-200 200])

%% correct cloud trials

xr_cc = vertcat(dpi_cellArrays.RightX{2*find(ccloudTrialIdx)});
yr_cc = vertcat(dpi_cellArrays.RightY{2*find(ccloudTrialIdx)});
xl_cc = vertcat(dpi_cellArrays.LeftX{2*find(ccloudTrialIdx)});
yl_cc = vertcat(dpi_cellArrays.LeftY{2*find(ccloudTrialIdx)});
% 
xr_cc = sgolayfilt(xr_cc, 2, 21);
xl_cc = sgolayfilt(xl_cc, 2, 21);
yr_cc = sgolayfilt(yr_cc, 2, 21);
yl_cc = sgolayfilt(yl_cc, 2, 21);

LPA = vertcat(dpi_cellArrays.LeftPupilArea{2*find(ccloudTrialIdx)});
RPA = vertcat(dpi_cellArrays.RightPupilArea{2*find(ccloudTrialIdx)});

blinkIdx_Right =  RPA < 0.5.*max(RPA);
blinkIdx_Left =  LPA < 0.5.*max(LPA);

xr_cc_calib = B_xr(1) + B_xr(2)*xr_cc + B_xr(3).*xr_cc.*yr_cc+ B_xr(4)*xr_cc.^2;
yr_cc_calib = B_yr(1) + B_yr(2)*yr_cc + B_yr(3).*xr_cc.*yr_cc+ B_yr(4)*yr_cc.^2;

xl_cc_calib = B_xl(1) + B_xl(2)*xl_cc + B_xl(3).*xl_cc.*yl_cc+ B_xl(4)*xl_cc.^2;
yl_cc_calib = B_yl(1) + B_yl(2)*yl_cc + B_yl(3).*xl_cc.*yl_cc+ B_yl(4)*yl_cc.^2;

t_dpi_plexon_cc = vertcat(dpi_cellArrays.dpi_ts_PlexonTime_cellArray{2*find(ccloudTrialIdx)});

% filter

xr_cc_calib_diff = [0; diff(xr_cc_calib)];
xl_cc_calib_diff = [0; diff(xl_cc_calib)];

xr_cc_calib_saccIdx = abs(xr_cc_calib_diff) > prctile(abs(xr_cc_calib_diff), saccThresh);
xl_cc_calib_saccIdx = abs(xl_cc_calib_diff) > prctile(abs(xl_cc_calib_diff), saccThresh);

xr_cc_calib_outlierIdx =  abs(xr_cc_calib) > 60; %abs(xr_cc_calib - median(xr_cc_calib(~blinkIdx_Right))) <prctile(abs(xr_cc_calib(~blinkIdx_Right) - median(xr_cc_calib(~blinkIdx_Right))), 95);
xl_cc_calib_outlierIdx =  abs(xl_cc_calib) > 60;%abs(xl_cc_calib - median(xl_cc_calib(~blinkIdx_Left))) <prctile(abs(xl_cc_calib(~blinkIdx_Left) - median(xl_cc_calib(~blinkIdx_Left))), 95);

yr_cc_calib_diff = [0; diff(yr_cc_calib)];
yl_cc_calib_diff = [0; diff(yl_cc_calib)];

yr_cc_calib_saccIdx = abs(yr_cc_calib_diff) > prctile(abs(yr_cc_calib_diff), saccThresh);
yl_cc_calib_saccIdx = abs(yl_cc_calib_diff) > prctile(abs(yl_cc_calib_diff), saccThresh);

yr_cc_calib_outlierIdx =  abs(yr_cc_calib) > 60;%abs(yr - median(yr)) < prctile(abs(yr - median(yr)), inlierThresh);
yl_cc_calib_outlierIdx =  abs(yl_cc_calib) > 60; %abs(yl - median(yl)) < prctile(abs(yl - median(yl)), inlierThresh);

xr_cc_calib(blinkIdx_Right | xr_cc_calib_outlierIdx | xr_cc_calib_saccIdx) = nan;
yr_cc_calib(blinkIdx_Right | yr_cc_calib_outlierIdx | yr_cc_calib_saccIdx) = nan;
xl_cc_calib(blinkIdx_Left | xl_cc_calib_outlierIdx | xl_cc_calib_saccIdx) = nan;
yl_cc_calib(blinkIdx_Left | yl_cc_calib_outlierIdx | yl_cc_calib_saccIdx) = nan;

X = mean([xr_cc_calib, xl_cc_calib], 2, 'omitmissing');
Y = mean([yr_cc_calib, yl_cc_calib], 2, 'omitmissing');

% we need to patch the nan values, ill do it by interpolation

sampleIdx = 1:numel(X);

goodX = ~isnan(X);
goodY = ~isnan(Y);

X = interp1(sampleIdx(goodX), X(goodX), sampleIdx, 'linear', 'extrap');
Y = interp1(sampleIdx(goodY), Y(goodY), sampleIdx, 'linear', 'extrap');

X_cellArray = binByStimIntervals(t_dpi_plexon_cc, X, stimTiming.stimIntervals);
Y_cellArray = binByStimIntervals(t_dpi_plexon_cc, Y, stimTiming.stimIntervals);

% now downsample eye data

upSampleFactor = 2;

X_frameRate_cellArray = ...
    cellfun(@(t,x, start, stop, n) interp1(t, x, linspace(start, stop, n), 'linear'),...
    dpi_cellArrays.dpi_ts_PlexonTime_cellArray(2*find(ccloudTrialIdx)),...
    X_cellArray(2*find(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStartTimes(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStopTimes(ccloudTrialIdx)),...
    num2cell(upSampleFactor * stimTiming.numFrames(ccloudTrialIdx)), ...
    'UniformOutput', false);

Y_frameRate_cellArray = ...
    cellfun(@(t,y, start, stop, n) interp1(t, y, linspace(start, stop, n), 'linear'),...
    dpi_cellArrays.dpi_ts_PlexonTime_cellArray(2*find(ccloudTrialIdx)),...
    Y_cellArray(2*find(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStartTimes(ccloudTrialIdx)),...
    num2cell(stimTiming.stimStopTimes(ccloudTrialIdx)),...
    num2cell(upSampleFactor * stimTiming.numFrames(ccloudTrialIdx)), ...
    'UniformOutput', false);