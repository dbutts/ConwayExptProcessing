addpath(genpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/'));
addpath(genpath('/home/bizon/Processing'));

dirpath = '~/Data/V1_Fovea/Sprout/260313';

filenameP = '260313_163301_Sprout';
plexon_dir = dir(fullfile(dirpath, '**', [filenameP '.pl2']));
plexon_fname = fullfile(plexon_dir(1).folder, plexon_dir(1).name);
pl2 = PL2ReadFileIndex(plexon_fname);

%% Hardcoded values
fps = 60; % display frames per second
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45; 
minSpikes = 2000;

%% Load kofiko data

tic;
fprintf('Loading Kofiko trial data\n')

% Find the mat files produced by Kofiko

kofiko_dir = dir(fullfile(dirpath, '**', [filenameP '*.mat']));
kofiko_fname_re = [filenameP '(_\d*)?.mat'];
valid_kofiko_fname_idx = cellfun(@(x) ~isempty(x), ...
    cellfun(@(x) regexp(x, kofiko_fname_re), {kofiko_dir(:).name}, 'UniformOutput',false));
kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
kofiko_fname = {kofiko_dir(valid_kofiko_fname_idx).name};

% Get the Kofiko mat file numbers from their names
kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'),...
    kofiko_fname, 'UniformOutput',false);
kofiko_file_numStr{cellfun(@isempty, kofiko_file_numStr)} = {'nan'};
kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);
kofiko_file_num(isnan(kofiko_file_num)) = max(kofiko_file_num)+1;

% Load the g_strctLocalExperimentRecording structs which contain
% information about each trials
g_strctLocalExperimentRecording = [];
for i = 1:length(kofiko_file_num)
    idx = kofiko_file_num(i);
    temp_full_fname = fullfile(kofiko_folder{i}, kofiko_fname{i});
    temp_strct = load(temp_full_fname, 'g_strctLocalExperimentRecording');
    temp_cell = temp_strct.g_strctLocalExperimentRecording;
    g_strctLocalExperimentRecording{idx} = temp_cell;
end

% Load the file which should be the "main" Kofiko mat file (no number at
% the end)

mainKofikoFolderIdx = strcmpi(kofiko_fname, [filenameP '.mat']);
load(fullfile(kofiko_folder{mainKofikoFolderIdx}, kofiko_fname{mainKofikoFolderIdx}), 'g_strctDAQParams', ...
    'g_astrctAllParadigms', 'g_strctLog', 'g_strctEyeCalib', ...
    'g_strctAppConfig', 'g_strctSharedParadigmData', ...
    'g_strctStimulusServer', 'g_strctSystemCodes');

% Concatenate the g_strctLocalExperimentRecorindg structs and remove empty
% cells
trialData = vertcat(g_strctLocalExperimentRecording{:});
trialData(cellfun(@isempty, trialData)) = [];

% Unify field names across g_strctLocalExperimentRecording to facilitate
% using cellfun

fieldNames = cellfun(@(x) fieldnames(x), trialData, 'UniformOutput', false);
uniqueFieldNames = unique(vertcat(fieldNames{:}));
missingFieldNames = cellfun(@(x) setdiff(uniqueFieldNames, fieldnames(x)),...
    trialData, 'UniformOutput',false);

for i = 1:numel(trialData)
    for  j = 1:numel(missingFieldNames{i})
        trialData{i}.(missingFieldNames{i}{j}) = [];
    end
end

toc;

%% Time alignment

tic;
fprintf('Loading Kofiko eye signals and converting to plexon time\n')

% Read kofiko strobes from plexon and kofiko files
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);
kofikoSyncStrobesTS = transpose(g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == g_strctSystemCodes.m_iSync));
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iSync);

% Linear regression to get kofiko time stamps into plexon time
B = [ones(size(kofikoSyncStrobesTS)) kofikoSyncStrobesTS]\plexonSyncStrobesTS;

%% Plexon eye data

[adfreq, n, ts, fn, sync_ad] = plx_ad_v(plexon_fname, 'AI01');
[~, ~, ~, ~, leftEyePupil_plexon] = plx_ad_v(plexon_fname, 'AI03');
[~, ~, ~, ~, rightEyePupil_plexon] = plx_ad_v(plexon_fname, 'AI04');
[~, ~, ~, ~, rightEyeX_plexon] = plx_ad_v(plexon_fname, 'AI05');
[~, ~, ~, ~, rightEyeY_plexon] = plx_ad_v(plexon_fname, 'AI06');
[~, ~, ~, ~, leftEyeX_plexon] = plx_ad_v(plexon_fname, 'AI07');
[~, ~, ~, ~, leftEyeY_plexon] = plx_ad_v(plexon_fname, 'AI08');

t_plexon = 0:(1/adfreq):(length(sync_ad)-1)/adfreq;

%% Kofiko eye data

% Get the screen dimensions
ScreenSizeX_pix = g_strctStimulusServer.m_aiScreenSize(3);
ScreenSizeY_pix = g_strctStimulusServer.m_aiScreenSize(4);

% Convert Kofiko eye signal timestamps to Plexon time
Kofiko_ET_TS = g_strctEyeCalib.EyeRaw.TimeStamp';
Kofiko_ET_TS_PlexonTime = [ones(size(Kofiko_ET_TS)) Kofiko_ET_TS]*B;

% Get Kofiko raw eye signals
Kofiko_EyeRawX = g_strctEyeCalib.EyeRaw.Buffer(:,1);
Kofiko_EyeRawY = g_strctEyeCalib.EyeRaw.Buffer(:,2);
% Get Kofiko eye signal gains
Kofiko_GainX = g_strctEyeCalib.GainX.Buffer;
Kofiko_GainY = g_strctEyeCalib.GainY.Buffer;
% Get Kofiko eye signal gain timestamps & convert to plexon time
Kofiko_GainX_TS = g_strctEyeCalib.GainX.TimeStamp;
Kofiko_GainX_TS_PlexonTime = [ones(size(Kofiko_GainX_TS')) Kofiko_GainX_TS']*B;

Kofiko_GainY_TS = g_strctEyeCalib.GainY.TimeStamp;
Kofiko_GainY_TS_PlexonTime = [ones(size(Kofiko_GainY_TS')) Kofiko_GainY_TS']*B;

% Get Kofiko eye signal offsets
Kofiko_CenterX = g_strctEyeCalib.CenterX.Buffer;
Kofiko_CenterY = g_strctEyeCalib.CenterY.Buffer;
% Get Kofiko eye signal offset timestamps
Kofiko_CenterX_TS = g_strctEyeCalib.CenterX.TimeStamp;
Kofiko_CenterX_TS_PlexonTime = [ones(size(Kofiko_CenterX_TS')) Kofiko_CenterX_TS']*B;

Kofiko_CenterY_TS = g_strctEyeCalib.CenterY.TimeStamp;
Kofiko_CenterY_TS_PlexonTime = [ones(size(Kofiko_CenterY_TS')) Kofiko_CenterY_TS']*B;

% Determine the Kofiko gains and offsets for each eye signal time stamp

Kofiko_GainX_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_GainX_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_GainX)
    kidx = Kofiko_ET_TS >= Kofiko_GainX_TS(i);
    pidx = t_plexon >= Kofiko_GainX_TS_PlexonTime(i);

    Kofiko_GainX_forEachTimeStamp(kidx) = Kofiko_GainX(i);
    Kofiko_GainX_forEachPlexonSample(pidx) = Kofiko_GainX(i);
end

Kofiko_GainY_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_GainY_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_GainY)
    kidx = Kofiko_ET_TS >= Kofiko_GainY_TS(i);
    pidx = t_plexon >= Kofiko_GainY_TS_PlexonTime(i);

    Kofiko_GainY_forEachTimeStamp(kidx) = Kofiko_GainY(i);
    Kofiko_GainY_forEachPlexonSample(pidx) = Kofiko_GainY(i);
end

Kofiko_CenterX_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_CenterX_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_CenterX)
    kidx = Kofiko_ET_TS >= Kofiko_CenterX_TS(i);
    pidx = t_plexon >= Kofiko_CenterX_TS_PlexonTime(i);

    Kofiko_CenterX_forEachTimeStamp(kidx) = Kofiko_CenterX(i);
    Kofiko_CenterX_forEachPlexonSample(pidx) = Kofiko_CenterX(i);

end

Kofiko_CenterY_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_CenterY_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_CenterY)
    kidx = Kofiko_ET_TS >= Kofiko_CenterY_TS(i);
    pidx = t_plexon >= Kofiko_CenterY_TS_PlexonTime(i);

    Kofiko_CenterY_forEachTimeStamp(kidx) = Kofiko_CenterY(i);
    Kofiko_CenterY_forEachPlexonSample(pidx) = Kofiko_CenterY(i);
end

% Convert raw eye signals to screen coordinates (pixels) by:
% [RawSignal-SignalOffset]*Gain + ScreenCenter = FixationPoint
% SignalOffset is CenterX or CenterY

Kofiko_Xpix = Kofiko_GainX_forEachTimeStamp.*(Kofiko_EyeRawX - Kofiko_CenterX_forEachTimeStamp) + ScreenSizeX_pix/2;
Kofiko_Ypix = Kofiko_GainY_forEachTimeStamp.*(Kofiko_EyeRawY - Kofiko_CenterY_forEachTimeStamp) + ScreenSizeY_pix/2;

toc;

% apply calibtaion to plexon signals
rightEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*rightEyeX_plexon';
rightEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*rightEyeY_plexon';
leftEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*leftEyeX_plexon';
leftEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*leftEyeY_plexon';

%% trial analysis
tic;
fprintf('Extracting important trial variables\n')

% Find unique trials by finding unique Flip ON timestamps
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
tempImageFlipON_TS_Kofiko_unique = unique(tempImageFlipON_TS_Kofiko);
[~, uniqueTrialIdx, ~] = intersect(tempImageFlipON_TS_Kofiko, tempImageFlipON_TS_Kofiko_unique);
trialData= trialData(uniqueTrialIdx);

% Order trials by flip ON
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
[~, trial_idx_for_sorting] = sort(tempImageFlipON_TS_Kofiko , 'ascend');
trialData = trialData(trial_idx_for_sorting);

% Determine stimulus intervals
StimulusON_MS = cellfun(@(x) x.m_fStimulusON_MS, trialData);
ImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, close trialData);
stimStartTimes = [ones(size(ImageFlipON_TS_Kofiko)), ImageFlipON_TS_Kofiko]*B; % plexon time
stimStopTimes = stimStartTimes + StimulusON_MS./1e3; % Note that dropped frames are not accounted for!
stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames = cellfun(@(x) x.numFrames, trialData);
repFrames = cellfun(@(x) x.repframes, trialData, 'UniformOutput',false); %numFrames./(fps*StimulusON_MS*1e-3);
repFrames(cellfun(@isempty, repFrames)) = {0};
numFrames = numFrames./ [repFrames{:}]';
numFrames(~isfinite(numFrames)) = 0;

% Get trial numbers
trialNum = cellfun(@(x) x.m_iTrialNumber, trialData);
trialNumPerFrame = cellfun(@(x, y) repelem(x, y), num2cell(trialNum), num2cell(numFrames), 'UniformOutput', false);
% Get trial IDs (insane that there are trial numbers and IDs....)

TrialID = cellfun(@(x) x.TrialID, trialData, 'UniformOutput',false);
TrialID(cellfun(@isempty, TrialID)) = {nan};
TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), TrialID, num2cell(numFrames), 'UniformOutput', false);

% usebinary: 0 = high contrast cloud, 1 = binary cloud, 2 = contrast
% matched cloud
usebinary = cellfun(@(x) x.usebinary, trialData, 'UniformOutput',false);
usebinary(cellfun(@isempty, usebinary)) = {nan};
useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), usebinary, num2cell(numFrames), 'UniformOutput', false);

% Cloud spatial scale
spatialscale = cellfun(@(x) x.spatialscale, trialData, 'UniformOutput', false);
spatialscale(cellfun(@isempty, spatialscale)) = {nan};
spatialScalePerFrame =  cellfun(@(x, y) repelem(x, y),spatialscale, num2cell(numFrames), 'UniformOutput', false);

% Stimulus blocks
BlockID = cellfun(@(x) x.BlockID, trialData, 'UniformOutput', false);
BlockID(cellfun(@isempty, BlockID)) = {nan};
BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), BlockID, num2cell(numFrames), 'UniformOutput', false);

% Matrix of unique conditions defined by cloud parameters and block
uniqueCloudConditions =  unique([[usebinary{:}]' [spatialscale{:}]' [BlockID{:}]'], 'rows');
uniqueCloudConditions(any(isnan(uniqueCloudConditions),2),:) = [];

% Get trial/stimulus types
trialType = cellfun(@(x) x.m_strTrialType, trialData, 'UniformOutput',false);
DualstimPrimaryuseRGBCloud=cellfun(@(x) x.DualstimPrimaryuseRGBCloud, trialData, 'UniformOutput',false);
DualstimPrimaryuseRGBCloud(cellfun(@isempty, DualstimPrimaryuseRGBCloud)) = {nan};
DualstimPrimaryuseRGBCloudPerFrame = cellfun(@(x, y) repelem(x, y),DualstimPrimaryuseRGBCloud, num2cell(numFrames), 'UniformOutput', false);


% Get stimulus area on each trial
StimulusArea = cellfun(@(x) x.m_aiStimulusArea, trialData, 'UniformOutput', false);
StimulusArea(cellfun(@isempty, StimulusArea)) = {nan};
StimulusAreaPerFrame =  cellfun(@(x, y) repelem(x, y),StimulusArea, num2cell(numFrames), 'UniformOutput', false);

% Stimulus rect per trials

StimulusRect = cellfun(@(x) x.m_aiStimulusRect, trialData, 'UniformOutput', false);
TiledStimulusRect = cellfun(@(x) x.m_aiTiledStimulusRect, trialData, 'UniformOutput', false);

% Get stimulus sequence on each trial
stimseq = cellfun(@(x) x.stimseq, trialData, 'UniformOutput',false);

% Get actual sequence accounting for repFrames
stimseq= cellfun(@(x,y) x(1:y:end), stimseq, repFrames, 'UniformOutput', false);

% Get fixation locations
pt2iFixationSpot = cellfun(@(x) x.m_pt2iFixationSpot, trialData, 'UniformOutput',false);
pt2iFixationSpot(cellfun(@isempty,pt2iFixationSpot)) = {[nan nan]}; % just in case
X_fixationSpot = cellfun(@(x) x(1), pt2iFixationSpot, 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), pt2iFixationSpot, 'UniformOutput',false);

% Get Kofiko determination of monkey fixation
MonkeyFixated = cellfun(@(x) logical(x.m_bMonkeyFixated), trialData);

% To get kofiko eye date for trial trialNum(t):
% Kofiko_ET_TS_PlexonTime(Kofiko_ET_TS_PlexonTime_Bin == stimON_Kofiko_ET_TS_PlexonTime_Bin(t))

% useLeye, useReye

UseLeye = cellfun(@(x) x.UseLeye, trialData, 'UniformOutput', false);
UseLeye(cellfun(@isempty, UseLeye)) = {nan};
UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), UseLeye, num2cell(numFrames), 'UniformOutput', false);

UseReye = cellfun(@(x) x.UseReye, trialData, 'UniformOutput', false);
UseReye(cellfun(@isempty, UseReye)) = {nan};
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), UseReye, num2cell(numFrames), 'UniformOutput', false);

toc;

%% Bin Kofiko eye signal timestamps by stimulus intervals (i.e., by trial)

[~,~,Kofiko_ET_TS_PlexonTime_Bin] = histcounts(Kofiko_ET_TS_PlexonTime, stimIntervals); 
%stimON_Kofiko_ET_TS_PlexonTime_Bin = uniqueKofiko_ET_TS_PlexonTime_Bin(1:2:end);

% I use accumarray to create cell array where even numbered cells give eye
% signal timestamps for stimulus on periods (stimulus onset to offset),
% and odd numbered (starting at 3) cells give eye signal timestamps between stimulus on
% periods (stimulus offset to subsequent onset).The first cell corresponds
% to invalid timestamps (i.e. outside of stimIntervals)

nBins = numel(stimIntervals) - 1;

Kofiko_ET_TS_PlexonTime_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:) + 1, ...
    Kofiko_ET_TS_PlexonTime(:), ...
    [nBins + 1, 1], ...   % force size (extra 1 for bin 0 → invalid)
    @(x){x}, ...
    {[]} ...              % fill empty bins with empty cells
);

% Do as above for eye signal X and Y coordinates:
Kofiko_Xpix_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:)+1, ...
    Kofiko_Xpix(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});

Kofiko_Ypix_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:)+1, ...
    Kofiko_Ypix(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});

%% Bin Plexon eye signal by stimulus intervals

[~,~, t_plexon_bin] = histcounts(t_plexon, stimIntervals);

t_plexon_cellArray = accumarray(t_plexon_bin(:)+1, t_plexon(:), [], @(x){x});

% pupils
leftEyePupil_plexon_cellArray = accumarray(t_plexon_bin(:)+1, leftEyePupil_plexon(:), [], @(x){x});
rightEyePupil_plexon_cellArray = accumarray(t_plexon_bin(:)+1, rightEyePupil_plexon(:), [], @(x){x});

% right eye
rightEyeX_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, rightEyeX_plexon_calib(:), [], @(x){x});
rightEyeY_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, rightEyeY_plexon_calib(:), [], @(x){x});

% left eye 
leftEyeX_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, leftEyeX_plexon_calib(:), [], @(x){x});
leftEyeY_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, leftEyeY_plexon_calib(:), [], @(x){x});

%% Determine trials with good fixation

goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < maxFixationErrorPix)./numel(x)...
    > minFixationDuration, Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < maxFixationErrorPix)./numel(y) ...
    > minFixationDuration, Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

goodFixationIdx = MonkeyFixated | (goodFixationX & goodFixationY);

dualStimIdx = strcmpi(trialType, 'Dual Stim');
trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), DualstimPrimaryuseRGBCloud);
areaOverZeroIdx = cellfun(@(x) x>0, StimulusArea);

isTrialOfInterest = dualStimIdx & ...
                      goodFixationIdx &...
                      trlonset_diffs > 4 &...
                      CCidx & ...
                      areaOverZeroIdx;