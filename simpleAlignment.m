

%% Set paths

% Add ConwayExptProcessing Dependencies to path
addpath(genpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/'));
addpath(genpath('/home/bizon/Processing'));

% "dirpath" should contain:
%   - All Kofiko mat files as well as the plexon file (which will be
%   searched for recursively)
%   - Plexon (.pl2) file
%
% "ks_path" should contain:
%   -Kilosort outputs (spike_times.npy, spike_clusters.npy, and
%   cluster_info.tsv)

dirpath = '~/Data/V1_Fovea/Jocamo/260403/';
%dirpath = '/mnt/isilon/DATA/monkey_ephys/Jocamo/2025_Singleprobe/250529';

dirpath = '/mnt/isilon/DATA/monkey_ephys/Jocamo/2022to23_ArrayExpts/';

%ks_path = '/mnt/isilon/PROJECTS/V1_Fovea/processing/250529_152043_Jacomo/';
%ks_path = '/home/bizon/Data/V1_Fovea/Jocamo/260403/260403_141736_Jacomo/kilosort_laminar_1to64';
ks_path = '/home/bizon/Data/V1_Fovea/Jocamo/220715/220715_131937_Jacomo';


stimpath = '/home/bizon/Processing/Cloudstims_calib_01_2022/'; % or 04_2024


savepath = '/home/bizon/Data/V1_Fovea/Jocamo/220715/220715_131937_Jacomo/Analysis';



% File prefix for Kofiko and plexon files
monkey_name = 'Jocamo';

filenameP = '220715_131937_Jacomo';
plexon_dir = dir(fullfile(dirpath, '**', [filenameP '.pl2']));
plexon_fname = fullfile(plexon_dir(1).folder, plexon_dir(1).name);
pl2 = PL2ReadFileIndex(plexon_fname);
%% flags
saving = 1;
compute_stas = 0;
plotting = 0;
%% Hardcoded values
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45;
minSpikes = 5000;
targ_ETstimtype = 0;
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
[~, ~, ~, arc_ad] = plx_ad_v(plexon_fname, 'AI02');
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

% create PlexET_ad structure
PlexET_times = t_plexon;
PlexET_ad_calib(:,1) = sync_ad;
PlexET_ad_calib(:,2) = arc_ad;
PlexET_ad_calib(:,3) = leftEyePupil_plexon;
PlexET_ad_calib(:,4) = rightEyePupil_plexon;
PlexET_ad_calib(:,5) = rightEyeX_plexon;
PlexET_ad_calib(:,6) = rightEyeY_plexon;
PlexET_ad_calib(:,7) = leftEyeX_plexon;
PlexET_ad_calib(:,8) = leftEyeY_plexon;

% save(

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
ImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
stimStartTimes = [ones(size(ImageFlipON_TS_Kofiko)), ImageFlipON_TS_Kofiko]*B; % plexon time
stimStopTimes = stimStartTimes + StimulusON_MS./1e3; % Note that dropped frames are not accounted for!
stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames = cellfun(@(x) x.numFrames, trialData);
repFrames = cellfun(@(x) x.repframes, trialData, 'UniformOutput',false);
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

%% Extract fixinfo

calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), trialType);

ETdata.fixloc = vertcat(pt2iFixationSpot{isCalibrationTrial});
ETdata.plxonset = stimStartTimes(isCalibrationTrial);

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

%% Bin stimulus sequences by trial
tic;
fprintf('Loading stimuli\n')
stimulus_cellArray = cell(size(trialData,1),1);

for i = 1:size(uniqueCloudConditions,1)
    thisUsebinary = uniqueCloudConditions(i,1);
    thisSpatialscale = uniqueCloudConditions(i,2);
    thisBlockID = uniqueCloudConditions(i,3);

    trialIdx = [usebinary{:}]' == thisUsebinary & [spatialscale{:}]' == thisSpatialscale & [BlockID{:}]' == thisBlockID;
    if thisUsebinary == 0 % full contrast
        load(fullfile(stimpath, sprintf('Cloudstims_Chrom_size60_scale%d_%02d.mat', thisSpatialscale, thisBlockID)));
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
    elseif thisUsebinary == 1
        load(fullfile(stimpath,sprintf('Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat', thisSpatialscale, thisBlockID)))
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
    elseif thisUsebinary == 2 % matched contrast
        load(fullfile(stimpath,sprintf('Cloudstims_ContrastMatched_size60_scale%d_%02d.mat', thisSpatialscale, thisBlockID)));
        DensenoiseChromcloud_DKlspace(:,:,:,2:3)=int8(127*(DensenoiseChromcloud_DKlspace(:,:,:,2:3)));
        DensenoiseChromcloud_DKlspace(:,:,:,1)=int8((1/LumScale)*127*(DensenoiseChromcloud_DKlspace(:,:,:,1)));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
    end

    % reshape stimuli into matrix where each row is stimulus, each col is
    % framez

    % make first dim frames

    stimulus_cellArray(trialIdx) = cellfun(@(x) permute(x, [3 1 2 4]), stimulus_cellArray(trialIdx), 'UniformOutput', false);
    stimulus_cellArray(trialIdx) = cellfun(@(x) reshape(x, [size(x,1), prod(size(x,2:4))]),stimulus_cellArray(trialIdx), 'UniformOutput', false);
    stimulus_cellArray(trialIdx) = cellfun(@transpose, stimulus_cellArray(trialIdx), 'UniformOutput', false);
end

% find dualstim elemetns

trialTypeOfInterestIdx = strcmpi(trialType, 'Dual Stim');
trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), DualstimPrimaryuseRGBCloud);
areaOverZeroIdx = cellfun(@(x) x>0, StimulusArea);

isTrialOfInterest = trialTypeOfInterestIdx & ...
    goodFixationIdx &...
    trlonset_diffs > 4 &...
    CCidx & ...
    areaOverZeroIdx;

stimulus_matrix = horzcat(stimulus_cellArray{isTrialOfInterest});

stimFrameBins = cellfun(@(x, y, z) linspace(x,y,z+1), num2cell(stimStartTimes), num2cell(stimStopTimes), num2cell(numFrames), 'UniformOutput', false);
valid_stimFrameBinIdx = cellfun(@(x) numel(x)>=2, stimFrameBins);
toc;

%% resample eye signal at frame rate

% NOTE: Without 'extrap' argument, there will be NANs:

Kofiko_Xpix_frameRate_cellArray = ...
    cellfun(@(t,x, start, stop, fixspot, n) interp1(t, x - fixspot, linspace(start, stop, n), 'linear'),...
    Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    Kofiko_Xpix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimStartTimes(isTrialOfInterest)),...
    num2cell(stimStopTimes(isTrialOfInterest)),...
    X_fixationSpot(isTrialOfInterest),...
    num2cell(numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);

Kofiko_Ypix_frameRate_cellArray = ...
    cellfun(@(t,y, start, stop, fixspot, n) interp1(t, y - fixspot, linspace(start, stop, n), 'linear'),...
    Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    Kofiko_Ypix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimStartTimes(isTrialOfInterest)),...
    num2cell(stimStopTimes(isTrialOfInterest)),...
    Y_fixationSpot(isTrialOfInterest),...
    num2cell(numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);
%% load in spike data
tic;
fprintf('Loaidng spike data and computing Robs\n');
% drop test check

numDigitsInLastSpkChan = ceil(log10(length(pl2.SpikeChannels)));

numAIchans = sum(cellfun(@(x) contains(x.Name, 'AI'), pl2.AnalogChannels));
numDigitsInLastAIchan = ceil(log10(numAIchans));
disp('Drop test check starting')

[fs, n, ts, fn, ~] = plx_ad_v(plexon_fname,...
    ['SPKC' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);

[fs_aux, n_aux, ts_aux, fn_aux, ~] = plx_ad_v(plexon_fname, ['AI' num2str(1, ['%0' num2str(numDigitsInLastAIchan), '.f'])] );

spikeChannel1SignalDurSec = n/fs; % samples / (samples/sec)
dpiSyncSignalDurSec = n_aux/fs_aux;

droptestcheck = n_aux/fs_aux - n/fs;
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present

if abs(droptestcheck)>0.1
    warning("Danger - Plexon might have dropped frames! Check pl2 file.")
    opts.spk_offset = droptestcheck;
else
    opts.spk_offset = 0;
end

disp('Drop test check complete')
% end drop test check %%%%
%% cluster label strings

goodStr = 'good ';
muaStr = 'mua  ';
noiseStr = 'noise';
blankStr = '     ';
%% Load and organize spike data

% Find folders with kilosort output
spike_times_dir = dir(fullfile(ks_path, '**/spike_times.npy'));
spike_clusters_dir = dir(fullfile(ks_path, '**/spike_clusters.npy'));
cluster_info_dir = dir(fullfile(ks_path, '**/cluster_info.tsv'));


spike_times_folders = {spike_times_dir(:).folder};
spike_clusters_folders = {spike_clusters_dir(:).folder};
cluster_info_folders = {cluster_info_dir(:).folder};



% For each folder with kilosort outputs
num_ks_batch = length(spike_times_dir);

[~, ks_folders, ~] =  cellfun(@fileparts, spike_times_folders, 'UniformOutput', false);
tokens = regexp(ks_folders, '^[^_]+_([^_]+)', 'tokens');
array_labels = cellfun(@(t) t{1}{1}, tokens, 'UniformOutput', false);
unique_array_labels = unique(array_labels);
chan_offset = 0;
cluster_offset = 0;

chan_offsets{1} = chan_offset;



for ks_batch = 1:num_ks_batch
    this_array_label = array_labels{ks_batch};

    try
        % Get label of kilosort batch (often corresponding to array name and range of channels processed)

        % Read in kilosort outputs
        spk_times = readNPY(fullfile(spike_times_folders{ks_batch}, 'spike_times.npy'));
        spk_times = double(spk_times)./fs + opts.spk_offset; % convert to seconds
        spk_clusters = readNPY(fullfile(spike_times_folders{ks_batch}, 'spike_clusters.npy'))...
            + cluster_offset;
        cluster_KSLabel = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_KSLabel.tsv'));

        chan_map = readNPY(fullfile(spike_times_folders{ks_batch}, 'channel_map.npy'));
        chan_map = chan_map + chan_offset; % make channel numbers unique within array

        if isfile(fullfile(spike_times_folders{ks_batch}, 'cluster_info.tsv'))
            cluster_info = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_info.tsv'));
            cluster_id = cluster_info.cluster_id + cluster_offset;
            group = cluster_info.group;
            n_spikes = cluster_info.n_spikes;
            chan_best = cluster_info.ch + double(chan_offset);

        else
            cluster_group = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_group.tsv'));
            cluster_group.cluster_id = cluster_group.cluster_id + cluster_offset;

            
            % account for blank units which may not be in cluster_group
            blank_cluster_id = setdiff(unique(spk_clusters), cluster_group.cluster_id);
            temp_cluster_id = [cluster_group.cluster_id; blank_cluster_id];
            temp_group = cluster_group.group;
            temp_group(end+1:end+length(blank_cluster_id),:) = ' ';
            [temp_cluster_id_sorted, I] = sort(temp_cluster_id);
            temp_group_sorted = temp_group(I,:);
            cluster_id = temp_cluster_id_sorted;
            group = temp_group_sorted;

            % find best channel of each cluster
            
            templates = readNPY(fullfile(spike_times_folders{ks_batch}, 'templates.npy'));
            n_spikes = accumarray(spk_clusters+1, spk_clusters, [], @numel);
            n_spikes = n_spikes(n_spikes>0);
            [~,I]= max(sum(templates.^2,2),[],3);
            chan_best = chan_map(I); % best channel for each unique cluster

        end


        % Find indices of units labeled "good", "mua", or ""
        isGood = cellfun(@(x) strcmpi(deblank(x), 'good'), cellstr(group));
        isMua = cellfun(@(x) strcmpi(deblank(x), 'mua'), cellstr(group));
        isBlank = cellfun(@(x) isempty(deblank(x)), cellstr(group));
        hasMinSpikes = n_spikes > minSpikes;

        % get cluster IDs
        SU_clusterIDs{ks_batch} = cluster_id(isGood & hasMinSpikes);
        MU_clusterIDs{ks_batch} = cluster_id((isMua | isBlank) & hasMinSpikes);
        allUnit_clusterIDs{ks_batch} =[SU_clusterIDs{ks_batch}; MU_clusterIDs{ks_batch}];

        % get channels
        SU_chanNums{ks_batch} = chan_best(isGood & hasMinSpikes);
        MU_chanNums{ks_batch} = chan_best((isMua | isBlank) & hasMinSpikes);
        allUnit_chanNums{ks_batch} =  [SU_chanNums{ks_batch}; MU_chanNums{ks_batch}];

        % get ks_batch
        SU_ks_batch{ks_batch} = ks_batch.*ones(size(SU_clusterIDs{ks_batch}));
        MU_ks_batch{ks_batch} = ks_batch.*ones(size(MU_clusterIDs{ks_batch}));
        allUnit_ks_batch{ks_batch} = [SU_ks_batch{ks_batch}; MU_ks_batch{ks_batch}];

        % Get rid of spike times and clusterIDs that correspond to bad
        % units
        %
        spk_times = spk_times(ismember(spk_clusters, allUnit_clusterIDs{ks_batch}));
        spk_clusters = spk_clusters(ismember(spk_clusters, allUnit_clusterIDs{ks_batch}));

        ks_batchForEachSpk= repelem(ks_batch, numel(spk_times))';

        % bin spikes by trial
        [~,~,spk_times_bin{ks_batch}] = histcounts(spk_times, stimIntervals);

        % each cell gives spike times for stim ON and stim OFF periods

        spk_times_cellArray{ks_batch} =accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            spk_times(:), ...
            [nBins + 1, 1], ...
            @(x){x}, ...
            {[]});

        clusterIDForEachSpk_cellArray{ks_batch} = accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            spk_clusters(:), ...
            [nBins + 1, 1], ...
            @(x){x}, ...
            {[]});

        ks_batchForEachSpk_cellArray{ks_batch} = accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            ks_batchForEachSpk(:), ...
            [nBins + 1, 1], ...)
            @(x){x}, ...
            {[]});

        clusterIDForEachSpk_stimON_cellArray = clusterIDForEachSpk_cellArray{ks_batch}(2:2:end);
        spk_times_stimON_cellArray =  spk_times_cellArray{ks_batch}(2:2:end);
        %    ks_batchForEachSpk_stimON_cellArray = ks_batchForEachSpk_stimON_cellArray{ks_batch}(2:2:end);

        stimFrameNumForEachSpk_cellArray{ks_batch} = cell(size(trialData,1),1);

        [~,~,stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx)] =...
            cellfun(@(x,y) histcounts(x, 'BinEdges', y), ...
            spk_times_stimON_cellArray(valid_stimFrameBinIdx),...
            stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
        %spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

        Robs{ks_batch} = single(zeros(numel(allUnit_clusterIDs{ks_batch}),  sum(numFrames(isTrialOfInterest))));

        for unit = 1: numel(allUnit_clusterIDs{ks_batch})
            unitID = allUnit_clusterIDs{ks_batch}(unit);

            % each cell is a trial, and contains frame indicies where spikes
            % occurred (repeated frames == multiple spikes on that frame)

            stimFrameNumForEachSpk_thisUnit_cellArray = cell(size(trialData,1),1);
            stimFrameNumForEachSpk_thisUnit_cellArray(valid_stimFrameBinIdx) = ...
                cellfun(@(x,y) x(y==unitID), ...
                stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx), ...
                clusterIDForEachSpk_stimON_cellArray(valid_stimFrameBinIdx),...
                'UniformOutput', false);

            spksPerFrame_cellArray =  cell(size(trialData,1),1);

            spksPerFrame_cellArray(isTrialOfInterest) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
                stimFrameNumForEachSpk_thisUnit_cellArray(isTrialOfInterest), num2cell(numFrames(isTrialOfInterest)) ,'UniformOutput',false);

            Robs{ks_batch}(unit,:) = [spksPerFrame_cellArray{:}];

        end

    catch ME
        warning("Could not process %s", ks_folders{ks_batch})
        disp(ME.message)
    end
    % update chan_offset
    if ks_batch < num_ks_batch % if we still got a batch ahead
        next_array_label = array_labels{ks_batch+1};
        if strcmpi(next_array_label,  this_array_label)
            chan_offset = max(chan_map)+1; % maximum channel number;
        else
            chan_offset = 0; % new array, reset offset to 0
        end
        chan_offsets{ks_batch+1} = chan_offset;
    end
    % update cluster offset, and keep track of them
    cluster_offset = max(cluster_id)+1;
    cluster_offsets{ks_batch} = cluster_offset;

end

nSU = sum(cellfun(@length, SU_clusterIDs));
nMU = sum(cellfun(@length, MU_clusterIDs));

toc;

%% Format output like PackageCloud

allRobs = vertcat(Robs{:});

%ETgains
ETgains = [Kofiko_GainX(end), Kofiko_GainY(end)];

%ETstim_location
ETstim_location =  [trialData{find(isTrialOfInterest,1, 'last')}.secondarystim_bar_rect;...
    trialData{find(isTrialOfInterest,1, 'last')}.tertiarystim_bar_rect];

% ETtrace
ETtrace = [[Kofiko_Xpix_frameRate_cellArray{:}]; [Kofiko_Ypix_frameRate_cellArray{:}]];

%ETtrace_raw

ETtrace_raw = transpose([vertcat(rightEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyePupil_plexon_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyePupil_plexon_cellArray{2*find(isTrialOfInterest)})]);

%Robs
RobsSU = allRobs(1:nSU,:);

%RobsMU
RobsMU = allRobs(nSU+1:end,:);

%RobsMU_probe_ID
RobsMU_probe_ID = vertcat(MU_chanNums{:});

%RobsMU_rating
RobsMU_rating = [];
%blockID

% Robs_probe_ID
Robs_probe_ID = vertcat(SU_chanNums{:});

%Robs rating

Robs_rating =[];

blockID = [BlockIDPerFrame{isTrialOfInterest}];

%block_inds
block_inds = [cumsum(numFrames(isTrialOfInterest))' - (numFrames(find(isTrialOfInterest,1,'first'))-1); cumsum(numFrames(isTrialOfInterest))'];

%cloud_area
cloud_area = [StimulusAreaPerFrame{isTrialOfInterest}];
%cloud_binary
cloud_binary = [useBinaryPerFrame{isTrialOfInterest}];

%cloud_scale
cloud_scale = [spatialScalePerFrame{isTrialOfInterest}];

%datafilts

datafilts = ones(size(RobsSU));

%datafiltsMU
datafiltsMU = ones(size(RobsMU));

%dt
dt = 0.0160; % why not 0.0167?

%electrode_info
electrode_info =[];

%exptdate
exptdate = filenameP(1:6);
%exptname
exptname = filenameP;

%fix_location
fix_location = unique(vertcat(pt2iFixationSpot{isTrialOfInterest}), 'rows');

%fix_size
fix_size = trialData{find(isTrialOfInterest, 1, 'last')}.m_fFixationSizePix -1;

%pixel_size
pixel_size = 1;

%sacc_inds
sacc_inds=[];

%ss
clusterIDs= cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), clusterIDForEachSpk_cellArray, 'UniformOutput', false);
clusterIDs = transpose(vertcat(clusterIDs{:}));
clusterIDs = vertcat(clusterIDs);

assert(numel(unique(clusterIDs)) == nSU + nMU);

%spike_ts
spike_ts = cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), spk_times_cellArray, 'UniformOutput', false);
spike_ts = transpose(vertcat(spike_ts{:}));
spike_ts = vertcat(spike_ts);

%stim
stim  = reshape(stimulus_matrix, 60,60,3,[]);

%stimET
stimET = [];
%stim_area

stim_area = [StimulusArea{isTrialOfInterest}]';

%stim_location
modal_stim_area = mode(stim_area);
medianTiledStimulusRect =  median(cat(3, TiledStimulusRect{:}),3);

stim_location =  [medianTiledStimulusRect(:,1:2) medianTiledStimulusRect(:,1:2)+60];

%stim_location_deltas: first two elements of first row of TiledStimulusRect
%minus median stim_location x and y

stim_location_deltas = cellfun(@(x) x(1,1:2) - stim_location(1,1:2), TiledStimulusRect(isTrialOfInterest), 'UniformOutput',false);
stim_location_deltas = vertcat(stim_location_deltas{:});

%stimscale
stimscale = (stim_location(3)-stim_location(1))/60;

%stimtype
stimtype = [DualstimPrimaryuseRGBCloudPerFrame{isTrialOfInterest}];

%trialID
trialID = [TrialIDPerFrame{isTrialOfInterest}];

%trial_start_ts
trial_start_ts = [stimStartTimes(isTrialOfInterest)]';

%useLeye
useLeye = [UseLeyePerFrame{isTrialOfInterest}];

%useReye
useReye = [UseReyePerFrame{isTrialOfInterest}];

%valid_data
totalFrames = sum(numFrames(isTrialOfInterest));
tvec=1:totalFrames;
Block_offsetinds = block_inds(2,:);

bad_inds_block=sort([block_inds(1,:), block_inds(1,:)+1, block_inds(1,:)+2,block_inds(1,:)+3,block_inds(1,:)+4,block_inds(1,:)+5,block_inds(1,:)+6],1);

ETdist_thresh=40;
bad_inds_fix = unique([find(abs(ETtrace(1,:))>ETdist_thresh),find(abs(ETtrace(2,:))>ETdist_thresh)]);
bad_inds_all = unique([bad_inds_block,bad_inds_fix,bad_inds_fix-1, bad_inds_fix-2, bad_inds_fix+1, bad_inds_fix+2]); % remove indices immediately preceding and following eye movement artifacts
use_inds_fix=setdiff(tvec,bad_inds_all);

% find sequences less than 10 due to eye movement removal, and exclude them
% to avoid clogging up the modeling pipeline with tiny snippets
diffs=diff([1,use_inds_fix]);
[~,X]=find(diff(diffs)<10);
for k= X, use_inds_fix(diffs(k):diffs(k+1)-1)=0; end
use_inds_fix(use_inds_fix==0)=[];

valid_data = use_inds_fix;

% added by me: ks batch

spikeSortingBatch = vertcat(SU_ks_batch{:});
spikeSortingBatchMU = vertcat(MU_ks_batch{:});

cluster = vertcat(SU_clusterIDs{:});
clusterMU = vertcat(MU_clusterIDs{:});

%% reorganize cluster ids by array
% make cell array identical to pks_times_cellArray where each element is ks
% folder number

[clusterIDs_sorted,sortByClusterID] = sort(clusterIDs);
spike_ts_sorted = spike_ts(sortByClusterID);

spikeIDs = clusterIDs_sorted;

%% add fields to data struct
data.ETgains = ETgains;
data.ETstim_location = ETstim_location;
data.ETtrace = ETtrace;
data.ETtrace_raw = ETtrace_raw;
data.Robs = RobsSU;
data.RobsMU = RobsMU;
data.RobsMU_probe_ID = RobsMU_probe_ID;
data.RobsMU_rating = RobsMU_rating;
data.Robs_probe_ID = Robs_probe_ID;
data.Robs_rating = [];
data.blockID = blockID;
data.block_inds = block_inds;
data.cloud_area = cloud_area;
data.cloud_binary = cloud_binary;
data.cloud_scale = cloud_scale;
data.datafilts = datafilts;
data.datafiltsMU = datafiltsMU;
data.dt = dt;
data.electrode_info = electrode_info;
data.exptdate = exptdate;
data.exptname = exptname;
data.fix_location = fix_location;
data.fix_size = fix_size;
data.pixel_size = pixel_size;
data.sacc_inds = sacc_inds;
data.spikeIDs = spikeIDs;
data.spike_ts = spike_ts_sorted;
data.stim = stim;
data.stimET = stimET;
data.stim_area = stim_area;
data.stim_location = stim_location;
data.stim_location_deltas = stim_location_deltas;
data.stimscale = stimscale;
data.stimtype = stimtype;
data.trialID = trialID;
data.trial_start_ts = trial_start_ts;
data.useLeye = useLeye;
data.useReye = useReye;
data.valid_data = valid_data;

data.spikeSortingBatch = spikeSortingBatch;
data.spikeSortingBatchMU = spikeSortingBatchMU;

data.ks_folders = ks_folders;
data.chan_offsets = chan_offsets;
data.cluster_offsets = cluster_offsets;
data.array_labels = array_labels;

[C, IA, IC] = unique(array_labels);
data.arrayPerSU = IC(data.spikeSortingBatch);
data.arrayPerMU = IC(data.spikeSortingBatchMU);

data.cluster = cluster;
data.clusterMU = clusterMU;

%% Compute STAs (optional)
if compute_stas
    tic;
    fprintf('Computing STAs\n');
    STA = cell(size(Robs));
    for ks_batch = 1:length(Robs)
        for unit = 1:size(Robs{ks_batch},1)
            r = transpose(Robs{ks_batch }(unit,:));
            spkTimes = find(r);
            nSpks = numel(spkTimes);
            nLags = 10;
            maxLag = nLags - 1;
            tempSTA = zeros(size(stimulus_matrix,1), nLags, 'like', stimulus_matrix);
            countsPerLag = zeros(1,nLags);
            for spk = 1:nSpks
                t = spkTimes(spk);
                if t <= maxLag
                    continue;
                end
                cols = (t - (0:maxLag));
                tempSTA = single(tempSTA) + (single(stimulus_matrix(:,cols))./127)*r(t);
                countsPerLag = countsPerLag + r(t);
            end
            countsPerLag(countsPerLag==0)=1;
            tempSTA = bsxfun(@rdivide, single(tempSTA), single(countsPerLag));
            STA{ks_batch}(unit,:,:) = tempSTA;
        end
        STA{ks_batch} = reshape(STA{ks_batch}, size(STA{ks_batch},1), 60, 60, 3, nLags);

    end
    toc;
else
end
%% Plotting

lags = 2:5;
if compute_stas && plotting
    fprintf('Now plotting\n')
    for ks_batch = 1:length(STA)
        for unit = 1:size(STA{ks_batch},1)
            figure; i = 1;
            for chromatic_channel = 1:3
                for lag = lags
                    subplot(3, length(lags), i );
                    imagesc(circshift(squeeze(STA{ks_batch}(unit,:,:,chromatic_channel,lag+1)), [30 30])), colormap gray
                    i = i+1;
                end
            end

            unitChanNum = allUnit_chanNums{ks_batch}(unit);
            unitClusterID = allUnit_clusterIDs{ks_batch}(unit);

            sgtitle(sprintf('%s, channel: %i, clusterID: %i',  replace(ks_folders{ks_batch}, '_', ' '),unitChanNum, unitClusterID))
        end
    end
end

%% saving
if saving

    if ~isdir(savepath)
        mkdir(savepath);
    end

    switch targ_ETstimtype
        case 0; curETstimtype = 'NA';
        case 1; curETstimtype='1D';
        case 7; curETstimtype='CC';
    end


    switch_stimtype = unique(vertcat(DualstimPrimaryuseRGBCloud{isTrialOfInterest}));
    switch switch_stimtype
        case 0; curstimstype='GT';
        case 3; curstimstype='HL';
        case 6; curstimstype='HC';
        case 8; curstimstype='CC';
    end

    unique_array_labels = unique(array_labels);
    array_label_filepart = [cellfun(@(x) [x '_'], unique_array_labels(1:end-1), 'UniformOutput', false) unique_array_labels(end)];
    array_label_filepart = horzcat(array_label_filepart{:});

    useofflinesorting =1;
    FullExpt_ET_filename = sprintf( '%s_FullExpt_ET.mat', filenameP );
    data_filename=[monkey_name '_' exptname(1:6) '_' array_label_filepart '_' curstimstype '_ET' curETstimtype '_v09.mat'];
    fixinfo_filename=[filenameP '_fixinfo.mat'];
    save(fullfile(savepath, data_filename),  '-struct', 'data', '-v7.3'); % save packaged cloud data
    save(fullfile(savepath, fixinfo_filename), '-struct', 'ETdata', '-v7.3') % save fixinfo
    save(fullfile(savepath, FullExpt_ET_filename), 'PlexET_ad_calib', 'PlexET_times', '-v7.3'); % save FullExpt_ET

end

