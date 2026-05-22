%% Set paths

function FR = FRfromCalibrationTrials


% Add ConwayExptProcessing Dependencies to path
% whats needed from this are functions to read npy and tsv files

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
% "stimpath" should contain:
%   -Mat files with cloud stimuli

dirpath = '/mnt/isilon/DATA/monkey_ephys/Jocamo/2025_Singleprobe/250529';
ks_path ='/mnt/isilon/PROJECTS/V1_Fovea/processing/250529_152043_Jacomo';
stimpath = '/home/bizon/Processing/Cloudstims_calib_04_2024'; % or 04_2024
savepath = '' ;%'/home/bizon/Data/V1_Fovea/Sprout/260505/260505_150945_Sprout/Analysis';%' ;%'/mnt/isilon/users/greenemj/V1_Fovea/Jocamo/220715/Analysis';

rig = 'C';

%'/home/bizon/Data/V1_Fovea/Jocamo/250529_152043_Jacomo/Analysis';%
% save time with these assertions
assert(isdir(dirpath) & isdir(ks_path) & isdir(stimpath), 'Check that dirpath, ks_path, and stimpath exist!');

% File prefix for Kofiko and plexon files
monkey_name = 'Jocamo';
filenameP = '250529_152043_Jacomo';%';
plexon_dir = dir(fullfile(dirpath, '**', [filenameP '.pl2']));
plexon_fname = fullfile(plexon_dir(1).folder, plexon_dir(1).name);
pl2 = PL2ReadFileIndex(plexon_fname);
%% flags
saving = 0;
compute_stas = 0;
plotting = 0;
%% Hardcoded values
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45;
minSpikes = 2000;
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
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber ==  g_strctSystemCodes.m_iSync);
% Linear regression to get kofiko time stamps into plexon time
B = [ones(size(kofikoSyncStrobesTS)) kofikoSyncStrobesTS]\plexonSyncStrobesTS;

%% Plexon eye data


if strcmpi(rig, 'B')
    % Rig B
    sync_ch = 'AI01';
    arc_ch = 'AI02';
    leftEyePupil_ch = 'AI07';
    rightEyePupil_ch = 'AI08';
    rightEyeX_ch = 'AI05';
    rightEyeY_ch = 'AI06';
    leftEyeX_ch ='AI03';
    leftEyeY_ch ='AI04';
elseif strcmpi(rig, 'C')
    % % Rig C
    sync_ch = 'AI01';
    arc_ch = 'AI02';
    leftEyePupil_ch = 'AI03';
    rightEyePupil_ch = 'AI04';
    rightEyeX_ch = 'AI05';
    rightEyeY_ch = 'AI06';
    leftEyeX_ch ='AI07';
    leftEyeY_ch ='AI08';

end


[adfreq, n, ts, fn, sync_ad] = plx_ad_v(plexon_fname, sync_ch);
[~, ~, ~, arc_ad] = plx_ad_v(plexon_fname, arc_ch);
[~, ~, ~, ~, leftEyePupil_plexon] = plx_ad_v(plexon_fname, leftEyePupil_ch);
[~, ~, ~, ~, rightEyePupil_plexon] = plx_ad_v(plexon_fname, rightEyePupil_ch);
[~, ~, ~, ~, rightEyeX_plexon] = plx_ad_v(plexon_fname,rightEyeX_ch);
[~, ~, ~, ~, rightEyeY_plexon] = plx_ad_v(plexon_fname, rightEyeY_ch);
[~, ~, ~, ~, leftEyeX_plexon] = plx_ad_v(plexon_fname, leftEyeX_ch);
[~, ~, ~, ~, leftEyeY_plexon] = plx_ad_v(plexon_fname, leftEyeY_ch);

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
PlexET_ad_calib = [];
PlexET_ad_calib(:,1) = sync_ad;
PlexET_ad_calib(:,2) = arc_ad;
PlexET_ad_calib(:,3) = leftEyePupil_plexon;
PlexET_ad_calib(:,4) = rightEyePupil_plexon;
PlexET_ad_calib(:,5) = rightEyeX_plexon;
PlexET_ad_calib(:,6) = rightEyeY_plexon;
PlexET_ad_calib(:,7) = leftEyeX_plexon;
PlexET_ad_calib(:,8) = leftEyeY_plexon;

% save(

%% Trial analysis
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
% spatialscale = cellfun(@(x) x.spatialscale, trialData, 'UniformOutput', false);
% spatialscale(cellfun(@isempty, spatialscale)) = {nan};

try
spatialscale_buffer = g_astrctAllParadigms{1}.DualstimScale.Buffer;
spatialscale_TS = g_astrctAllParadigms{1}.DualstimScale.TimeStamp;
spatialscale = nan(size(trialData));
for i = 1:numel(spatialscale_buffer)
    spatialscale(ImageFlipON_TS_Kofiko >= spatialscale_TS(i)) = spatialscale_buffer(i);
end
catch
    spatialscale = cellfun(@(x) x.spatialscale, trialData, 'UniformOutput', false);
    spatialscale(cellfun(@isempty, spatialscale)) = {nan};
    spatialscale = vertcat(spatialscale{:});

end


%spatialScalePerFrame =  cellfun(@(x, y) repelem(x, y),spatialscale, num2cell(numFrames), 'UniformOutput', false);


% Stimulus blocks
BlockID = cellfun(@(x) x.BlockID, trialData, 'UniformOutput', false);
BlockID(cellfun(@isempty, BlockID)) = {nan};
BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), BlockID, num2cell(numFrames), 'UniformOutput', false);

% Matrix of unique conditions defined by cloud parameters and block
uniqueCloudConditions =  unique([[usebinary{:}]' spatialscale [BlockID{:}]'], 'rows');
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



%% Bin stimulus sequences by trial

isTrialOfInterest = ismember(trialType,calibrationTasks);

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
fprintf('Loaidng spike data and computing FR\n');
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
assert(~isempty(spike_times_dir));

[~, ks_folders, ~] =  cellfun(@fileparts, spike_times_folders, 'UniformOutput', false);
tokens = regexp(ks_folders, '^[^_]+_([^_]+)', 'tokens');
array_labels = cellfun(@(t) t{1}{1}, tokens, 'UniformOutput', false);
unique_array_labels = unique(array_labels);
chan_offset = 0;
cluster_offset = 0;

%chan_offsets(1) = chan_offset;

for ks_batch = 1:num_ks_batch
    this_array_label = array_labels{ks_batch};
    disp(ks_batch)
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
            % n_spikes = accumarray(spk_clusters+1, spk_clusters, [], @numel);
            % n_spikes = n_spikes(n_spikes>0);

            n_spikes = accumarray(spk_clusters+1, 1, [], @sum);
            n_spikes = n_spikes(cluster_id + 1);
            
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

        

        FR{ks_batch} = single(zeros(numel(allUnit_clusterIDs{ks_batch}),  1));

        for unit = 1: numel(allUnit_clusterIDs{ks_batch})
            unitID = allUnit_clusterIDs{ks_batch}(unit);

            FR{ks_batch}(unit,:) = sum(spk_clusters==unitID);

        end

    catch ME
        warning("Could not process %s", ks_folders{ks_batch})
        disp(ME.message)
    end

    %%%% Added 5/9/26

    nSU{ks_batch} = length(SU_clusterIDs{ks_batch});
    nMU{ks_batch} = length(MU_clusterIDs{ks_batch});
    FRSU{ks_batch} = FR{ks_batch}(1:nSU{ks_batch},:);
    FRMU{ks_batch} = FR{ks_batch}(nSU{ks_batch}+1:nSU{ks_batch}+nMU{ks_batch},:);
    %%%%

    % update chan_offset
    if ks_batch < num_ks_batch % if we still got a batch ahead
        next_array_label = array_labels{ks_batch+1};
        if strcmpi(next_array_label,  this_array_label)
            chan_offset = max(chan_map)+1; % maximum channel number;
        else
            chan_offset = 0; % new array, reset offset to 0
        end
     %   chan_offsets(ks_batch+1) = chan_offset;
    end
    % update cluster offset, and keep track of them
    cluster_offset = max(cluster_id)+1;
    cluster_offsets(ks_batch+1) = cluster_offset;
    chan_offsets(ks_batch+1) = chan_offset;

end

nSU = sum(cellfun(@length, SU_clusterIDs));
nMU = sum(cellfun(@length, MU_clusterIDs));

toc;

%% Format output like PackageCloud

%allFR = vertcat(FR{:});

%%%% Added 5/9/26
allFRSU = vertcat(FRSU{:});
allFRMU = vertcat(FRMU{:});
allFR = vertcat(allFRSU, allFRMU);
%%%%

%ETgains
ETgains = [Kofiko_GainX(end), Kofiko_GainY(end)];


% ETtrace
ETtrace = [[Kofiko_Xpix_frameRate_cellArray{:}]; [Kofiko_Ypix_frameRate_cellArray{:}]];

%ETtrace_raw

ETtrace_raw = transpose([vertcat(rightEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyePupil_plexon_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyePupil_plexon_cellArray{2*find(isTrialOfInterest)})]);

%FR
FRSU = allFR(1:nSU,:);

%FRMU
FRMU = allFR(nSU+1:end,:);
allFR = allFR/t_plexon(end);


