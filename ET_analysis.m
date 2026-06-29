kofiko_subfolder = '/home/bizon/Data/V1_Fovea/Sprout/260601/260601_150343_Sprout';
plexon_fname = [kofiko_subfolder '.pl2']; %/home/bizon/Data/V1_Fovea/Sprout/260505/260505_150945_Sprout.pl2';
mainKofiko_fname = [kofiko_subfolder '.mat'];% '/home/bizon/Data/V1_Fovea/Sprout/260505/260505_150945_Sprout.mat';
%dpi_fname = '/home/bizon/Data/V1_Fovea/Sprout/260601/RawDDPI-2026Jun01-150214/RawDDPI-2026Jun01-150214.txt';
[~, filenameP] = fileparts(mainKofiko_fname);

pl2 = PL2ReadFileIndex(plexon_fname);

rig = 'C';
%% Hardcoded values
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45;
minSpikes = 2000;
targ_ETstimtype = 0;
ETdist_thresh=40;

%% Load kofiko data
tic;
fprintf('Loading Kofiko trial data\n')

computerLocation = 1;
if computerLocation < 10
    kofiko_dir = dir(fullfile(kofiko_subfolder, [filenameP '*.mat']));
    kofiko_fname_re = [filenameP '(_\d*)?.mat'];
    valid_kofiko_fname_idx = cellfun(@(x) ~isempty(x), ...
        cellfun(@(x) regexp(x, kofiko_fname_re), {kofiko_dir(:).name}, 'UniformOutput',false));
    kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
    kofiko_fname = {kofiko_dir(valid_kofiko_fname_idx).name};

    % Get the Kofiko mat file numbers from their names
    kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'),...
        kofiko_fname, 'UniformOutput',false);

    kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);
else
    kofiko_dir = dir(fullfile(ks_path, '*.mat'));
    all_names = {kofiko_dir.name};
    kofiko_fname_re = ['^' filenameP '(?:_\d+)?\.mat$'];
    matches = regexp(all_names, kofiko_fname_re, 'once');
    valid_kofiko_fname_idx = ~cellfun(@isempty, matches);

    kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
    kofiko_fname  = {kofiko_dir(valid_kofiko_fname_idx).name};

    % Get the Kofiko mat file numbers from their names
    kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'), kofiko_fname, 'UniformOutput',false);
    if sum(cellfun(@isempty, kofiko_file_numStr)) > 0
    	kofiko_file_numStr{cellfun(@isempty, kofiko_file_numStr)} = {'nan'};
    end
    kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);
    kofiko_file_num(isnan(kofiko_file_num)) = max(kofiko_file_num)+1;

end



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

% add main kofiko file

temp_strct = load(mainKofiko_fname, 'g_strctLocalExperimentRecording');
temp_cell = temp_strct.g_strctLocalExperimentRecording;
g_strctLocalExperimentRecording{end+1} = temp_cell;

% Load the file which should be the "main" Kofiko mat file (no number at the end)

%mainKofikoFolderIdx = strcmpi(kofiko_fname, [filenameP '.mat']);
%load(fullfile(kofiko_folder{mainKofikoFolderIdx}, kofiko_fname{mainKofikoFolderIdx}), 'g_strctDAQParams', ...
load(mainKofiko_fname, 'g_strctDAQParams', ...
    'g_astrctAllParadigms', 'g_strctLog', 'g_strctEyeCalib', ...
    'g_strctAppConfig', 'g_strctSharedParadigmData', ...
    'g_strctStimulusServer', 'g_strctSystemCodes');

% Concatenate the g_strctLocalExperimentRecorindg structs and remove empty cells
trialData = vertcat(g_strctLocalExperimentRecording{:});
trialData(cellfun(@isempty, trialData)) = [];

% Unify field names across g_strctLocalExperimentRecording to facilitate using cellfun
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



%% %%%%%%%% Time alignment %%%%%%%%%%
tic;
fprintf('Loading Kofiko eye signals and converting to plexon time\n')

% Read kofiko strobes from plexon and kofiko files
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);

kofikoSyncStrobesTS = transpose(g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == g_strctSystemCodes.m_iSync));
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iSync);

% juice times in plexon time
REW_ON_ts = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iJuiceON);
REW_OFF_ts = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iJuiceOFF);

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

t_plexon = (0:n-1)/adfreq;
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

% Apply calibration to plexon signals
rightEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*rightEyeX_plexon';
rightEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*rightEyeY_plexon';
leftEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*leftEyeX_plexon';
leftEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*leftEyeY_plexon';

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

allFieldNames = fieldnames(trialData{1});

for i = 1:numel(allFieldNames)
    temp = cellfun(@(x) x.(allFieldNames{i}), trialData, 'UniformOutput', false);
    sz = cellfun(@size, temp, 'UniformOutput', false);
    sz(cellfun(@isempty, sz)) = [];
    sz(cellfun(@(x) any(x==0), sz)) = [];
    try
    sz = mode(vertcat(sz{:}), 1);
    temp(cellfun(@isempty, temp)) = {nan(sz)};
    catch
    end
    vars.(allFieldNames{i}) = temp;
end

% Determine stimulus intervals
stimStartTimes = [ones(size([vars.m_fImageFlipON_TS_Kofiko{:}]')),...
    [vars.m_fImageFlipON_TS_Kofiko{:}]']*B;

stimStopTimes = stimStartTimes + [vars.m_fStimulusON_MS{:}]'/1e3;

stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames =  min([vars.numFrames{:}]', [vars.numFrames{:}]' ./ [vars.repframes{:}]');

% Expand relevant variables across frame

% Get trial numbers
trialNumPerFrame = cellfun(@(x, y) repelem(x, y), vars.TrialNum, num2cell(numFrames), 'UniformOutput', false);

% Get trial IDs (insane that there are trial numbers and IDs....)
TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), vars.TrialID, num2cell(numFrames), 'UniformOutput', false);

useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), vars.usebinary, num2cell(numFrames), 'UniformOutput', false);

try
    spatialscale = [vars.spatialscale{:}];
catch
    spatialscale_buffer = g_astrctAllParadigms{1}.DualstimScale.Buffer;
    spatialscale_TS = g_astrctAllParadigms{1}.DualstimScale.TimeStamp;
    spatialscale = nan(size(trialData));
    for i = 1:numel(spatialscale_buffer)
        spatialscale(ImageFlipON_TS_Kofiko >= spatialscale_TS(i)) = spatialscale_buffer(i);
    end

end

BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), vars.BlockID, num2cell(numFrames), 'UniformOutput', false);

% Matrix of unique conditions defined by cloud parameters and block
uniqueCloudConditions =  unique([[vars.usebinary{:}]' spatialscale' [vars.BlockID{:}]'], 'rows');
uniqueCloudConditions(any(isnan(uniqueCloudConditions),2),:) = [];

DualstimPrimaryuseRGBCloudPerFrame = cellfun(@(x, y) repelem(x, y), vars.DualstimPrimaryuseRGBCloud, num2cell(numFrames), 'UniformOutput', false);

StimulusAreaPerFrame =  cellfun(@(x, y) repelem(x, y), vars.m_aiStimulusArea, num2cell(numFrames), 'UniformOutput', false);

stimseq= cellfun(@(x,y) x(1:y:end), vars.stimseq, num2cell(cellfun(@(x) max(0, x), vars.repframes)), 'UniformOutput', false);


%s

X_fixationSpot = cellfun(@(x) x(1), vars.m_pt2iFixationSpot, 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), vars.m_pt2iFixationSpot, 'UniformOutput',false);

UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), vars.UseLeye, num2cell(numFrames), 'UniformOutput', false);
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), vars.UseReye, num2cell(numFrames), 'UniformOutput', false);

toc;

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

trialTypeOfInterest = 'Dual Stim';
trialTypeOfInterestIdx = strcmpi( vars.m_strTrialType, trialTypeOfInterest);

trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), vars.DualstimPrimaryuseRGBCloud);
areaOverZeroIdx = cellfun(@(x) x>0, vars.m_aiStimulusArea);


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

goodFixationIdx = vertcat(vars.m_bMonkeyFixated{:}) | (goodFixationX & goodFixationY);


isTrialOfInterest = trialTypeOfInterestIdx & ...
    goodFixationIdx &...
    trlonset_diffs > 4 &...
    CCidx & ...
    areaOverZeroIdx & ...
    vertcat(vars.DualstimPrimaryuseRGBCloud{:}) == 8 &... % clouds 
    vertcat(vars.UseLeye{:}) == 1 & vertcat(vars.UseReye{:}) == 1;

% %% openiris data
% 
% read in ddpi file
ET = readtable(dpi_fname);

% ddpi timestamps
t_dpi = ET.RightSeconds;
t_dpi = t_dpi - t_dpi(1); % time stamps in s

% ddpi sync signal
sync_dpi = ET.Int0;
sync_dpi = sync_dpi - min(sync_dpi);
sync_dpi = sync_dpi & 1;

% p1
cr_x_right = ET.RightCR1X;
cr_y_right = ET.RightCR1Y;

% p4
p4_x_right = ET.RightCR4X;
p4_y_right = ET.RightCR4Y;

% raw ddpi eye signal
dpi_raw_right = [cr_x_right - p4_x_right, cr_y_right - p4_y_right];

% p1
cr_x_left = ET.LeftCR1X;
cr_y_left = ET.LeftCR1Y;

% p4
p4_x_left = ET.LeftCR4X;
p4_y_left = ET.LeftCR4Y;

% raw ddpi eye signal
dpi_raw_left = [cr_x_left - p4_x_left, cr_y_left - p4_y_left];

%% Align DDPI to plexon

% find rising and falling edges in ddpi computer

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

% delays

delays_dpi = (t_rising_dpi(2:end) - t_falling_dpi(1:end-1))./dt_dpi;
delays_dpi_int = int8(delays_dpi);

t_plexon = 0:(1/adfreq):(length(sync_ad)-1)/adfreq;

% find rising and falling edges of pulse in plexon
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

% now compare openiris and plexon eye signals
% 
% figure, plot(t_dpi_plexon, dpi_raw(:,1)); title('dpi')
% figure, plot(t_plexon, rightEyeX_plexon); title('plexon')

%% Bin dpi eye signal by stimulus intervals
% [~,~, t_dpi_plexon_bin] = histcounts(t_dpi_plexon, stimIntervals);
% 
% t_dpi_plexon_cellArray = accumarray(t_dpi_plexon_bin(:)+1, t_dpi_plexon(:), [], @(x){x});
% 
% % right eye
% rightEyeX_dpi_cellArray = accumarray(t_dpi_plexon_bin(:)+1, dpi_raw_right(:,1), [], @(x){x});
% rightEyeY_dpi_cellArray = accumarray(t_dpi_plexon_bin(:)+1, dpi_raw_right(:,2), [], @(x){x});
% 
% % left eye
% 
% leftEyeX_dpi_cellArray = accumarray(t_dpi_plexon_bin(:)+1, dpi_raw_left(:,1), [], @(x){x});
% leftEyeY_dpi_cellArray = accumarray(t_dpi_plexon_bin(:)+1, dpi_raw_left(:,2), [], @(x){x});

% Apply calibration to dpi signals
% rightEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*rightEyeX_plexon';
% rightEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*rightEyeY_plexon';
% leftEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*leftEyeX_plexon';
% leftEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*leftEyeY_plexon';


%% mean eye position per cloud trial + std

% means
rightEyeX_plexon_calib_meanPerTrial = cellfun(@(x) mean(x, 'omitnan'), rightEyeX_plexon_calib_cellArray(2:2:end));

rightEyeY_plexon_calib_meanPerTrial = cellfun(@(x) mean(x, 'omitnan'), rightEyeY_plexon_calib_cellArray(2:2:end));

leftEyeX_plexon_calib_meanPerTrial = cellfun(@(x) mean(x, 'omitnan'), leftEyeX_plexon_calib_cellArray(2:2:end));

leftEyeY_plexon_calib_meanPerTrial = cellfun(@(x) mean(x, 'omitnan'), leftEyeY_plexon_calib_cellArray(2:2:end));

% stds

rightEyeX_plexon_calib_stdPerTrial = cellfun(@(x) std(x, 'omitnan'), rightEyeX_plexon_calib_cellArray(2:2:end));

rightEyeY_plexon_calib_stdPerTrial = cellfun(@(x) std(x, 'omitnan'), rightEyeY_plexon_calib_cellArray(2:2:end));

leftEyeX_plexon_calib_stdPerTrial = cellfun(@(x) std(x, 'omitnan'), leftEyeX_plexon_calib_cellArray(2:2:end));

leftEyeY_plexon_calib_stdPerTrial = cellfun(@(x) std(x, 'omitnan'), leftEyeY_plexon_calib_cellArray(2:2:end));



meanX =mean([rightEyeX_plexon_calib_meanPerTrial, leftEyeX_plexon_calib_meanPerTrial], 2, 'omitnan');

meanY =mean([rightEyeY_plexon_calib_meanPerTrial, leftEyeY_plexon_calib_meanPerTrial], 2, 'omitnan');

meanPos = mean([meanX, meanY],2, 'omitnan');

%% calibration  fixation points
calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), vars.m_strTrialType) &...
                        goodFixationIdx;
                        %vertcat(vars.m_bMonkeyFixated{:});

isCalibrationTrial(find(isCalibrationTrial) > find(isTrialOfInterest,1)) = 0;
% unique fps during calibration
fps = vertcat(vars.m_pt2iFixationSpot{:});
uniquefps = unique(vertcat(vars.m_pt2iFixationSpot{isCalibrationTrial}), 'rows');



%% plotting

% mean eye position across eyes and axes
figure, plot(meanPos(isTrialOfInterest));

figure; hold on

subplot(2,2,1); hold on
% Left eye

for i = 1:size(uniquefps,1)
    scatter(vertcat(leftEyeX_plexon_calib_cellArray{2*find(isCalibrationTrial & all(fps == uniquefps(i,:),2))}),...
        vertcat(leftEyeY_plexon_calib_cellArray{2*find(isCalibrationTrial & all(fps == uniquefps(i,:),2))}))
    

end
for i = 1:size(uniquefps,1)
    plot(uniquefps(i,1) - ScreenSizeX_pix/2, uniquefps(i,2) - ScreenSizeY_pix/2 ,'ko', 'MarkerFaceColor', 'k')
end
set(gca, 'ydir', 'reverse')

subplot(2,2,2); hold on

% Right Eye

for i = 1:size(uniquefps,1)
    scatter(vertcat(rightEyeX_plexon_calib_cellArray{2*find(isCalibrationTrial & all(fps == uniquefps(i,:),2))}),...
        vertcat(rightEyeY_plexon_calib_cellArray{2*find(isCalibrationTrial & all(fps == uniquefps(i,:),2))}))

     
end

for i = 1:size(uniquefps,1)
    plot(uniquefps(i,1) - ScreenSizeX_pix/2, uniquefps(i,2) - ScreenSizeY_pix/2 ,'ko', 'MarkerFaceColor', 'k')
end
set(gca, 'ydir', 'reverse')



subplot(2,2,3); hold on
% histogram2(leftEyeX_plexon_calib_meanPerTrial(isTrialOfInterest),...
%     leftEyeY_plexon_calib_meanPerTrial(isTrialOfInterest), "FaceColor", "flat", 'XBinEdges', -120:3:120, 'YBinEdges', -120:3:120)

edges = -45:2:45;
centers = edges(1:end-1) + diff(edges)/2;
%N_leftEye = histcounts2(leftEyeY_plexon_calib_meanPerTrial(isTrialOfInterest),leftEyeX_plexon_calib_meanPerTrial(isTrialOfInterest),edges, edges);
N_leftEye = histcounts2(vertcat(leftEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}), edges, edges);

imagesc(centers, centers,N_leftEye'); hold on;
%axis xy
plot(0, 0, 'rx')
set(gca, 'ydir', 'reverse');

colorbar;
axis equal
xlabel ('X (arcmin)')
ylabel('Y (arcmin)')
title('Left eye')
xlim([min(edges), max(edges)])
ylim([min(edges), max(edges)])

subplot(2,2,4); hold on
% histogram2(rightEyeX_plexon_calib_meanPerTrial(isTrialOfInterest),...
%     rightEyeY_plexon_calib_meanPerTrial(isTrialOfInterest), "FaceColor", "flat", 'XBinEdges', -120:3:120, 'YBinEdges', -120:3:120)

%N_rightEye = histcounts2(rightEyeX_plexon_calib_meanPerTrial(isTrialOfInterest),rightEyeY_plexon_calib_meanPerTrial(isTrialOfInterest),  edges, edges);
N_rightEye = histcounts2(vertcat(rightEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}), edges, edges);


imagesc(centers, centers, N_rightEye'); hold on
%axis xy
plot(0,0,'rx')
set(gca, 'ydir', 'reverse');

colorbar;
axis equal
xlabel ('X (arcmin)')
ylabel('Y (arcmin)')
title('Right eye')
xlim([min(edges), max(edges)])
ylim([min(edges), max(edges)])




