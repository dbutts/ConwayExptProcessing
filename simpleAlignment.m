%% Set paths

% "dirpath" should contain:
%   - All Kofiko mat files as well as the plexon file (which will be
%   searched for recursively)
%   - Plexon (.pl2) file

%dirpath = '~/Data/V1_Fovea/Sprout/260205';
dirpath = '/mnt/isilon/DATA/monkey_ephys/Jocamo/2022to23_ArrayExpts/';

% "ks_path" should contain:
%   -Kilosort outputs (spike_times.npy, spike_clusters.npy, and
%   cluster_info.tsv)

%ks_path = '/mnt/isilon/PROJECTS/V1_Fovea/processing/250529_152043_Jacomo';
ks_path = '/home/bizon/Data/V1_Fovea/Jocamo/220715/220715_131937_Jacomo/';
% Should point to stimuli mat files
stimpath = '/home/bizon/Processing/Cloudstims_calib_01_2022/';

% Add ConwayExptProcessing Dependencies to path
addpath(genpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/'));
addpath(genpath('/home/bizon/Processing'));

% File prefix for Kofiko and plexon files
%filenameP = '250529_152043_Jacomo';
filenameP = '220715_131937_Jacomo';
plexon_dir = dir(fullfile(dirpath, '**', [filenameP '.pl2']));
plexon_fname = fullfile(plexon_dir(1).folder, plexon_dir(1).name);

%% Hardcoded values
fps = 60; % display frames per second
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
% Get Kofiko eye signal gain timestamps
Kofiko_GainX_TS = g_strctEyeCalib.GainX.TimeStamp;
Kofiko_GainY_TS = g_strctEyeCalib.GainY.TimeStamp;
% Get Kofiko eye signal offsets
Kofiko_CenterX = g_strctEyeCalib.CenterX.Buffer;
Kofiko_CenterY = g_strctEyeCalib.CenterY.Buffer;
% Get Kofiko eye signal offset timestamps
Kofiko_CenterX_TS = g_strctEyeCalib.CenterX.TimeStamp;
Kofiko_CenterY_TS = g_strctEyeCalib.CenterY.TimeStamp;

% Determine the Kofiko gains and offsets for each eye signal time stamp

Kofiko_GainX_forEachTimeStamp = zeros(size(Kofiko_ET_TS));
for i = 1:numel(Kofiko_GainX)
    idx = Kofiko_ET_TS >= Kofiko_GainX_TS(i);
    Kofiko_GainX_forEachTimeStamp(idx) = Kofiko_GainX(i);
end

Kofiko_GainY_forEachTimeStamp = zeros(size(Kofiko_ET_TS));
for i = 1:numel(Kofiko_GainY)
    idx = Kofiko_ET_TS >= Kofiko_GainY_TS(i);
    Kofiko_GainY_forEachTimeStamp(idx) = Kofiko_GainY(i);
end

Kofiko_CenterX_forEachTimeStamp = zeros(size(Kofiko_ET_TS));
for i = 1:numel(Kofiko_CenterX)
    idx = Kofiko_ET_TS >= Kofiko_CenterX_TS(i);
    Kofiko_CenterX_forEachTimeStamp(idx) = Kofiko_CenterX(i);
end

Kofiko_CenterY_forEachTimeStamp = zeros(size(Kofiko_ET_TS));
for i = 1:numel(Kofiko_CenterY)
    idx = Kofiko_ET_TS >= Kofiko_CenterY_TS(i);
    Kofiko_CenterY_forEachTimeStamp(idx) = Kofiko_CenterY(i);
end

% Convert raw eye signals to screen coordinates (pixels) by:
% [RawSignal-SignalOffset]*Gain + ScreenCenter = FixationPoint
% SignalOffset is CenterX or CenterY

Kofiko_Xpix = Kofiko_GainX_forEachTimeStamp.*(Kofiko_EyeRawX - Kofiko_CenterX_forEachTimeStamp) + ScreenSizeX_pix/2;
Kofiko_Ypix = Kofiko_GainY_forEachTimeStamp.*(Kofiko_EyeRawY - Kofiko_CenterY_forEachTimeStamp) + ScreenSizeY_pix/2;

toc;
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


% Get stimulus area on each trial
StimulusArea = cellfun(@(x) x.m_aiStimulusArea, trialData, 'UniformOutput', false);
StimulusArea(cellfun(@isempty, StimulusArea)) = {nan};

% Get stimulus sequence on each trial
stimseq = cellfun(@(x) x.stimseq, trialData, 'UniformOutput',false);

% Get actual sequence accounting for repFrames
stimseq= cellfun(@(x,y) x(1:y:end), stimseq, repFrames, 'UniformOutput', false);

% Get fixation locations
pt2iFixationSpot = cellfun(@(x) x.m_pt2iFixationSpot, trialData, 'UniformOutput',false);
pt2iFixationSpot(cellfun(@isempty,pt2iFixationSpot)) = {[0 0]}; % just in case
X_fixationSpot = cellfun(@(x) x(1), pt2iFixationSpot, 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), pt2iFixationSpot, 'UniformOutput',false);

% Get Kofiko determination of monkey fixation
MonkeyFixated = cellfun(@(x) logical(x.m_bMonkeyFixated), trialData);

% To get kofiko eye date for trial trialNum(t):
% Kofiko_ET_TS_PlexonTime(Kofiko_ET_TS_PlexonTime_Bin == stimON_Kofiko_ET_TS_PlexonTime_Bin(t))

toc;

%% Bin eye signal timestamps by stimulus intervals (i.e., by trial)

[~,~,Kofiko_ET_TS_PlexonTime_Bin] = histcounts(Kofiko_ET_TS_PlexonTime, stimIntervals); % wont really need this ever, I think
uniqueKofiko_ET_TS_PlexonTime_Bin= unique(Kofiko_ET_TS_PlexonTime_Bin);
%stimON_Kofiko_ET_TS_PlexonTime_Bin = uniqueKofiko_ET_TS_PlexonTime_Bin(1:2:end);

% I use accumarray to create cell array where even numbered cells give eye
% signal timestamps for stimulus on periods (stimulus onset to offset),
% and odd numbered (starting at 3) cells give eye signal timestamps between stimulus on
% periods (stimulus offset to subsequent onset).The first cell corresponds
% to invalid timestamps (i.e. outside of stimIntervals)
Kofiko_ET_TS_PlexonTime_cellArray =accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_ET_TS_PlexonTime(:), [], @(x){x});

% Do as above for eye signal X and Y coordinates:
Kofiko_Xpix_cellArray = accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_Xpix(:), [], @(x){x});
Kofiko_Ypix_cellArray = accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_Ypix(:), [], @(x){x});

%% Determine trials with good fixation
fix_thresh = 0.6; 

goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < 45)./numel(x)...
    > fix_thresh, Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < 45)./numel(y) ...
    > fix_thresh, Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

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
    LumScale = 0.1085;
    if thisUsebinary == 0
        load(fullfile(stimpath, sprintf('Cloudstims_Chrom_size60_scale%d_%02d.mat', thisSpatialscale, thisBlockID)));
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
    elseif thisUsebinary == 1
        load(fullfile(stimpath,sprintf('Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat', thisSpatialscale, thisBlockID)))
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
    elseif thisUsebinary == 2
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

dualStimIdx = strcmpi(trialType, 'Dual Stim');
trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), DualstimPrimaryuseRGBCloud);
areaOverZeroIdx = cellfun(@(x) x>0, StimulusArea);

trialsOfInterestIdx = dualStimIdx & ...
                      goodFixationIdx &...
                      trlonset_diffs > 4 &...
                      CCidx & ...
                      areaOverZeroIdx;

stimulus_matrix = horzcat(stimulus_cellArray{trialsOfInterestIdx});

stimFrameBins = cellfun(@(x, y, z) linspace(x,y,z+1), num2cell(stimStartTimes), num2cell(stimStopTimes), num2cell(numFrames), 'UniformOutput', false);
valid_stimFrameBinIdx = cellfun(@(x) numel(x)>=2, stimFrameBins);
toc;
%% load in spike data
tic;
fprintf('Loaidng spike data and computing Robs\n');
% drop test check
pl2 = PL2ReadFileIndex(plexon_fname);
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
    disp('Experiment kilosorted and ready for curation!')
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
spike_times_seconds_dir = dir(fullfile(ks_path, '**/spike_times_seconds.npy'));
spike_clusters_dir = dir(fullfile(ks_path, '**/spike_clusters.npy'));
cluster_info_dir = dir(fullfile(ks_path, '**/cluster_info.tsv'));

spike_times_seconds_folders = {spike_times_seconds_dir(:).folder};
spike_clusters_folders = {spike_clusters_dir(:).folder};
cluster_info_folders = {cluster_info_dir(:).folder};

% For each folder with kilosort outputs
for ks_batch = 1:length(spike_times_seconds_dir)

    try
         % Get label of kilosort batch (often corresponding to array name and range of channels processed)
        [~,ks_folders{ks_batch},~] = fileparts(spike_times_seconds_folders{ks_batch});

        % Read in kilosort outputs
        spk_times = readNPY(fullfile(spike_times_seconds_folders{ks_batch}, 'spike_times_seconds.npy')) + opts.spk_offset;
        spk_clusters = readNPY(fullfile(spike_times_seconds_folders{ks_batch}, 'spike_clusters.npy'));
        spk_info = tdfread(fullfile(spike_times_seconds_folders{ks_batch}, 'cluster_info.tsv'));

        % Find indices of units labeled "good", "mua", or ""
        isGood = cellfun(@(x) strcmpi(deblank(x), 'good'), cellstr(spk_info.group));
        isMua = cellfun(@(x) strcmpi(deblank(x), 'mua'), cellstr(spk_info.group));
        isBlank = cellfun(@(x) isempty(deblank(x)), cellstr(spk_info.group));
        minSpikes = 2000;
        hasMinSpikes = spk_info.n_spikes > minSpikes;

        SU_IDs{ks_batch} = spk_info.cluster_id(isGood & hasMinSpikes);
        MU_IDs{ks_batch} = spk_info.cluster_id((isMua | isBlank) & hasMinSpikes);
        allUnit_IDs{ks_batch} =[SU_IDs{ks_batch}; MU_IDs{ks_batch}];

        % bin spikes by trial
        [~,~,spk_times_bin{ks_batch}] = histcounts(spk_times, stimIntervals);

        % each cell gives spike times for stim ON and stim OFF periods
        spk_times_cellArray{ks_batch} =accumarray(spk_times_bin{ks_batch}(:)+1, spk_times(:), [], @(x){x});
        clusterIDForEachSpk_cellArray{ks_batch} = accumarray(spk_times_bin{ks_batch}(:)+1, spk_clusters(:), [], @(x){x});

        clusterIDForEachSpk_stimON_cellArray = clusterIDForEachSpk_cellArray{ks_batch}(2:2:end); 
        spk_times_stimON_cellArray =  spk_times_cellArray{ks_batch}(2:2:end);

        spksPerFrame_cellArray =  cell(size(trialData,1),1);

        stimFrameNumForEachSpk_cellArray = cell(size(trialData,1),1);

        [~,~,stimFrameNumForEachSpk_cellArray(valid_stimFrameBinIdx)] =...
            cellfun(@(x,y) histcounts(x, 'BinEdges', y), ...
            spk_times_stimON_cellArray(valid_stimFrameBinIdx),...
            stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
        %spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

        Robs{ks_batch} = single(zeros(numel(allUnit_IDs{ks_batch}),  sum(numFrames(trialsOfInterestIdx))));

        for unit = 1: numel(allUnit_IDs{ks_batch})
            unitID = allUnit_IDs{ks_batch}(unit);

            % each cell is a trial, and contains frame indicies where spikes
            % occurred (repeated frames == multiple spikes on that frame)


            stimFrameNumForEachSpk_thisUnit_cellArray = cell(size(trialData,1),1);
            stimFrameNumForEachSpk_thisUnit_cellArray(valid_stimFrameBinIdx) = ...
                cellfun(@(x,y) x(y==unitID), ...
                stimFrameNumForEachSpk_cellArray(valid_stimFrameBinIdx), ...
                clusterIDForEachSpk_stimON_cellArray(valid_stimFrameBinIdx),...
                'UniformOutput', false);

  
            spksPerFrame_cellArray(trialsOfInterestIdx) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
                stimFrameNumForEachSpk_thisUnit_cellArray(trialsOfInterestIdx), num2cell(numFrames(trialsOfInterestIdx)) ,'UniformOutput',false);

            Robs{ks_batch}(unit,:) = [spksPerFrame_cellArray{:}];

        end
    catch ME
        warning("Could not process %s", ks_folders{ks_batch})
        disp(ME.message)
    end

end

toc;

%% Compute STAs
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
%% Plotting 

for ks_batch = 1:length(STA)
    for unit = 1:size(STA{ks_batch},1)
        figure;

        subplot(1,4,1); title('Lum lag 2')
        imagesc(circshift(squeeze(STA{ks_batch}(unit,:,:,1,3)), [30 30])), colormap gray

        subplot(1,4,2); title('Lum lag 3')
        imagesc(circshift(squeeze(STA{ks_batch}(unit,:,:,1,4)), [30 30])), colormap gray

        subplot(1,4,3); title('Lum lag 4')
        imagesc(circshift(squeeze(STA{ks_batch}(unit,:,:,1,5)), [30 30])), colormap gray

        subplot(1,4,4); title('Lum lag 5')
        imagesc(circshift(squeeze(STA{ks_batch}(unit,:,:,1,6)), [30 30])), colormap gray


        sgtitle(sprintf('%s, cluster ID: %i', replace(ks_folders{ks_batch}, '_', ' '), allUnit_IDs{ks_batch}(unit)))
    end
end
