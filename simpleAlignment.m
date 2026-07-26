%%

% 3/2026 mjg -- wrote it
% cell arrays:

% For cell arrays of size (2*# trials) x 1, even numbered cells
% correspond to stimulus ON periods, odd numbered cells (starting at 3)
% correspond to interstimulus periods.

% Kofiko_ET_TS_PlexonTime_cellArray [(2*# trials) x 1]: Kofiko eye tracking timestamps. Even numbered
% cells correspond to stimulus on periods, odd numbered cells starting at 3
% correspond to interstimulus periods.

% Kofiko_Xpix_cellArray [(2*# trials) x 1]: Calibrated Kofiko horizontal eye trace (screen coordinates)

% Kofiko_Ypix_cellArray [(2*# trials) x 1]: Calibrated Kofiko vertical eye trace (screen coordinates)

% t_plexon_cellArray [(2*# trials) x 1]: Plexon sample times

% leftEyePupil_plexon_cellArray [(2*# trials) x 1]: Plexon left eye pupil
% trace

% rightEyePupil_plexon_cellArray: same as above for right eye.

% rightEyeX_plexon_calib_cellArray [(2*# trials) x 1]: Plexon right eye
% horizontal trace with Kofiko gains and offsets applied

% rightEyeY_plexon_calib_cellArray [(2*# trials) x 1]: same as above for
% vertical trace.

% leftEyeX_plexon_calib_cellArray [(2*# trials) x 1]: Plexon left eye
% horizontal trace with Kofiko gains and offsets applied

% leftEyeY_plexon_calib_cellArray [(2*# trials) x 1]: same as above for
% vertical trace.

% stimulus_cellArray [# trials x 1]: Each cell is a (stimulus width*stimulus
% height*# chromatic channels) x # frames matrix that contains the cloud
% pixel values (-127 to 127).

% Kofiko_Xpix_frameRate_cellArray [# trials of interest x 1]: Calibrated Kofiko
% horizontal eye trace downsampled to frame rate.

% Kofiko_Xpix_frameRate_cellArray [# trials of interest x 1]: Calibrated Kofiko
% vertical eye trace downsampled to frame rate.

% spk_times_cellArray [# spike sorting batches x 1]: Each cell contains a
% (2*# trials) x 1 cell array. Even cells give spikes that occurred during stimulus ON.
% Odd trials starting at 3 give spikes during interstimulus periods.

% clusterIDForEachSpk_cellArray [# spike sorting batches x 1]: Each cell contains a
% (2*# trials) x 1 cell array. Even cells give IDs corresponding to spikes that occurred during stimulus ON.
% Odd trials starting at 3 give IDs during interstimulus periods.

% ks_batchForEachSpk_cellArray [# spike sorting batches x 1]: Each cell
% contains a (2*# trials) x 1 cell array. Even cells give spike sorting batch corresponding to spikes that occurred during stimulus ON.
% Odd trials starting at 3 give batch during interstimulus periods.

%  stimFrameNumForEachSpk_cellArray  [# spike sorting batches x 1]: Each cell
% contains a (2*# trials) x 1 cell array. Even cells give frames shown during spikes that occurred during stimulus ON.
% Odd trials starting at 3 give frames during interstimulus periods.

% spksPerFrame_cellArray [# trials x 1]: Each cell contains # frames x 1
% vector giving spike counts per frame.


%% Set paths
setPathsAndFlags;
addpath(genpath(fullfile(codedir, 'ConwayExptProcessing')))

%% Initialization for everyone given setup

%addpath(genpath(fullfile(codedir, 'ConwayExptProcessing', 'npy-matlab')))

pl2 = PL2ReadFileIndex(plexon_fname);
fs = pl2.SpikeChannels{1}.SamplesPerSecond;
% Hardcoded values
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45;
minSpikes = 2000;
targ_ETstimtype = 0;
ETdist_thresh=40;
nLags = 10;
%% Load kofiko data
g_strcts = loadKofikoStructs(mainKofiko_fname);
trial = loadKofikoTrialData(kofiko_subfolder,mainKofiko_fname,filenameP);
monkey_name = g_strcts.g_strctAppConfig.m_strctSubject.m_strName;
%% %%%%%%%% Time alignment %%%%%%%%%%
tic;
fprintf('Loading Kofiko eye signals and converting to plexon time\n')
B = KofikoPlexonStrobesLinReg(g_strcts, plexon_fname);

%% Plexon eye data

% make vector of plexon sample times, starting at 0
temp = vertcat(pl2.AnalogChannels{:});
analogChanNames = {temp.Name};
numDigitsInLastAIchan = ceil(log10(sum(contains(analogChanNames, 'AI'))));
chanName =  ['AI' num2str(1, ['%0' num2str(numDigitsInLastAIchan) '.f'])];
[adfreq, n, ~, ~, ~] = plx_ad_v(plexon_fname, chanName);
t_plexon = (0:n-1)/adfreq;

[Kofiko_ET_TS, Kofiko_Xpix, Kofiko_Ypix, KofikoGains, KofikoOffsets, KofikoGains_Plexon, KofikoOffsets_Plexon]...
    = loadKofikoEyeData(g_strcts,t_plexon, B);

Kofiko_ET_TS_PlexonTime = [ones(size(Kofiko_ET_TS)) Kofiko_ET_TS]*B;

chanNums = 1:8;
PlexET_ad_calib = loadPlexonEyeData(plexon_fname, chanNums, rig, plexonAnalogScale, KofikoGains_Plexon);


%% %%%%%%%%%%%%%% Trial analysis %%%%%%%%%%%%%%
tic;

% Determine stimulus intervals
stimStartTimes = [ones(size([trial.m_fImageFlipON_TS_Kofiko]')),...
    [trial.m_fImageFlipON_TS_Kofiko]']*B;

stimStopTimes = stimStartTimes + [trial.m_fStimulusON_MS]'/1e3;

stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames =  min([trial.numFrames], [trial.numFrames] ./ [trial.repframes]);

% Expand relevant variables across frame

% Get trial IDs (insane that there are trial numbers and IDs....)
TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), {trial.TrialID}', num2cell(numFrames)', 'UniformOutput', false);

useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), {trial.usebinary}', num2cell(numFrames)', 'UniformOutput', false);

try
    spatialscale = {trial.spatialscale}';
    
catch
    spatialscale_buffer = g_astrctAllParadigms{1}.DualstimScale.Buffer;
    spatialscale_TS = g_strcts.g_astrctAllParadigms{1}.DualstimScale.TimeStamp;
    spatialscale = nan(size(trial));
    for i = 1:numel(spatialscale_buffer)
        spatialscale(ImageFlipON_TS_Kofiko >= spatialscale_TS(i)) = spatialscale_buffer(i);
    end

end

BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.BlockID}', num2cell(numFrames)', 'UniformOutput', false);

X_fixationSpot = cellfun(@(x) x(1), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);

UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseLeye}', num2cell(numFrames)', 'UniformOutput', false);
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseReye}', num2cell(numFrames)', 'UniformOutput', false);

toc;

%% %%%%%%%%%%%%% Extract fixinfo %%%%%%%%%%%%%
calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), {trial.m_strTrialType});

ETdata.fixloc = vertcat(trial(isCalibrationTrial).m_pt2iFixationSpot);
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
    [nBins + 1, 1], ...   % force size (extra 1 for bin 0)
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
pupil1_plexon_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,3), [], @(x){x});
pupil2_plexon_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,4), [], @(x){x});

% formerly "right eye"
eyeX2_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,5), [], @(x){x});
eyeY2_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,6), [], @(x){x});

% formerly "left eye"
eyeX1_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,7), [], @(x){x});
eyeY1_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, PlexET_ad_calib(:,8), [], @(x){x});

%% Determine trials with good fixation
goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < maxFixationErrorPix)./numel(x)...
    > minFixationDuration, Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < maxFixationErrorPix)./numel(y) ...
    > minFixationDuration, Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

goodFixationIdx = vertcat(trial.m_bMonkeyFixated) | (goodFixationX & goodFixationY);

%% Bin stimulus sequences by trial
tic;
fprintf('Loading stimuli\n')

[stim1_cellArray, stim2_cellArray, stim3_cellArray] = makeStimMatrix(stimpath, trial);

% find dualstim elemetns
trialTypeOfInterest = 'Dual Stim';
trialTypeOfInterestIdx = strcmpi( {trial.m_strTrialType}, trialTypeOfInterest);

trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), {trial.DualstimPrimaryuseRGBCloud});
areaOverZeroIdx = cellfun(@(x) x>0, {trial.m_aiStimulusArea});

isTrialOfInterest = trialTypeOfInterestIdx & ...
    goodFixationIdx &...
    trlonset_diffs > 4 &...
    CCidx & ...
    areaOverZeroIdx & ...
    vertcat(trial.DualstimPrimaryuseRGBCloud) == 8 ; % clouds

stimulus_matrix = horzcat(stimulus_cellArray{isTrialOfInterest});
stimulusET_matrix = horzcat(stimulusET_cellArray{isTrialOfInterest});

stimFrameBins = cellfun(@(t_start, t_stop, nFrames) linspace(t_start,t_stop,nFrames+1),...
    num2cell(stimStartTimes), num2cell(stimStopTimes),...
    num2cell(numFrames), 'UniformOutput', false);

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
fprintf('Loading spike data and computing Robs\n');
% drop test check
spk_offset = droptestcheck(plexon_fname);

% end drop test check %%%%
%% cluster label strings

goodStr = 'good ';
muaStr = 'mua  ';
noiseStr = 'noise';
blankStr = '     ';

%% %%%%%%%%%%%%% Load and organize spike data %%%%%%%%%%%%%
% if using kilosort
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
%unique_array_labels = unique(array_labels);
chan_offset = 0;
cluster_offset = 0;
%chan_offsets(1) = chan_offset;

nSU = [];
nMU = [];
RobsSU =[];
RobsMU =[];

Robs =  cell(1,num_ks_batch);
for ks_batch = 1:num_ks_batch
    this_array_label = array_labels{ks_batch};
    disp(ks_batch)

    % Get label of kilosort batch (often corresponding to array name and range of channels processed)

    % Read in kilosort outputs
    spk_times = readNPY(fullfile(spike_times_folders{ks_batch}, 'spike_times.npy'));
    spk_times = double(spk_times)./fs + spk_offset; % convert to seconds
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

    % Get rid of spike times and clusterIDs that correspond to bad units
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

    stimFrameNumForEachSpk_cellArray{ks_batch} = cell(size(trial,1),1);

    [~,~,stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx)] =...
        cellfun(@(x,y) histcounts(x, 'BinEdges', y), ...
        spk_times_stimON_cellArray(valid_stimFrameBinIdx),...
        stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
    %spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

    Robs{ks_batch} = single(zeros(numel(allUnit_clusterIDs{ks_batch}),  sum(numFrames(isTrialOfInterest))));

    for unit = 1:numel(allUnit_clusterIDs{ks_batch})
        unitID = allUnit_clusterIDs{ks_batch}(unit);

        % each cell is a trial, and contains frame indicies where spikes
        % occurred (repeated frames == multiple spikes on that frame)

        stimFrameNumForEachSpk_thisUnit_cellArray = cell(size(trial,1),1);
        stimFrameNumForEachSpk_thisUnit_cellArray(valid_stimFrameBinIdx) = ...
            cellfun(@(x,y) x(y==unitID), ...
            stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx), ...
            clusterIDForEachSpk_stimON_cellArray(valid_stimFrameBinIdx),...
            'UniformOutput', false);

        spksPerFrame_cellArray =  cell(size(trial,1),1);

        spksPerFrame_cellArray(isTrialOfInterest) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
            stimFrameNumForEachSpk_thisUnit_cellArray(isTrialOfInterest), num2cell(numFrames(isTrialOfInterest)) ,'UniformOutput',false);

        Robs{ks_batch}(unit,:) = [spksPerFrame_cellArray{:}];
    end

    nSU{ks_batch} = length(SU_clusterIDs{ks_batch});
    nMU{ks_batch} = length(MU_clusterIDs{ks_batch});
    RobsSU{ks_batch} = Robs{ks_batch}(1:nSU{ks_batch},:);
    RobsMU{ks_batch} = Robs{ks_batch}(nSU{ks_batch}+1:nSU{ks_batch}+nMU{ks_batch},:);

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


%% Process LFPs
if ~skipLFP
    disp('Processing LFPs')
    tic;
    if (exptdate/10000) < 24  % year less than 2024
        LFPchans{1} = 1:24; % laminar
        LFPchans{2} = [33,40,46,47,52,53,54,59,65,67,71,81,83,89,90,95,98,102,103,109,112,131,138,139,145,146,152,158]; % Nform channels that worked
        LFPchans{3} = 161:256; % utah
    else
        LFPchans{1} = length(pl2.SpikeChannels); % one array so far
    end

    [LFP_adfreq, LFP_n, LFP_ts, ~, ~] = plx_ad_v(plexon_fname, ['FP' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
    for ii=1:length(LFPchans)
        nchans = length(LFPchans{ii});
        %LFP_ad = zeros(nchans, LFP_n);
        LFPs{ii} = zeros(nchans, LFP_n);
        for ch = 1:nchans
            %[~,~,~,~, LFP_ad(i,:)] = plx_ad_v(plexon_fname, ['FP' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
            [~,~,~,~, LFPs{ii}(ch,:)] = plx_ad_v(plexon_fname, ['FP' num2str(LFPchans{ii}(ch), ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
        end
    end
    LFP_times=(0:LFP_n-1)/LFP_adfreq;
    toc;
end

%% %%%%%%%%%%%% Format output like PackageCloud %%%%%%%%%%%%
%allRobs = vertcat(Robs{:});
allRobsSU = vertcat(RobsSU{:});
allRobsMU = vertcat(RobsMU{:});
allRobs = vertcat(allRobsSU, allRobsMU);

%ETgains
ETgains = [Kofiko_GainX(end), Kofiko_GainY(end)];

%ETstim_location
ETstim_location =  [trial{find(isTrialOfInterest,1, 'last')}.secondarystim_bar_rect;...
    trial{find(isTrialOfInterest,1, 'last')}.tertiarystim_bar_rect];

% ETtrace
ETtrace = [[Kofiko_Xpix_frameRate_cellArray{:}]; [Kofiko_Ypix_frameRate_cellArray{:}]];

%ETtrace_raw
ETtrace_raw = transpose([vertcat(eyeX2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(eyeY2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(eyeX1_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(eyeY2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(pupil2_plexon_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(pupil1_plexon_cellArray{2*find(isTrialOfInterest)})]);

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
cloud_scale = transpose(repelem(spatialscale(isTrialOfInterest), numFrames(isTrialOfInterest)));
% spatialScalePerFrame(isTrialOfInterest);

%datafilts
datafilts = ones(size(RobsSU));

%datafiltsMU
datafiltsMU = ones(size(RobsMU));

%dt
dt = 1/60.0; %0.0160; % why not 0.0167?

%electrode_info
electrode_info =[];

%exptdate and exptname
exptdate = filenameP(1:6);
exptname = filenameP;

%fix_location
fix_location = unique(vertcat(trial(isTrialOfInterest).m_pt2iFixationSpot), 'rows');

%fix_size
fix_size = trial{find(isTrialOfInterest, 1, 'last')}.m_fFixationSizePix -1;

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
spike_ts_raw = cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), spk_times_cellArray, 'UniformOutput', false);
spike_ts_raw = transpose(vertcat(spike_ts_raw{:}));
spike_ts_raw = vertcat(spike_ts_raw);

% this will get you spike times relative to trial start, ie in range 0 to 4:
trlsecs = unique([trial(isTrialOfInterest).m_fStimulusON_MS])/1e3;
trialStart = stimStartTimes(isTrialOfInterest);
trialStop  = stimStopTimes(isTrialOfInterest);
spike_ts_raw = spike_ts_raw(:);
trialStart   = trialStart(:);
trialStop    = trialStop(:);

trialIdx = discretize(spike_ts_raw, [trialStart; inf]);

valid = ~isnan(trialIdx) & trialIdx >= 1 & trialIdx <= numel(trialStart);

valid2 = false(size(valid));
valid2(valid) = spike_ts_raw(valid) >= trialStart(trialIdx(valid)) & ...
    spike_ts_raw(valid) <= trialStop(trialIdx(valid));

valid = valid & valid2;

spike_ts_raw = spike_ts_raw(valid);
trialIdx     = trialIdx(valid);
spike_ts     = spike_ts_raw - trialStart(trialIdx) + (trialIdx - 1) * trlsecs;

%stim
stim = reshape(stimulus_matrix, 60,60,3,[]);

%stimET
if stimET
    stimET = reshape(stimulusET_matrix, 60,60,3,[]);
end

%stim_area
stim_area = [trial(isTrialOfInterest).m_aiStimulusArea]';

%stim_location
modal_stim_area = mode([trial(isTrialOfInterest).m_aiStimulusArea]');
if computerLocation < 10
    % this line does not work in Dan lab
    medianTiledStimulusRect = median(cat(3, trial.m_aiTiledStimulusRect),3, 'omitmissing');
else
    % previous code that worked
    if isfield(trial{1}, 'm_aiTiledStimulusRect')
        TiledStimulusRect = cellfun(@(x) x.m_aiTiledStimulusRect, trial, 'UniformOutput', false);
    else
        StimulusRect = cellfun(@(x) x.m_aiStimulusRect, trial, 'UniformOutput', false);
        TiledStimulusRect = StimulusRect; %%%%% Only for old datasets before tiling (e.g., 7/07/22)
    end
    medianTiledStimulusRect = median(cat(3, TiledStimulusRect{:}),3);
end

stim_location = [medianTiledStimulusRect(:,1:2) medianTiledStimulusRect(:,1:2)+60];

%stim_location_deltas: first two elements of first row of TiledStimulusRect
%minus median stim_location x and y

if (computerLocation < 10) || isfield(trial{1}, 'm_aiTiledStimulusRect')
    stim_location_deltas = cellfun(@(x) x(1,1:2) - stim_location(1,1:2), {trial(isTrialOfInterest).m_aiTiledStimulusRect}, 'UniformOutput',false);
else
    stim_location_deltas = cellfun(@(x) x(1,1:2) - stim_location(1,1:2), TiledStimulusRect(isTrialOfInterest), 'UniformOutput',false);
end
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

%% remap cluster ids from 1 to number of clusters

uniqueClusterIDs = vertcat(vertcat(SU_clusterIDs{:}), vertcat(MU_clusterIDs{:}));

spikeIDs_unsorted = accumarray(transpose(1:numel(clusterIDs)), clusterIDs', [], @(x) find(uniqueClusterIDs==x));
[spikeIDs, sortBySpikeID] = sort(spikeIDs_unsorted);
spike_ts_sorted = spike_ts(sortBySpikeID);


%% Make list of reward times on each trial

% juice times in plexon time
REW_ON_ts = events.timeStamps(events.strobeNumber == g_strcts.g_strctSystemCodes.m_iJuiceON);
REW_OFF_ts = events.timeStamps(events.strobeNumber == g_strcts.g_strctSystemCodes.m_iJuiceOFF);

reward_on_ts = [];
reward_off_ts = [];

if computerLocation < -100
    [~, ~, rew_off_bin] = histcounts(REW_OFF_ts, stimIntervals);
    REW_OFF_cellArray = accumarray(rew_off_bin(:)+1,...
        REW_OFF_ts(:),...
        [nBins + 1,1],...
        @(x){x},...
        {[]});

    [~, ~, rew_on_bin] = histcounts(REW_ON_ts, stimIntervals);
    REW_ON_cellArray = accumarray(rew_on_bin(:)+1,...
        REW_ON_ts(:),...
        [nBins + 1,1],...
        @(x){x},...
        {[]});

    REW_ON_cellArray = cellfun(@(rewardTime, trialStartTime) rewardTime - trialStartTime, REW_ON_cellArray, REW_OFF_cellArray);

    reward_on_ts = vertcat(REW_ON_cellArray{2*find(isTrialOfInterest)});
    reward_off_ts = vertcat(REW_OFF_cellArray{2*find(isTrialOfInterest)});
else
    for tr = 1:length(trial_start_ts)
        t = trial_start_ts(tr);
        % Identify reward times for the current trial
        rewardON = REW_ON_ts(REW_ON_ts >= t & REW_ON_ts < t+4.0);
        rewardOFF = REW_OFF_ts(REW_OFF_ts >= t & REW_OFF_ts < t+4.0);
        for ii=1:length(rewardON)
            reward_on_ts(end+1) = rewardON(ii)-t + 4*(tr-1);
        end
        for ii=1:length(rewardOFF)
            reward_off_ts(end+1) = rewardOFF(ii)-t + 4*(tr-1);
        end
    end
end

%% Add fields to data struct
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
data.reward_on_ts = reward_on_ts;
data.reward_off_ts = reward_off_ts;

data.valid_data = valid_data;

data.spikeSortingBatch = spikeSortingBatch;
data.spikeSortingBatchMU = spikeSortingBatchMU;

data.ks_folders = ks_folders;
data.chan_offsets = chan_offsets;
data.cluster_offsets = cluster_offsets;

%% modify array_labels so its python readable
unique_array_labels = unique(array_labels);

[C, IA, IC] = unique(array_labels);
data.arrayPerSU = IC(data.spikeSortingBatch);
data.arrayPerMU = IC(data.spikeSortingBatchMU);

%array_labels = cellfun(@(num,lab) [num2str(num) lab], num2cell(1:num_ks_batch), array_labels, 'UniformOutput', false);
%data.array_labels = horzcat(array_labels{:});
%data.array_labels = array_labels;
unique_labels = horzcat(cellfun(@(num,lab) ['_' num2str(num) lab], num2cell(1:length(C)), unique_array_labels, 'UniformOutput', false));
data.array_labels = horzcat(unique_labels{:});

data.cluster = cluster;
data.clusterMU = clusterMU;


%% Saving
if saving
    disp('Saving')
    if ~isdir(savepath)
        mkdir(savepath);
    end

    switch targ_ETstimtype
        case 0; curETstimtype = 'NA';
        case 1; curETstimtype='1D';
        case 7; curETstimtype='CC';
    end

    switch_stimtype = unique(vertcat(trial(isTrialOfInterest).DualstimPrimaryuseRGBCloud));
    switch switch_stimtype
        case 0; curstimstype='GT';
        case 3; curstimstype='HL';
        case 6; curstimstype='HC';
        case 8; curstimstype='CC';
    end


    array_label_filepart = [cellfun(@(x) [x '_'], unique_array_labels(1:end-1), 'UniformOutput', false) unique_array_labels(end)];
    array_label_filepart = horzcat(array_label_filepart{:});

    useofflinesorting = 1;  % why is this flag all the way down here?
    %FullExpt_ET_filename = sprintf( '%s_FullExpt_ET.mat', filenameP );
    FullExpt_ET_filename = sprintf( 'K%s_FullExpt_ET.mat', exptdate );

    if computerLocation < 10
        data_filename=[monkey_name '_' exptname(1:6) '_' array_label_filepart '_' curstimstype '_ET' curETstimtype '_v10.mat'];
        fixinfo_filename=[filenameP '_fixinfo.mat'];
    else
        data_filename = ['K' exptdate '_' curstimstype '_ET' curETstimtype '_v10.mat'];
        fixinfo_filename = sprintf( 'K%s_fixinfo.mat', exptdate );
    end
    LFPfilename = sprintf( 'K%s_LFPs.mat', exptdate );


    %data_filename=[ 'K' exptdate '_' curstimstype '_ET' curETstimtype '_v10.mat'];
    %fixinfo_filename=['K' exptdate '_fixinfo.mat'];
    save(fullfile(savepath, data_filename),  '-struct', 'data', '-v7.3'); % save packaged cloud data
    save(fullfile(savepath, fixinfo_filename), '-struct', 'ETdata', '-v7.3') % save fixinfo
    save(fullfile(savepath, FullExpt_ET_filename), 'PlexET_ad_calib', 'PlexET_times', '-v7.3'); % save FullExpt_ET

    if ~skipLFP
        trial_start_inds = floor(trial_start_ts*1000);
        LFP_ad = LFPs{1};
        if length(LFPs) == 1  % then this is the only array
            save(fullfile(savepath, LFPfilename), 'LFP_ad', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
        elseif length(LFPs) == 3  % then old-school (2022 expt date)
            LFPa2 = LFPs{2};
            LFPa3 = LFPs{3};
            %save(fullfile(savepath, [filenameP '_LFP.mat']), 'LFP_ad', 'LFPa2', 'LFPa3', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
            save(fullfile(savepath, LFPfilename), 'LFP_ad', 'LFPa2', 'LFPa3', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
        else
            disp('Have not programmed in this LFP-array config yet')
        end
    end

end


% trial_start_inds = round(trialstart_plx'.*1000);