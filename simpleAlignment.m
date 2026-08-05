
%function data = simpleAlignment
%%
% 3/2026 mjg -- wrote it
% 7/2026 mjg -- major refactoring

%% Set paths
setPathsAndFlags;
addpath(genpath(fullfile(codedir, 'ConwayExptProcessing')))

%% Initialization for everyone given setup

%addpath(genpath(fullfile(codedir, 'ConwayExptProcessing', 'npy-matlab')))

pl2 = PL2ReadFileIndex(plexon_fname);
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

[Kofiko_ET_TS, Kofiko_Xpix,...
    Kofiko_Ypix, KofikoGains,...
    KofikoOffsets, KofikoGains_Plexon,...
    KofikoOffsets_Plexon]...
    = loadKofikoEyeData(g_strcts,t_plexon, B);

Kofiko_ET_TS_PlexonTime = [ones(size(Kofiko_ET_TS)) Kofiko_ET_TS]*B;

chanNums = 1:8;
PlexET_ad_calib = loadPlexonEyeData(plexon_fname, chanNums, rig, plexonAnalogScale, KofikoGains_Plexon);


%% %%%%%%%%%%%%%% Trial analysis %%%%%%%%%%%%%%
tic;

% Determine stimulus intervals
stimTiming = getStimTiming(trial, B);

% Expand relevant variables across frame

% Get trial IDs (insane that there are trial numbers and IDs....)
TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), {trial.TrialID}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), {trial.usebinary}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

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

BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.BlockID}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
StimulusAreaPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.m_aiStimulusArea}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
DualstimPrimaryuseRGBCloudPerFrame = cellfun(@(x, y) repelem(x, y),{trial.DualstimPrimaryuseRGBCloud}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

DiscDiameterPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.m_iDiscDiameter}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

X_fixationSpot = cellfun(@(x) x(1), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);

UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseLeye}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseReye}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

toc;
%% %%%%%%%%%%%%% Load and organize spike data %%%%%%%%%%%%%
% if using kilosort
onlineSortingOnly = false;
spkData = organizeSpikeDataByTrial(stimTiming.stimIntervals,plexon_fname, minSpikes, ks_path, onlineSortingOnly);

%% %%%%%%%%%%%%% Extract fixinfo %%%%%%%%%%%%%
calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), {trial.m_strTrialType});

ETdata.fixloc = vertcat(trial(isCalibrationTrial).m_pt2iFixationSpot);
ETdata.plxonset = stimTiming.stimStartTimes(isCalibrationTrial);

%% Bin Kofiko eye signal timestamps by stimulus intervals (i.e., by trial)
[~,~,Kofiko_ET_TS_PlexonTime_Bin] = histcounts(Kofiko_ET_TS_PlexonTime, stimTiming.stimIntervals);
%stimON_Kofiko_ET_TS_PlexonTime_Bin = uniqueKofiko_ET_TS_PlexonTime_Bin(1:2:end);

% I use accumarray to create cell array where even numbered cells give eye
% signal timestamps for stimulus on periods (stimulus onset to offset),
% and odd numbered (starting at 3) cells give eye signal timestamps between stimulus on
% periods (stimulus offset to subsequent onset).The first cell corresponds
% to invalid timestamps (i.e. outside of stimIntervals)

nBins = numel(stimTiming.stimIntervals) - 1;

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
[~,~, t_plexon_bin] = histcounts(t_plexon, stimTiming.stimIntervals);

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

%% Load stim matrices
tic;
fprintf('Loading stimuli\n')
[stim1_cellArray, stim2_cellArray, stim3_cellArray] = makeStimMatrix(stimpath, trial, LumScale);

discProbeTrialIdx =  strcmpi( {trial.m_strTrialType}', 'Disc Probe');
ccloudTrialIdx = strcmpi( {trial.m_strTrialType}', 'Dual Stim') & vertcat(trial.DualstimPrimaryuseRGBCloud) == 8;
hartleyTrialIdx = strcmpi( {trial.m_strTrialType}', 'Dual Stim') & vertcat(trial.DualstimPrimaryuseRGBCloud) == 6;

ccloudTrialIdx = ccloudTrialIdx & goodFixationIdx;
discProbeTrialIdx = discProbeTrialIdx & goodFixationIdx;
hartleyTrialIdx = hartleyTrialIdx & goodFixationIdx;

isTrialOfInterestIndices = {ccloudTrialIdx, discProbeTrialIdx, hartleyTrialIdx};
trialLabels = {'ccloud', 'discProbe', 'hartley'};

% stim1_matrix = horzcat(stim1_cellArray{isTrialOfInterest});
% stim2_matrix = horzcat(stim2_cellArray{isTrialOfInterest});

%% Build Robs
%Robs_strct = struct;
for i = 1:numel(isTrialOfInterestIndices)

    isTrialOfInterest = isTrialOfInterestIndices{i};
    stim1_matrix = horzcat(stim1_cellArray{isTrialOfInterest});
    stim2_matrix = horzcat(stim2_cellArray{isTrialOfInterest});
    stim3_matrix = horzcat(stim3_cellArray{isTrialOfInterest});

    %% resample eye signal at frame rate
    % NOTE: Without 'extrap' argument, there will be NANs:

    Kofiko_Xpix_frameRate_cellArray = ...
        cellfun(@(t,x, start, stop, fixspot, n) interp1(t, x - fixspot, linspace(start, stop, n), 'linear'),...
        Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
        Kofiko_Xpix_cellArray(2*find(isTrialOfInterest)),...
        num2cell(stimTiming.stimStartTimes(isTrialOfInterest)),...
        num2cell(stimTiming.stimStopTimes(isTrialOfInterest)),...
        X_fixationSpot(isTrialOfInterest),...
        num2cell(stimTiming.numFrames(isTrialOfInterest)), ...
        'UniformOutput', false);

    Kofiko_Ypix_frameRate_cellArray = ...
        cellfun(@(t,y, start, stop, fixspot, n) interp1(t, y - fixspot, linspace(start, stop, n), 'linear'),...
        Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
        Kofiko_Ypix_cellArray(2*find(isTrialOfInterest)),...
        num2cell(stimTiming.stimStartTimes(isTrialOfInterest)),...
        num2cell(stimTiming.stimStopTimes(isTrialOfInterest)),...
        Y_fixationSpot(isTrialOfInterest),...
        num2cell(stimTiming.numFrames(isTrialOfInterest)), ...
        'UniformOutput', false);

    % Robs
    if isfield(spkData, 'spkDataOffline')
        Robs_strct(i) = buildRobs(spkData.spkDataOffline, stimTiming, isTrialOfInterest);
        Robs_strct_online(i) = buildRobs(spkData.spkDataOnline, stimTiming, isTrialOfInterest);
    else
        Robs_strct(i) = buildRobs(spkData.spkDataOnline, stimTiming, isTrialOfInterest);
    end


    %% %%%%%%%%%%%% Format output like PackageCloud %%%%%%%%%%%%
    %allRobs = vertcat(Robs{:});
    % allRobsSU = vertcat(RobsSU{:});
    % allRobsMU = vertcat(RobsMU{:});
    % allRobs = vertcat(allRobsSU, allRobsMU);

    %ETgains
    ETgains = [g_strcts.g_strctEyeCalib.GainX.Buffer(end), g_strcts.g_strctEyeCalib.GainY.Buffer(end)];

    %ETstim_location
    ETstim_location =  [trial(find(isTrialOfInterest,1, 'last')).secondarystim_bar_rect;...
        trial(find(isTrialOfInterest,1, 'last')).tertiarystim_bar_rect];

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

    nSU = numel(Robs_strct(i).SU_clusters);
    nMU = numel(Robs_strct(i).MU_clusters);

    RobsSU = Robs_strct(i).Robs(1:nSU,:);

    %RobsMU
    RobsMU = Robs_strct(i).Robs(nSU+1:end,:);

    %RobsMU_probe_ID
    RobsMU_probe_ID = Robs_strct(i).MU_chans;

    %RobsMU_rating
    RobsMU_rating = [];
    %blockID

    % Robs_probe_ID
    Robs_probe_ID = Robs_strct(i).SU_chans;

    %Robs rating
    Robs_rating =[];
    blockID = [BlockIDPerFrame{isTrialOfInterest}];

    %block_inds
    block_inds = [cumsum(stimTiming.numFrames(isTrialOfInterest))' - (stimTiming.numFrames(find(isTrialOfInterest,1,'first'))-1); cumsum(stimTiming.numFrames(isTrialOfInterest))'];

    %cloud_area
    cloud_area = [StimulusAreaPerFrame{isTrialOfInterest}];

    disc_diameter = [DiscDiameterPerFrame{isTrialOfInterest}];

    %cloud_binary
    cloud_binary = [useBinaryPerFrame{isTrialOfInterest}];

    %cloud_scale
    cloud_scale = transpose(repelem(spatialscale(isTrialOfInterest), stimTiming.numFrames(isTrialOfInterest)));
    % spatialScalePerFrame(isTrialOfInterest);

    %datafilts
    datafilts = ones(size(RobsSU));

    %datafiltsMU
    datafiltsMU = ones(size(RobsMU));

    %dt
    dt = 1/60.0;

    %electrode_info
    electrode_info =[];

    %exptdate and exptname
    exptdate = filenameP(1:6);
    exptname = filenameP;

    %fix_location
    fix_location = unique(vertcat(trial(isTrialOfInterest).m_pt2iFixationSpot), 'rows');

    %fix_size
    fix_size = [trial(find(isTrialOfInterest, 1, 'last')).m_fFixationSizePix] -1;

    %pixel_size
    pixel_size = 1;

    %sacc_inds
    sacc_inds=[];

    %ss
    clusterIDs= cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), {spkData.spkDataOnline.spk_clusters_cellArray}, 'UniformOutput', false);
    clusterIDs = transpose(vertcat(clusterIDs{:}));
    clusterIDs = vertcat(clusterIDs);

    assert(numel(unique(clusterIDs)) == nSU + nMU);

    %spike_ts
    spike_ts_raw = cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), {spkData.spkDataOnline.spk_times_cellArray}, 'UniformOutput', false);
    spike_ts_raw = transpose(vertcat(spike_ts_raw{:}));
    spike_ts_raw = vertcat(spike_ts_raw);

    % this will get you spike times relative to trial start, ie in range 0 to 4:
    trlsecs = unique([trial(isTrialOfInterest).m_fStimulusON_MS])/1e3;
    trialStart = stimTiming.stimStartTimes(isTrialOfInterest);
    trialStop  = stimTiming.stimStopTimes(isTrialOfInterest);
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

    stimW = 60;
    stimH =60;

    stim = stim1_matrix;
    if strcmpi(trialLabels{i}, 'ccloud')
        stim = reshape(stim, stimH,stimW,3,[]);
    end

    stim2 = stim2_matrix;
    stim3 = stim3_matrix;

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
    trial_start_ts = [stimTiming.stimStartTimes(isTrialOfInterest)]';

    %useLeye
    useLeye = [UseLeyePerFrame{isTrialOfInterest}];

    %useReye
    useReye = [UseReyePerFrame{isTrialOfInterest}];

    %valid_data
    totalFrames = sum(stimTiming.numFrames(isTrialOfInterest));
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

    %% remap cluster ids from 1 to number of clusters

    uniqueClusterIDs = [Robs_strct(i).SU_clusters; Robs_strct(i).MU_clusters];

    spikeIDs_unsorted = accumarray(transpose(1:numel(clusterIDs)), clusterIDs', [], @(x) find(uniqueClusterIDs==x));
    [spikeIDs, sortBySpikeID] = sort(spikeIDs_unsorted);
    spike_ts_sorted = spike_ts(sortBySpikeID);


    %% Make list of reward times on each trial

    % juice times in plexon time
    [events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);
    REW_ON_ts = events.timeStamps(events.strobeNumber == g_strcts.g_strctSystemCodes.m_iJuiceON);
    REW_OFF_ts = events.timeStamps(events.strobeNumber == g_strcts.g_strctSystemCodes.m_iJuiceOFF);

    reward_on_ts = [];
    reward_off_ts = [];

    if computerLocation < -100
        [~, ~, rew_off_bin] = histcounts(REW_OFF_ts, stimTiming.stimIntervals);
        REW_OFF_cellArray = accumarray(rew_off_bin(:)+1,...
            REW_OFF_ts(:),...
            [nBins + 1,1],...
            @(x){x},...
            {[]});

        [~, ~, rew_on_bin] = histcounts(REW_ON_ts, stimTiming.stimIntervals);
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
    data(i).ETgains = ETgains;
    data(i).ETstim_location = ETstim_location;
    data(i).ETtrace = ETtrace;
    data(i).ETtrace_raw = ETtrace_raw;
    data(i).Robs = RobsSU;
    data(i).RobsMU = RobsMU;
    data(i).RobsMU_probe_ID = RobsMU_probe_ID;
    data(i).RobsMU_rating = RobsMU_rating;
    data(i).Robs_probe_ID = Robs_probe_ID;
    data(i).Robs_rating = [];
    data(i).blockID = blockID;
    data(i).block_inds = block_inds;
    data(i).cloud_area = cloud_area;
    data(i).cloud_binary = cloud_binary;
    data(i).cloud_scale = cloud_scale;
    data(i).datafilts = datafilts;
    data(i).datafiltsMU = datafiltsMU;
    data(i).disc_diameter = disc_diameter; % added 8/4/26 by mjg
    data(i).dt = dt;
    data(i).electrode_info = electrode_info;
    data(i).exptdate = exptdate;
    data(i).exptname = exptname;
    data(i).fix_location = fix_location;
    data(i).fix_size = fix_size;
    data(i).pixel_size = pixel_size;
    data(i).sacc_inds = sacc_inds;
    data(i).spikeIDs = spikeIDs;
    data(i).spike_ts = spike_ts_sorted;
    data(i).stim = stim;
    data(i).stim2 = stim2;
    data(i).stim3 = stim3;
    data(i).stim_area = stim_area;
    data(i).stim_location = stim_location;
    data(i).stim_location_deltas = stim_location_deltas;
    data(i).stimscale = stimscale;
    data(i).stimtype = stimtype;
    data(i).trialID = trialID;
    data(i).trial_start_ts = trial_start_ts;
    data(i).useLeye = useLeye;
    data(i).useReye = useReye;
    data(i).reward_on_ts = reward_on_ts;
    data(i).reward_off_ts = reward_off_ts;
    data(i).valid_data = valid_data;
    data(i).primaryStimType = trialLabels{i};


end

%% modify array_labels so its python readable
% unique_array_labels = unique(array_labels);
%
% [C, IA, IC] = unique(array_labels);
% data(i).arrayPerSU = IC(data(i).spikeSortingBatch);
% data(i).arrayPerMU = IC(data(i).spikeSortingBatchMU);
%
% %array_labels = cellfun(@(num,lab) [num2str(num) lab], num2cell(1:num_ks_batch), array_labels, 'UniformOutput', false);
% %data(i).array_labels = horzcat(array_labels{:});
% %data(i).array_labels = array_labels;
% unique_labels = horzcat(cellfun(@(num,lab) ['_' num2str(num) lab], num2cell(1:length(C)), unique_array_labels, 'UniformOutput', false));
% data(i).array_labels = horzcat(unique_labels{:});
%
% data(i).cluster = cluster;
% data(i).clusterMU = clusterMU;


%% Process LFPs
numDigitsInLastSpkChan = ceil(log10(length(pl2.SpikeChannels)));
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

%end