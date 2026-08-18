function data = makeDataStruct(filenameP, isTrialOfInterest,g_strcts, trial, ...
    PlexET_cellArrays, KofikoET_cellArrays, stimTiming, stim_mats,...
    spkDataXline, Robs_strct, events, ETdist_thresh, trialLabel, computerLocation)

%% %%%%%%%%%%%% Format output like PackageCloud %%%%%%%%%%%%
% Get trial IDs (insane that there are trial numbers and IDs....)

% Dual stim specific

if any(strcmpi(trialLabel, {'ccloud', 'hartley'}))

TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), {trial.TrialID}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), {trial.usebinary}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

try
    spatialscale = {trial.spatialscale}';

catch
    spatialscale_buffer = g_strcts.g_astrctAllParadigms{1}.DualstimScale.Buffer;
    spatialscale_TS = g_strcts.g_astrctAllParadigms{1}.DualstimScale.TimeStamp;
    spatialscale = nan(size(trial));
    for i = 1:numel(spatialscale_buffer)
        spatialscale(ImageFlipON_TS_Kofiko >= spatialscale_TS(i)) = spatialscale_buffer(i);
    end

end

BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.BlockID}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
StimulusAreaPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.m_aiStimulusArea}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
DualstimPrimaryuseRGBCloudPerFrame = cellfun(@(x, y) repelem(x, y),{trial.DualstimPrimaryuseRGBCloud}', num2cell(stimTiming.numFrames), 'UniformOutput', false);

else
    TrialIDPerFrame = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
    useBinaryPerFrame = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
    spatialscale = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
    BlockIDPerFrame = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
    StimulusAreaPerFrame = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
    DualstimPrimaryuseRGBCloudPerFrame = num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));

end

% Disc probe specific
if strcmpi(trialLabel, 'discProbe')
DiscDiameterPerFrame =  cellfun(@(x, y) repelem(x, y), {trial.m_iDiscDiameter}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
else
    DiscDiameterPerFrame =num2cell(nan(1, sum(stimTiming.numFrames(isTrialOfInterest))));
end



X_fixationSpot = cellfun(@(x) x(1), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);

UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseLeye}', num2cell(stimTiming.numFrames), 'UniformOutput', false);
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), {trial.UseReye}', num2cell(stimTiming.numFrames), 'UniformOutput', false);


Kofiko_Xpix_frameRate_cellArray = ...
    cellfun(@(t,x, start, stop, fixspot, n) interp1(t, x - fixspot, linspace(start, stop, n), 'linear'),...
    KofikoET_cellArrays.Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    KofikoET_cellArrays.Kofiko_Xpix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimTiming.stimStartTimes(isTrialOfInterest)),...
    num2cell(stimTiming.stimStopTimes(isTrialOfInterest)),...
    X_fixationSpot(isTrialOfInterest),...
    num2cell(stimTiming.numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);

Kofiko_Ypix_frameRate_cellArray = ...
    cellfun(@(t,y, start, stop, fixspot, n) interp1(t, y - fixspot, linspace(start, stop, n), 'linear'),...
    KofikoET_cellArrays.Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    KofikoET_cellArrays.Kofiko_Ypix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimTiming.stimStartTimes(isTrialOfInterest)),...
    num2cell(stimTiming.stimStopTimes(isTrialOfInterest)),...
    Y_fixationSpot(isTrialOfInterest),...
    num2cell(stimTiming.numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);

%ETgains
ETgains = [g_strcts.g_strctEyeCalib.GainX.Buffer(end), g_strcts.g_strctEyeCalib.GainY.Buffer(end)];

%ETstim_location
ETstim_location =  [trial(find(isTrialOfInterest,1, 'last')).secondarystim_bar_rect;...
    trial(find(isTrialOfInterest,1, 'last')).tertiarystim_bar_rect];

% ETtrace
ETtrace = [[Kofiko_Xpix_frameRate_cellArray{:}]; [Kofiko_Ypix_frameRate_cellArray{:}]];

%ETtrace_raw
ETtrace_raw = transpose([vertcat(PlexET_cellArrays.eyeX2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(PlexET_cellArrays.eyeY2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(PlexET_cellArrays.eyeX1_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(PlexET_cellArrays.eyeY2_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(PlexET_cellArrays.pupil2_plexon_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(PlexET_cellArrays.pupil1_plexon_cellArray{2*find(isTrialOfInterest)})]);

%Robs

nSU = numel(Robs_strct.SU_clusters);
nMU = numel(Robs_strct.MU_clusters);

RobsSU = Robs_strct.Robs(1:nSU,:);

%RobsMU
RobsMU = Robs_strct.Robs(nSU+1:end,:);

%RobsMU_probe_ID
RobsMU_probe_ID = Robs_strct.MU_chans;

%RobsMU_rating
RobsMU_rating = [];
%blockID

% Robs_probe_ID
Robs_probe_ID = Robs_strct.SU_chans;

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
clusterIDs= cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), {spkDataXline.spk_clusters_cellArray}, 'UniformOutput', false);
clusterIDs = transpose(vertcat(clusterIDs{:}));
clusterIDs = vertcat(clusterIDs);

assert(numel(unique(clusterIDs)) == nSU + nMU);

%spike_ts
spike_ts_raw = cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), {spkDataXline.spk_times_cellArray}, 'UniformOutput', false);
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

stim = stim_mats.stim1_matrix;
if any(strcmpi(trialLabel, {'ccloud', 'hartley'}))
    stim = reshape(stim, stimH,stimW,3,[]);
end

stim2 = stim_mats.stim2_matrix;
stim3 = stim_mats.stim3_matrix;

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

uniqueClusterIDs = [Robs_strct.SU_clusters; Robs_strct.MU_clusters];

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
data.disc_diameter = disc_diameter; % added 8/4/26 by mjg
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
data.stim2 = stim2;
data.stim3 = stim3;
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
data.primaryStimType = trialLabel;
data.SU_clusters = Robs_strct.SU_clusters;
data.MU_clusters = Robs_strct.MU_clusters;
data.SU_chans = Robs_strct.SU_chans;
data.MU_chans = Robs_strct.MU_chans;
data.SU_arrays= Robs_strct.SU_arrays;
data.MU_arrays = Robs_strct.MU_arrays;

end