% directory containing Kofiko mat files
dirpath ='/home/bizon/Data/V1_Fovea/Sprout/260205';
filenameP = '260205_153540_Sprout';

addpath(genpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/'));
addpath(genpath('/home/bizon/Processing'));


%% Read in plexon data

plexon_fname = fullfile(dirpath, [filenameP '.pl2']);
% plexon ddpi sync signal
[adfreq, n, ts, fn, sync_ad] = plx_ad_v(plexon_fname, 'AI01'); % ddpi sync

[~, ~, ~, ~, lefPupil_plexon] = plx_ad_v(plexon_fname, 'AI03'); % L pupil
[~, ~, ~, ~, rightPupil_plexon] = plx_ad_v(plexon_fname, 'AI04'); % R pupil
[~, ~, ~, ~, rightEyeX_plexon] = plx_ad_v(plexon_fname, 'AI05'); % R X, used for tracking
[~, ~, ~, ~, rightEyeY_plexon] = plx_ad_v(plexon_fname, 'AI06'); % R Y used for tracking
[~, ~, ~, ~,leftEyeX_plexon] = plx_ad_v(plexon_fname, 'AI07'); % L X
[ET_adfreq, ET_n, ET_ts, ET_fn, leftEyeY_plexon] = plx_ad_v(plexon_fname, 'AI08'); % L Y

t_plexon = 0:(1/adfreq):(length(sync_ad)-1)/adfreq;
%% load kofiko data

% find the mat files produced by Kofiko

kofiko_dir = dir(fullfile(dirpath, '**', [filenameP '*.mat']));
kofiko_fname_re = [filenameP '(_\d*)?.mat'];
valid_kofiko_fname_idx = cellfun(@(x) ~isempty(x), ...
    cellfun(@(x) regexp(x, kofiko_fname_re), {kofiko_dir(:).name}, 'UniformOutput',false));

kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
kofiko_fname = {kofiko_dir(valid_kofiko_fname_idx).name};
%sort

kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'),...
    kofiko_fname, 'UniformOutput',false);
kofiko_file_numStr{cellfun(@isempty, kofiko_file_numStr)} = {'nan'};
kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);
kofiko_file_num(isnan(kofiko_file_num)) = max(kofiko_file_num)+1;

% load the g_strctLocalExperimentRecording structs which contain
% information about each trials
g_strctLocalExperimentRecording = [];
for i = 1:length(kofiko_file_num)
    idx = kofiko_file_num(i);

    temp_full_fname = fullfile(kofiko_folder{i}, kofiko_fname{i});
    temp_strct = load(temp_full_fname, 'g_strctLocalExperimentRecording');
    temp_cell = temp_strct.g_strctLocalExperimentRecording;
    g_strctLocalExperimentRecording{idx} = temp_cell;
end

% load the file which should be the "main" Kofiko mat file (no number at
% the end)

load(fullfile(kofiko_folder{1}, kofiko_fname{1}), 'g_strctDAQParams', ...
    'g_astrctAllParadigms', 'g_strctLog', 'g_strctEyeCalib', ...
    'g_strctAppConfig', 'g_strctSharedParadigmData', ...
    'g_strctStimulusServer', 'g_strctSystemCodes');


% concatenate the g_strctLocalExperimentRecorindg structs and remove empty
% cells
g_strctLocalExperimentRecording = vertcat(g_strctLocalExperimentRecording{:});
g_strctLocalExperimentRecording(cellfun(@isempty, g_strctLocalExperimentRecording)) = [];

% unify field names across g_strctLocalExperimentRecording

fieldNames = cellfun(@(x) fieldnames(x), g_strctLocalExperimentRecording, 'UniformOutput', false);
uniqueFieldNames = unique(vertcat(fieldNames{:}));
missingFieldNames = cellfun(@(x) setdiff(uniqueFieldNames, fieldnames(x)),...
    g_strctLocalExperimentRecording, 'UniformOutput',false);

for i = 1:numel(g_strctLocalExperimentRecording)
    for  j = 1:numel(missingFieldNames{i})
        g_strctLocalExperimentRecording{i}.(missingFieldNames{i}{j}) = [];
    end
end


%% Time alignment

% read kofiko strobes from plexon and kofiko files
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);
kofikoSyncStrobesTS = transpose(g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == g_strctSystemCodes.m_iSync));
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iSync);

% linear regression to get kofiko time stamps into plexon time

B = [ones(size(kofikoSyncStrobesTS)) kofikoSyncStrobesTS]\plexonSyncStrobesTS;


%% Kofiko eye data

% get kofiko sync timestamps from control computer

% in seconds
t_kofiko_sync_control = g_strctDAQParams.LastStrobe.TimeStamp(...
    g_strctDAQParams.LastStrobe.Buffer == g_strctSystemCodes.m_iSync);

% in seconds
t_kofiko_sync_plexon = events.timeStamps(events.strobeNumber==g_strctSystemCodes.m_iSync);

b_control_plexon = [ones(size(t_kofiko_sync_control')) t_kofiko_sync_control']\t_kofiko_sync_plexon;

% kofiko eye signal
ScreenSizeX_pix = g_strctStimulusServer.m_aiScreenSize(3);
ScreenSizeY_pix = g_strctStimulusServer.m_aiScreenSize(4);

Kofiko_EyeRawX = g_strctEyeCalib.EyeRaw.Buffer(:,1);
Kofiko_EyeRawY = g_strctEyeCalib.EyeRaw.Buffer(:,2);
Kofiko_ET_TS = g_strctEyeCalib.EyeRaw.TimeStamp';
Kofiko_ET_TS_PlexonTime = [ones(size(Kofiko_ET_TS)) Kofiko_ET_TS]*B;

Kofiko_GainX = g_strctEyeCalib.GainX.Buffer;
Kofiko_GainY = g_strctEyeCalib.GainY.Buffer;

Kofiko_GainX_TS = g_strctEyeCalib.GainX.TimeStamp;
Kofiko_GainY_TS = g_strctEyeCalib.GainY.TimeStamp;

Kofiko_CenterX = g_strctEyeCalib.CenterX.Buffer;
Kofiko_CenterY = g_strctEyeCalib.CenterY.Buffer;
Kofiko_CenterX_TS = g_strctEyeCalib.CenterX.TimeStamp;
Kofiko_CenterY_TS = g_strctEyeCalib.CenterY.TimeStamp;

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

% The way to convert raw eye signal to pixel coordinates is:
% [RawSignal-SignalOffset]*Gain + ScreenCenter = FixationPoint
% SignalOffset is CenterX or CenterY

Kofiko_Xpix = Kofiko_GainX_forEachTimeStamp.*(Kofiko_EyeRawX - Kofiko_CenterX_forEachTimeStamp) + ScreenSizeX_pix/2;
Kofiko_Ypix = Kofiko_GainY_forEachTimeStamp.*(Kofiko_EyeRawY - Kofiko_CenterY_forEachTimeStamp) + ScreenSizeY_pix/2;


% we can compare the kofiko and plexon eye signals for alignment as well




%% trial analysis

% order trials in g_strctLocalExperimentRecording
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, g_strctLocalExperimentRecording);
[~, trial_idx_for_sorting] = sort(tempImageFlipON_TS_Kofiko , 'ascend');

% trialData cell array
trialData = g_strctLocalExperimentRecording(trial_idx_for_sorting);
trialNum = cellfun(@(x) x.m_iTrialNumber, trialData);

% usebinary: 0 = high contrast cloud, 1 = binary cloud, 2 = contrast
% matched cloud
usebinary = cellfun(@(x) x.usebinary, trialData, 'UniformOutput',false);

% cloud spatial scale
spatialscale = cellfun(@(x) x.spatialscale, trialData, 'UniformOutput', false);

% stimulus blocks
BlockID = cellfun(@(x) x.BlockID, trialData, 'UniformOutput', false);

% matrix of unique conditions defined by cloud parameters and block
uniqueCloudConditions =  unique([[usebinary{:}]' [spatialscale{:}]' [BlockID{:}]'], 'rows');

usebinary(cellfun(@isempty, usebinary)) = {nan};
spatialscale(cellfun(@isempty, spatialscale)) = {nan};
BlockID(cellfun(@isempty, BlockID)) = {nan};

trialType = cellfun(@(x) x.m_strTrialType, trialData, 'UniformOutput',false);

% determine stimulus intervals
StimulusON_MS = cellfun(@(x) x.m_fStimulusON_MS, trialData);
StimulusOFF_MS = cellfun(@(x) x.m_fStimulusOFF_MS, trialData);
ImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
stimStartTimes = [ones(size(ImageFlipON_TS_Kofiko)), ImageFlipON_TS_Kofiko]*B; % plexon time
stimStopTimes = stimStartTimes + StimulusON_MS./1e3;
stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% stim seq on each trial
stimseq = cellfun(@(x) x.stimseq, trialData, 'UniformOutput',false);
numFrames = cellfun(@(x) x.numFrames, trialData);
fps = 60;
repFrames = numFrames./(fps*StimulusON_MS*1e-3);
repFrames(~isfinite(repFrames)) = 0;
numFrames = numFrames./repFrames;
numFrames(~isfinite(numFrames)) = 0;

% get actual sequence accounting for repFrames
stimseq= cellfun(@(x,y) x(1:y:end), stimseq, num2cell(repFrames), 'UniformOutput', false);


% bin eye data by trial
% each bin tells you which stim interval ts belongs to
% only odd numbered bins, correspond to stim on, other bins are period
% between stim offset and next trial


% To get kofiko eye date for trial trialNum(t):
% Kofiko_ET_TS_PlexonTime(Kofiko_ET_TS_PlexonTime_Bin == stimON_Kofiko_ET_TS_PlexonTime_Bin(t))

[~,~,Kofiko_ET_TS_PlexonTime_Bin] = histcounts(Kofiko_ET_TS_PlexonTime, stimIntervals); % wont really need this ever, I think
uniqueKofiko_ET_TS_PlexonTime_Bin= unique(Kofiko_ET_TS_PlexonTime_Bin);
%stimON_Kofiko_ET_TS_PlexonTime_Bin = uniqueKofiko_ET_TS_PlexonTime_Bin(1:2:end);

% % get trial number for each time bin
% trialNumPerKofiko_ET_TS_PlexonTime_Bin =


%[~,~,bin_t_dpi_plexon] = histcounts(t_dpi_plexon, stim_intervals_plexon);
%[~,~,bin_t_kofiko_plexon] = histcounts(t_kofiko_plexon, stim_intervals_plexon);

% t_dpi_plexon
% which bins correspond to which trials?
% bin 0 should be times not within range of stim_intervals_plexon
% bin 1 should be times within first stim interval. If it's empty, dpi
% wasnt recording for this stim presentation.


Kofiko_ET_TS_PlexonTime_cellArray =accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_ET_TS_PlexonTime(:), [], @(x){x});
Kofiko_Xpix_cellArray = accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_Xpix(:), [], @(x){x});
Kofiko_Ypix_cellArray = accumarray(Kofiko_ET_TS_PlexonTime_Bin(:)+1, Kofiko_Ypix(:), [], @(x){x});

% Kofiko_Xpix_mean = cellfun(@mean, Kofiko_Xpix_cellArray(2:2:end));
% Kofiko_Xpix_iqr = cellfun(@iqr, Kofiko_Xpix_cellArray(2:2:end));
% Kofiko_Xpix_iqr(isnan(Kofiko_Xpix_iqr)) = Inf;
%
% Kofiko_Ypix_mean = cellfun(@mean, Kofiko_Ypix_cellArray(2:2:end));
% Kofiko_Ypix_iqr = cellfun(@iqr, Kofiko_Ypix_cellArray(2:2:end));
% Kofiko_Ypix_iqr(isnan(Kofiko_Ypix_iqr)) = Inf;
%
% goodFixationIdx = Kofiko_Xpix_iqr < 20 & Kofiko_Ypix_iqr < 20 & abs(Kofiko_Xpix_mean-ScreenSizeX_pix/2) < 30 & abs(Kofiko_Ypix_mean-ScreenSizeY_pix/2) < 30;

fix_thresh = 0.6;
MonkeyFixated = cellfun(@(x) logical(x.m_bMonkeyFixated), trialData);
goodFixationX = cellfun(@(x) sum(abs(x - ScreenSizeX_pix/2) < 45)./numel(x) > fix_thresh, Kofiko_Xpix_cellArray(2:2:end)) ;
goodFixationY = cellfun(@(x) sum(abs(x - ScreenSizeY_pix/2) < 45)./numel(x) > fix_thresh, Kofiko_Ypix_cellArray(2:2:end)) ;

goodFixationIdx = MonkeyFixated | (goodFixationX & goodFixationY);

%% load in spike data

load('/home/bizon/Data/V1_Fovea/Sprout/260205/preprocessingInfo.mat')

% drop test chech
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

ii = 1;
for a = 1:numel(arrayLabels)
    for c = 1:numel(curChannels{a})
        if contains(arrayLabels{a}, 'lam')
            ksFolders{ii} = fullfile(dirpath, filenameP, 'kilosorting_laminar');
            ks_nChans(ii) = length(str2num(curChannels{a}{c}));
            ks_arraylabels_rep{ii} = arrayLabels{a};
        else
            ksFolders{ii} = ['kilosorting_' arrayLabels{a} '_' num2str(min(str2num(curChannels{a}{c}))) 'to' num2str(max(str2num(curChannels{a}{c})))];
            ksFolders{ii} = fullfile(dirpath, filenameP, ksFolders{ii});
            ks_nChans(ii) = length(str2num(curChannels{a}{c}));
            ks_arraylabels_rep{ii} = arrayLabels{a};
            ii = ii + 1;
        end
    end
end

ks.filepath = ksFolders; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
%ks.filepath = [dirpath filenameP filesep 'kilosorting_stitched' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.arraylabels_rep = ks_arraylabels_rep;
ks.ks_nChans = ks_nChans;

for f = 1:numel(ks.filepath)
    spk_times{f} = readNPY(fullfile(ks.filepath{f}, 'spike_times_seconds.npy')) + opts.spk_offset;
    spk_clusters{f} = readNPY(fullfile(ks.filepath{f}, 'spike_clusters.npy')); % zero indexed
    spk_info{f} = tdfread(fullfile(ks.filepath{f}, 'cluster_info.tsv'));
end

% get array label for each individual channel
arrayLabelPerChan = repelem(ks.arraylabels_rep, ks.ks_nChans);
arrayLabelPerCluster = repelem(ks.arraylabels_rep, cellfun(@numel, spk_clusters));
plexonChanNum= 1:sum(ks.ks_nChans);
ksChanNum = arrayfun(@(x) 0:x, ks.ks_nChans-1, 'UniformOutput', false);
ksChanNum = horzcat(ksChanNum{:});

ks_spk_clusters = vertcat(spk_clusters{:});
% concatenate these cell arrays, carefully
% each cluster must be unique, so add maximum cluster number of
% previous cell
max_cluster_nums = cellfun(@max, spk_clusters);
max_cluster_nums_cumsum = cumsum(max_cluster_nums + 1);
max_cluster_nums_cumsum = [0 max_cluster_nums_cumsum(1:end-1)];

ks_chnOffset = [0 cumsum(ks.ks_nChans(1:end-1))];
for f = 1:numel(ks.filepath)
    spk_clusters{f} = spk_clusters{f} + max_cluster_nums_cumsum(f);
    spk_info{f}.cluster_id = unique(spk_clusters{f});
    spk_info{f}.ch = spk_info{f}.ch + ks_chnOffset(f);
end

spk_times = vertcat(spk_times{:});
spk_clusters = vertcat(spk_clusters{:});

% now concatenate fields of spk_info
spk_info_temp = [spk_info{:}];
spk_info_fields = fieldnames(spk_info_temp);
spk_info_new = struct;

for i = 1:numel(spk_info_fields)
    spk_info_new.(spk_info_fields{i}) = vertcat(spk_info_temp.(spk_info_fields{i}));
end

spk_info = spk_info_new;

spk_clustIDs = unique(spk_clusters);

% get indices of good, mua, noise

goodStr = 'good ';
muaStr = 'mua  ';
noiseStr = 'noise';
blankStr = '     ';

SU_IDs = spk_info.cluster_id(sum(spk_info.group == goodStr, 2) == numel(goodStr) & spk_info.n_spikes > 1000);
MU_IDs = spk_info.cluster_id((sum(spk_info.group == muaStr, 2) == numel(muaStr) | sum(spk_info.group == blankStr, 2) == numel(blankStr))  & spk_info.n_spikes > 1000);

SU_MU_IDs =[SU_IDs; MU_IDs];

% bin spikes by trial
[~,~,spk_times_Bin] = histcounts(spk_times, stimIntervals);
spk_times_cellArray =accumarray(spk_times_Bin(:)+1, spk_times(:), [], @(x){x});

% IMPORTANT: spk_times_Bin will work on spk_clusters as well!!!

spk_clusters_cellArray = accumarray(spk_times_Bin(:)+1, spk_clusters(:), [], @(x){x});
% for each unique cloud condition (binary, spatial scale, blockID) load the
% clouds, and organize stimuli into array where each cell is a trial

%% stimulus

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

% to get back stimulus j of trial i into 60 x 60 x 3 form, do reshape(stimulus_cellArray{i}(:,j), [60 60 3]);

% now, prepare bins to group spike times by frame

stimFrameBins = cellfun(@(x, y, z) linspace(x,y,z+1), num2cell(stimStartTimes), num2cell(stimStopTimes), num2cell(numFrames), 'UniformOutput', false);

spksPerFrame_cellArray =  cell(size(trialData,1),1);

% element of each cell is the frame on which spike occured -- cluster to
% which spike belongs given by spk_cluster_cellArray;

valid_stimFrameBinIdx = cellfun(@(x) numel(x)>=2, stimFrameBins);

[~,~,spkFrameIdx_cellArray(valid_stimFrameBinIdx)] = cellfun(@(x,y) histcounts(x, 'BinEdges', y), spk_times_cellArray(2*find(valid_stimFrameBinIdx)), stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

%nSpkPerFrame_cellArray(valid_stimFrameBinIdx) = cellfun(@(x,y) histcounts(x, 'BinEdges', y), spk_times_cellArray(2*find(valid_stimFrameBinIdx)), stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);

% find dualstim elemetns

dualStimIdx = strcmpi(trialType, 'Dual Stim');
DualstimPrimaryuseRGBCloud=cellfun(@(x) x.DualstimPrimaryuseRGBCloud, trialData, 'UniformOutput',false);
StimulusArea = cellfun(@(x) x.m_aiStimulusArea, trialData, 'UniformOutput', false);

highContrastIdx = cellfun(@(x) x==0, usebinary);

trlonset_diffs = [4; diff(stimStartTimes)];

trialsOfInterestIdx = dualStimIdx & goodFixationIdx & trlonset_diffs > 4 &...
    cellfun(@(x) any(x>0), StimulusArea) & cellfun(@(x) any(x==8), DualstimPrimaryuseRGBCloud);
% horizontally concatenate stimuli

stimulus_matrix = horzcat(stimulus_cellArray{trialsOfInterestIdx});


%stimulus_matrix = reshape(stimulus_matrix, 60, 60, 3, []);

% to get stas multiply simulus matrix by matrix of spike trains

Robs = single(zeros(numel(SU_MU_IDs),  sum(numFrames(trialsOfInterestIdx))));

for i = 1: numel(SU_MU_IDs)
    unitID = SU_MU_IDs(i);

    % each cell is a trial, and contains frame indicies where spikes
    % occurred (repeated frames == multiple spikes on that frame)

    spksPerFrameThisUnit_cellArray = cellfun(@(x,y) x(y == unitID), spkFrameIdx_cellArray, spk_clusters_cellArray(2:2:end) ,'UniformOutput', false);

    %compute spiks per bin (frame)

    % computes spikes/bin (frame) on trials where its possible (at least 2
    % frames)
    spksPerFrame_cellArray(trialsOfInterestIdx) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
        spksPerFrameThisUnit_cellArray(trialsOfInterestIdx), num2cell(numFrames(trialsOfInterestIdx)) ,'UniformOutput',false);

    Robs(i,:) = [spksPerFrame_cellArray{:}];

end


disp('Calculating and plotting STAst')
for i = 1:size(Robs,1)
    r = transpose(Robs(i,:));
    spkTimes = find(r);
    nSpks = numel(spkTimes);
    nLags = 10;
    maxLag = nLags - 1;

    tempSTA = zeros(size(stimulus_matrix,1), nLags, 'like', stimulus_matrix);
    countsPerLag = zeros(1,nLags);

    for j = 1:nSpks
        t = spkTimes(j);
        if t <= maxLag
            continue;
        end
        cols = (t - (0:maxLag));
        tempSTA = single(tempSTA) + (single(stimulus_matrix(:,cols))./127)*r(t);
        countsPerLag = countsPerLag + r(t);
    end
    countsPerLag(countsPerLag==0)=1;
    tempSTA = bsxfun(@rdivide, single(tempSTA), single(countsPerLag));

    STA(i,:,:) = tempSTA;
    %STA(i,:,:,:,:) = reshape(tempSTA, 60,60,3, nLags);

    % A = reshape(tempSTA(:,4), 60, 60, 3);
    % A = circshift(A, [30 30]);
    % figure, imagesc(A(:,:,1)); colormap gray
end

STA = reshape(STA, size(STA,1), 60, 60, 3, nLags);


% for SU_MU_IDs(13);
i =13;
PhyClus = unqiue(ks_spk_clusters(spk_clusters == SU_MU_IDs(i)));
PhyChan =  spk_info.ch(spk_info.cluster_id == SU_MU_IDs(i));