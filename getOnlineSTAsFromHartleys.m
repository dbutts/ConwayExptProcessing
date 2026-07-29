function [STA, Robs_strct] = getOnlineSTAsFromHartleys(plexon_fname, mainKofiko_fname, hartleySize)

stimpath = '/Volumes/lsr-conway/PROJECTS/V1_Fovea/stimuli/Cloudstims_calib_04_2024';
codedir = '/Users/greenemj/Git';
[a,filenameP,~] = fileparts(mainKofiko_fname);

kofiko_subfolder = fullfile(a,filenameP);

addpath(genpath(fullfile(codedir, 'ConwayExptProcessing')));

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

g_strcts = loadKofikoStructs(mainKofiko_fname);
trial = loadKofikoTrialData(kofiko_subfolder,mainKofiko_fname,filenameP);

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
%PlexET_ad_calib = loadPlexonEyeData(plexon_fname, chanNums, rig, plexonAnalogScale, KofikoGains_Plexon);

%% Trial analysis

% Determine stimulus intervals
stimStartTimes = [ones(size([trial.m_fImageFlipON_TS_Kofiko]')),...
    [trial.m_fImageFlipON_TS_Kofiko]']*B;

stimStopTimes = stimStartTimes + [trial.m_fStimulusON_MS]'/1e3;

stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames =  transpose(min([trial.numFrames], [trial.numFrames] ./ [trial.repframes]));

stimTiming.stimStartTimes = stimStartTimes;
stimTiming.stimStopTimes = stimStopTimes;
stimTiming.numFrames = numFrames;

X_fixationSpot = cellfun(@(x) x(1), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);

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


goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < maxFixationErrorPix)./numel(x)...
    > minFixationDuration, Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < maxFixationErrorPix)./numel(y) ...
    > minFixationDuration, Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

goodFixationIdx = vertcat(trial.m_bMonkeyFixated) | (goodFixationX & goodFixationY);

%% Bin stimulus sequences by trial

[stim1_cellArray, stim2_cellArray, stim3_cellArray] = makeStimMatrix(stimpath, trial, LumScale);

% find dualstim elemetns
trialTypeOfInterest = 'Dual Stim';
trialTypeOfInterestIdx = strcmpi( {trial.m_strTrialType}, trialTypeOfInterest)';

trlonset_diffs = [4; diff(stimStartTimes)];
areaIdx = cellfun(@(x) x==hartleySize, {trial.m_aiStimulusArea})';
modalStimRect = mode(cat(1,trial.m_aiStimulusRect));

modalStimRectIdx = all(cat(1,trial.m_aiStimulusRect) - modalStimRect == 0, 2);

isTrialOfInterest = trialTypeOfInterestIdx & ...
    goodFixationIdx &...
    trlonset_diffs > 4 &...
    areaIdx & ...
    vertcat(trial.DualstimPrimaryuseRGBCloud) == 6 &... % Hartleys
    cellfun(@(x) x>0, {trial.m_aiStimulusArea})';

assert(sum(isTrialOfInterest > 0));

stimulus_matrix = horzcat(stim1_cellArray{isTrialOfInterest});
stimulusET_matrix = horzcat(stim2_cellArray{isTrialOfInterest});

spkData = organizeSpikeDataByTrial(stimIntervals,plexon_fname, minSpikes);
Robs_strct = buildRobs(spkData.spkDataOnline, stimTiming, isTrialOfInterest);

stim = reshape(stimulus_matrix, hartleySize,hartleySize,3,[]);

STA = generate_stas(Robs_strct.Robs, stim,nLags);


end