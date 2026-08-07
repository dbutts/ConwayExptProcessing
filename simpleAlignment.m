
%function data = simpleAlignment
%%
% 3/2026 mjg -- wrote it
% 7/2026 mjg -- major refactoring

%% Set paths
setPathsAndFlags;
addpath(genpath(fullfile(codedir, 'ConwayExptProcessing')))

%% Initialization for everyone given setup

pl2 = PL2ReadFileIndex(plexon_fname);
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);
%% Hardcoded values
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

KofikoET_cellArrays.Kofiko_ET_TS_PlexonTime_cellArray = binByStimIntervals(Kofiko_ET_TS_PlexonTime,  Kofiko_ET_TS_PlexonTime, stimTiming.stimIntervals);
KofikoET_cellArrays.Kofiko_Xpix_cellArray = binByStimIntervals(Kofiko_ET_TS_PlexonTime,  Kofiko_Xpix, stimTiming.stimIntervals);
KofikoET_cellArrays.Kofiko_Ypix_cellArray = binByStimIntervals(Kofiko_ET_TS_PlexonTime, Kofiko_Ypix, stimTiming.stimIntervals);


%% Bin Plexon eye signal by stimulus intervals

PlexET_cellArrays.t_plexon_cellArray = binByStimIntervals(t_plexon, t_plexon, stimTiming.stimIntervals);

% pupils
PlexET_cellArrays.pupil1_plexon_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,3), stimTiming.stimIntervals);
PlexET_cellArrays.pupil2_plexon_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,4), stimTiming.stimIntervals);
% formerly "right eye"
PlexET_cellArrays.eyeX2_plexon_calib_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,5), stimTiming.stimIntervals);
PlexET_cellArrays.eyeY2_plexon_calib_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,6), stimTiming.stimIntervals);
% formerly "left eye"
PlexET_cellArrays.eyeX1_plexon_calib_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,7), stimTiming.stimIntervals);
PlexET_cellArrays.eyeY1_plexon_calib_cellArray = binByStimIntervals(t_plexon, PlexET_ad_calib(:,8), stimTiming.stimIntervals);


%% Determine trials with good fixation
X_fixationSpot = cellfun(@(x) x(1), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), {trial.m_pt2iFixationSpot}', 'UniformOutput',false);

goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < maxFixationErrorPix)./numel(x)...
    > minFixationDuration, KofikoET_cellArrays.Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < maxFixationErrorPix)./numel(y) ...
    > minFixationDuration, KofikoET_cellArrays.Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

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
    stim_mats.stim1_matrix = horzcat(stim1_cellArray{isTrialOfInterest});
    stim_mats.stim2_matrix = horzcat(stim2_cellArray{isTrialOfInterest});
    stim_mats.stim3_matrix = horzcat(stim3_cellArray{isTrialOfInterest});

    %% resample eye signal at frame rate
    % NOTE: Without 'extrap' argument, there will be NANs:

    % Robs
    if isfield(spkData, 'spkDataOffline')
        Robs_strct(i) = buildRobs(spkData.spkDataOffline, stimTiming, isTrialOfInterest);
        data(i) = makeDataStruct(filenameP,isTrialOfInterest,g_strcts, trial, PlexET_cellArrays, KofikoET_cellArrays, stimTiming, stim_mats, spkData.spkDataOffline,Robs_strct(i), events, ETdist_thresh, trialLabels{i}, computerLocation);
    else

    end

     Robs_strct_online(i) = buildRobs(spkData.spkDataOnline, stimTiming, isTrialOfInterest);
     data_online(i) = makeDataStruct(filenameP, isTrialOfInterest,g_strcts, trial, PlexET_cellArrays, KofikoET_cellArrays, stimTiming, stim_mats, spkData.spkDataOnline,Robs_strct_online(i), events, ETdist_thresh, trialLabels{i}, computerLocation);
    %% %%%%%%%%%%%% Format output like PackageCloud %%%%%%%%%%%%

end

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