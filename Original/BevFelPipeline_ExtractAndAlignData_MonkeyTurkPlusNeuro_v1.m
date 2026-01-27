%%
clear all
%%
filenames = {'210123_130251_Jacomo'}; iso_SUs=[5,6,8,9,10,13,16,17,21];
%filenames = {'210130_162547_Jacomo'}; iso_SUs=[2, 4:8, 11, 13, 14, 15, 19, 20, 24];
%filenames = {'210202_101158_Jacomo'}; iso_SUs=[8, 10, 17, 18, 19, 20, 23, 24];

eeday = 1;

%for eeday = 1 %can loop over multiple days to align behavior data
filenameP = filenames{eeday};
    
%% Extract Bevil's data
useofflinesorting = 1;
nChans=24; 
BehaviorOnly=0; %this flag skips aligning plexon data
skip_saccades=1; %this skips saccade extraction and just keeps the raw eye traces

sessionsToProcess = [];
ExperimentFolder = '/media/felix/Internal_1/Data/BevilColor/';
plxFilePath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '.pl2'];
matFilePath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/'];
sessionTimeOffsets = [];                                     
strExperimentPath = [matFilePath '/Analysis/'];

if ~exist(strExperimentPath,'dir');
    mkdir(strExperimentPath);
end
filenameK = matFilePath;
iSessions = 1;%

%%
global strctColorValues unitsInFile PlottingVars g_strctStatistics ExperimentRecording topLevelIndex 
ExperimentRecording=[];
%persistent plexonDataAlignedToThisTrial


experimentIndex = {};

currentDirectory = pwd;
cd(ExperimentFolder);

allPLXfiles = dir([pwd,filesep, '*.plx']);
g_strctStatistics.trackedXYZCoordinates = [];
g_strctStatistics.m_aiXYZCoorindateColorMatchingTable = [];
g_strctStatistics.preTrialWindow = -.5;
g_strctStatistics.postTrialWindow = 6.5; %this includes all parts of each trial
g_strctStatistics.bPlotRFFieldWithinStimRectangle = 1;
g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod = [-.050,2.2];


g_strctStatistics.numBins = 200;
g_strctStatistics.positionBinGrid = [12,20];
g_strctStatistics.bPlotRFFieldWithinStimRectangle = 1;
%g_strctStatistics.RFMappingTrialIntegrationTime = [.050,.250];
%g_strctStatistics.RFMappingTrialIntegrationTime = [];
g_strctStatistics.RFMappingTrialIntegrationTime(:,1) = -.05:.01:.34 ;
g_strctStatistics.RFMappingTrialIntegrationTime(:,2) = -.04:.01:.35 ;
g_strctStatistics.m_bPlotRFFieldByColor = 1;


PlottingVars.iNumPSTHBins = 200;
plotOrientation = 1;
numPositionBins = [7,5];
numOrientationBins = 20;
syncStrobeID = 32757;
eventChannelNumber = 257;
startRecordID = 32767;
stopRecordID = 32766;
startTrialID = 32700;
endTrialID = 32699;

%allMatFiles = dir([strExperimentPath, '*.mat']);
allMatFiles = dir('*.mat');
channel = 13;

%allMatFiles = what(strExperimentPath);
[~, indices] = sort(vertcat(allMatFiles(:).datenum));
allMatFiles = allMatFiles(indices);
%%
allMatFiles = dir([filenameK, '*.mat']);
[~, indices] = sort(vertcat(allMatFiles(:).datenum));
allMatFiles = allMatFiles(indices);
g_strctStatistics.ExperimentRecording = {};
for iFiles = 1:size(allMatFiles,1)
    if allMatFiles(iFiles).isdir
        continue
    end
    load([matFilePath,filesep,allMatFiles(iFiles).name])
    disp(sprintf('loading experiment file %s', allMatFiles(iFiles).name));
    if exist('g_strctLocalExperimentRecording') && size(g_strctLocalExperimentRecording{1},1) > 1
       warning('incorrect format detected in save structure in file  %s, skipping', allMatFiles(iFiles).name)
        %tmp = g_strctLocalExperimentRecording{1};
        continue
    end

    % details = whos([strExperimentPath,'\',allMatFiles(iFiles).name])
    %g_strctLocalg_strctStatistics.ExperimentRecording(cellfun('isempty',g_strctLocalg_strctStatistics.ExperimentRecording)) = [];
    g_strctStatistics.ExperimentRecording(cellfun('isempty',g_strctStatistics.ExperimentRecording)) = [];
    %g_strctLocalg_strctStatistics.ExperimentRecording = vertcat(g_strctLocalg_strctStatistics.ExperimentRecording,g_strctDynamicStimLog.TrialLog);
    %g_strctLocalg_strctStatistics.ExperimentRecording = vertcat(g_strctLocalg_strctStatistics.ExperimentRecording{cellfun(@isempty,g_strctLocalg_strctStatistics.ExperimentRecording),dataToSave');
    if (exist('g_strctLocalExperimentRecording') == 1)
        disp('case 1')
        g_strctStatistics.ExperimentRecording = vertcat(g_strctStatistics.ExperimentRecording,vertcat({g_strctLocalExperimentRecording{find(~cellfun(@isempty,g_strctLocalExperimentRecording))}})');
       clear g_strctLocalExperimentRecording
    elseif (exist('strctLocalExperimentRecording') == 1)
        disp('case 2')
        g_strctStatistics.ExperimentRecording = vertcat(g_strctStatistics.ExperimentRecording,vertcat({strctLocalExperimentRecording{find(~cellfun(@isempty,strctLocalExperimentRecording))}})');
        clear('strctLocalExperimentRecording')
    else
        disp('case 3')
        try
        g_strctStatistics.ExperimentRecording = vertcat(g_strctStatistics.ExperimentRecording,vertcat({dataToSave{find(~cellfun(@isempty,dataToSave))}})');
        catch
           warning(sprintf('file %s did not contain appropriately formatted trials', allMatFiles(iFiles).name)); 
        end
    end
end
ExperimentRecording = g_strctStatistics.ExperimentRecording;
clear g_strctStatistics.ExperimentRecording;


%%
%strctColorValues = fnCheckForColorFile();
cd(ExperimentFolder);

try
load(filenameK, 'g_strctDAQParams');
catch
    try
        load([ExperimentFolder, filenameP, '.mat'], 'g_strctDAQParams');
        
    catch
       [~, filename] = fileparts(filenameK); 
       load([filename, '.mat'], 'g_strctDAQParams');
    end
end

if ~BehaviorOnly
try
    [events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plxFilePath, eventChannelNumber);
%    [events] = PL2EventTs(plxFilePath, eventChannelNumber);
catch
    [~,sessionName] = fileparts(plxFilePath);
    sprintf('corrupt or missing information from Plexon file, experiment %s', sessionName)
    return;
end

firstStrobePlexonTS = events.timeStamps(find(events.strobeNumber(events.strobeNumber == syncStrobeID),1));
plexonStrobeIDX = find(events.strobeNumber == syncStrobeID);
lastStrobePlexonTS = events.timeStamps(find(events.strobeNumber(events.strobeNumber == syncStrobeID),1,'last'));
plexonStrobeAllTS = events.timeStamps(plexonStrobeIDX);
firstStrobeKofikoTS = g_strctDAQParams.LastStrobe.TimeStamp(find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID, 1,'first'));
lastStrobeKofikoTS = g_strctDAQParams.LastStrobe.TimeStamp(find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID, 1,'last'));
kofikoStrobeIDX = find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID);

plexontrialstartIDX = find(events.strobeNumber == startTrialID);
plexontrialstartTS = events.timeStamps(plexontrialstartIDX);
plexontrialendIDX = find(events.strobeNumber == endTrialID);
plexontrialendTS = events.timeStamps(plexontrialendIDX);

%haxxor correction right here
plexontrialtest=plexontrialendTS-plexontrialstartTS;
plexontrialstartTS(find(plexontrialtest<1.5))=[];
plexontrialendTS(find(plexontrialtest<1.5))=[];

kofikoStrobeAllTS = g_strctDAQParams.LastStrobe.TimeStamp(kofikoStrobeIDX);

if numel(kofikoStrobeIDX) ~= numel(plexonStrobeIDX)
    
    % sprintf('strobe ID mismatch, different number of sync timestamps detected')
    %sprintf('Will attempt to reconstruct sessions use sync timestamps from within recording sessions')
end

colorsInThisExperiment = zeros(1,3);
anglesInThisExperiment = zeros(1,1);
%%
%{
startRecordID = 32767;
stopRecordID = 32766;
%}

sessionStartPlexonIDX = find(events.strobeNumber == startRecordID);
sessionEndPlexonIDX = find(events.strobeNumber == stopRecordID);
sessionsStartPlexonTS = events.timeStamps(sessionStartPlexonIDX);
sessionEndPlexonTS = events.timeStamps(sessionEndPlexonIDX);
sessionsStartID = 1;
end

sessionsStartIDX = find(g_strctDAQParams.LastStrobe.Buffer == startRecordID);
sessionStartTS = g_strctDAQParams.LastStrobe.TimeStamp(sessionsStartIDX);% - firstStrobeKofikoTS;

sessionsEndIDX = find(g_strctDAQParams.LastStrobe.Buffer == stopRecordID);
sessionEndTS = g_strctDAQParams.LastStrobe.TimeStamp(sessionsEndIDX);% - firstStrobeKofikoTS;
numSessions = numel(sessionsStartIDX);

if isempty(sessionStartTS); sessionStartTS = g_strctDAQParams.LastStrobe.TimeStamp(1); end
if isempty(sessionEndTS); sessionEndTS = g_strctDAQParams.LastStrobe.TimeStamp(end); numSessions=1; end

%%
for iTrials = 1:size(ExperimentRecording,1)
    ExperimentRecording{iTrials, 2} = ExperimentRecording{iTrials, 1}.m_fImageFlipON_TS_Kofiko;
    trialIter(iTrials) = ExperimentRecording{iTrials, 1}.m_iTrialNumber;
end
allStartTS = [ExperimentRecording{:,2}];

%% if picking a particular subset of trials
targ_trials=[];
for tt=1:length(ExperimentRecording)
%     if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Dense Noise');
%         targ_trials=[targ_trials,tt];
%     end
    if isfield(ExperimentRecording{tt, 1}, 'm_iMkTurkTaskType') 
        if ~any(strfind(ExperimentRecording{tt, 1}.m_strctTrialOutcome.m_strResult, 'Aborted'))
            targ_trials=[targ_trials,tt]; 
        end
    end
end
ExperimentRecording=ExperimentRecording(targ_trials,:);
%% check plexon file
ntrials=length(ExperimentRecording);

topLevelIndex = [];
sessionIndex  = [];
sessionTags   = {};

cd ..
allPLXfiles = vertcat(dir([pwd,filesep, '*.plx']),dir([pwd,filesep, '*.pl2']));

thisSessionFile = plxFilePath;
if ~exist(plxFilePath) && isempty(allPLXfiles) && ~any(~arrayfun(@isempty,strfind({allPLXfiles(:).name}, [filenameP,'.plx']))) && ~any(~arrayfun(@isempty,strfind({allPLXfiles(:).name}, [filenameP,'.pl2'])))
    warning('could not find PLX file for this experiment')
else
  %  thisSessionFile = [pwd,'\',filenameP];
end


cd(ExperimentFolder)
unitsInFile = [];
[~,~,experimentFileExtension] = fileparts(thisSessionFile);
if ~strcmp(experimentFileExtension,'.plx') && ~strcmp(experimentFileExtension,'.pl2')
    thisSessionFile = [thisSessionFile, '.plx'];
end
%% align spiking data from plexon online sorting
if ~BehaviorOnly

numUnitsInSession = 0;
[tscounts, wfcounts, evcounts, contcounts] = plx_info([thisSessionFile], false);
    numUnits = find(sum(wfcounts,2));
    allNumUnits=[];
    topLevelIndex.m_iNumNeuronUnits = numUnits;
    
for channel=1:nChans
    numUnitsInSession = 0;
    for iUnit = 1:numel(numUnits)
        nameOfUnit = ['unit',num2str(numUnits(iUnit))];
        tempTS = 0;
        try
        [~, ~, tempTS, ~] = plx_waves_v([thisSessionFile], channel, iUnit);
%                [~, ~, tempTS, ~] = PL2Waves([thisSessionFile], channel, iUnit);
        catch
            fprintf('big oof \n')
        end
        if tempTS > 0
            numUnitsInSession = numUnitsInSession + 1
            unitsInFile = [unitsInFile, iUnit];
        [spikes(channel).(nameOfUnit).count, spikes(channel).(nameOfUnit).numWaves, spikes(channel).(nameOfUnit).timeStamps, spikes(channel).(nameOfUnit).Waves] = ...
                                                                        plx_waves_v([thisSessionFile], channel, iUnit);
%                 [spikes.(nameOfUnit).count, spikes.(nameOfUnit).numWaves, spikes.(nameOfUnit).timeStamps, spikes.(nameOfUnit).Waves] = ...
%                                                                                 PL2Waves([thisSessionFile], channel, iUnit);
        end
    end
    allNumUnits=[allNumUnits,numUnitsInSession];
end
numUnitsInSession=max(allNumUnits);
end
%% get ET data
load([ExperimentFolder, filenameP, '.mat'], 'g_strctEyeCalib');
load([ExperimentFolder, filenameP, '.mat'], 'g_strctStimulusServer');
%/{
% if loading ET data from plexon as well, uncomment this
[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
% ET_ad(1,:)=ET_ad(1,:)-median(ET_ad(1,:));
% ET_ad(2,:)=ET_ad(2,:)-median(ET_ad(2,:));
% ET_ad(1,:)=-(ET_ad(1,:)-ExperimentRecording{1, 1}.m_pt2iFixationSpot(1));
% ET_ad(2,:)=-(ET_ad(2,:)-ExperimentRecording{1, 1}.m_pt2iFixationSpot(2));
PlexET_times=[1:ET_n]/ET_adfreq;
%}
ET_ad = g_strctEyeCalib.EyeRaw.Buffer(:,1:2)';
ET_times = g_strctEyeCalib.EyeRaw.TimeStamp;
%% adjust ET data for fixation calibrations / recentering
Recalibs=unique(round([g_strctEyeCalib.CenterX.TimeStamp, g_strctEyeCalib.CenterY.TimeStamp, g_strctEyeCalib.GainX.TimeStamp, g_strctEyeCalib.GainY.TimeStamp]));
Recalibs=[Recalibs, ET_times(end)];
for rr=1:length(Recalibs)-1
    eyeCenterXID = find(Recalibs(rr) > g_strctEyeCalib.CenterX.TimeStamp,1,'last'); if isempty(eyeCenterXID); eyeCenterXID=1;end
    eyeCenterYID = find(Recalibs(rr) > g_strctEyeCalib.CenterY.TimeStamp,1,'last'); if isempty(eyeCenterYID);eyeCenterYID=1;end
    eyeGainXID = find(Recalibs(rr) > g_strctEyeCalib.GainX.TimeStamp,1,'last'); if isempty(eyeGainXID);eyeGainXID=1;end
    eyeGainYID = find(Recalibs(rr) > g_strctEyeCalib.GainY.TimeStamp,1,'last'); if isempty(eyeGainYID);eyeGainYID=1;end

    CalibEyeBufferIDX = find(g_strctEyeCalib.EyeRaw.TimeStamp > Recalibs(rr) & ...
        g_strctEyeCalib.EyeRaw.TimeStamp < Recalibs(rr+1));

    ET_ad(1,CalibEyeBufferIDX) = (ET_ad(1,CalibEyeBufferIDX) - g_strctEyeCalib.CenterX.Buffer(eyeCenterXID)) * g_strctEyeCalib.GainX.Buffer(eyeGainXID);% + (g_strctStimulusServer.m_aiScreenSize(3)/2)-ExperimentRecording{1, 1}.m_pt2iFixationSpot(1);
    ET_ad(2,CalibEyeBufferIDX) = (ET_ad(2,CalibEyeBufferIDX) - g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)) * g_strctEyeCalib.GainY.Buffer(eyeGainYID);% + (g_strctStimulusServer.m_aiScreenSize(4)/2)-ExperimentRecording{1, 1}.m_pt2iFixationSpot(2);
end
%{
% code to check ET traces
figure;
oL = length(ET_ad);
downsampled = interp1(1:oL, ET_ad(1,:), linspace(1,oL,round(oL/10)));
for tt=1000:1000:10e6
    plot(downsampled(1,tt:tt+999)'-mean(downsampled(1,tt:tt+999)')); hold on
    plot(diff(downsampled(1,tt:tt+999))); hold off
    pause
end
%}
%% SACCADE DETECTION (from detect_saccades_v2
if ~skip_saccades
fprintf('Detecting saccades\n');
eye_smooth=4; eye_dt=1;
sm_avg_eyepos = ET_ad';
sm_avg_eyepos(:,1) = smooth(sm_avg_eyepos(:,1),eye_smooth);
sm_avg_eyepos(:,2) = smooth(sm_avg_eyepos(:,2),eye_smooth);
eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]./eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]./eye_dt;
all_eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
% %
%parameters
et_params.eye_fs = 120;%ET_adfreq;
sac_thresh = 0.7; %threshold eye speed
peri_thresh = 1; %threshold eye speed for defining saccade boundary inds
min_isi = 0.1; max_isi = Inf; %min/max inter-saccade intervals

%find local maxima of eye speed signal exceeding sac_thresh
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

% %
%find times when speed signal crossed above and below the peri-saccade
%threshold
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
end

%% get rid of double-peaks
% isis = [Inf; diff(sac_start_inds)]/et_params.eye_fs;
% bad_isis = (isis < min_isi | isis > max_isi);
% bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
% saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
isi_cutoff=10;
nsacs=length(saccade_inds);
ssi=1; cur_ssi=2;
while ssi<nsacs
    cur_isi=saccade_inds(cur_ssi)-saccade_inds(cur_ssi-1);
    if cur_isi<isi_cutoff & all_eye_speed(saccade_inds(cur_ssi)) <= all_eye_speed(saccade_inds(cur_ssi-1));
        saccade_inds(cur_ssi)=[];
        sac_start_inds(cur_ssi) = []; 
        sac_stop_inds(cur_ssi) = [];
        ssi=ssi+1;
    elseif cur_isi<isi_cutoff & all_eye_speed(saccade_inds(cur_ssi)) >= all_eye_speed(saccade_inds(cur_ssi-1));
        saccade_inds(cur_ssi-1)=[];
        sac_start_inds(cur_ssi-1) = []; 
        sac_stop_inds(cur_ssi-1) = [];
        ssi=ssi+1;
    else
        cur_ssi=cur_ssi+1;
        ssi=ssi+1;
    end
end

% next_isi = [isis(2:end); Inf];
saccade_times = ET_times(saccade_inds); %saccade peak times
% sac_start_times = ET_times(sac_start_inds); %saccade start times
% sac_stop_times = ET_times(sac_stop_inds); %saccade end times
%sac_durs = sac_stop_times - sac_start_times; %saccade durations 
%sac_peakvel = all_eye_speed(saccade_inds); %peak eye vels

sprintf('detected %d saccades: %d per second \n',length(saccade_inds), length(saccade_inds)./(length(ET_ad)./100))
end
% %% %
% figure;
% for sact=600:length(saccade_inds)
% plot([-400:1:1200],ET_ad(:,saccade_inds(sact)-400:saccade_inds(sact)+1200)'.*1.5); hold on
% plot([-400:1:1200],all_eye_speed(saccade_inds(sact)-400:saccade_inds(sact)+1200),'k');
% vline(0); hold off
% ylim([-60 60])
% ylabel('Arcmin'); xlabel('Time (ms)')
% pause
% end

%%
sprintf('done with preprocessing \n')
%% get spikesorted timestamps
if ~BehaviorOnly

for chan=1:nChans
[LFP_adfreq, LFP_n, LFP_ts, LFP_fn, LFP_ad(:,chan)] = plx_ad_v(thisSessionFile, ['FP' num2str(chan,'%02.f')]);
end
LFP_ad=LFP_ad';
LFP_times=[1:LFP_n]/LFP_adfreq;
%%
% Aud: you can safely skip this part
if useofflinesorting==1
%     load([matFilePath, filenameP 'SortedSpikeTimes.mat']);
%     allNumUnits=ones(length(allchan_spktimes));
%     unitsInFile=ones(length(allchan_spktimes));
%     % extend this for multiple SUs per channel + MU
    load([strExperimentPath, filenameP '_Clusters.mat']);
    for channel=1:nChans
        curunits=unique(Clusters{channel}.spike_clusts);
        allNumUnits(channel)=length(find(curunits>1));
    end
end
%%
exptDataP = []; 
for channel=1:nChans    
    for iUnit = 1:allNumUnits(channel)
        nameOfUnit = ['unit',num2str(iUnit)];
        if useofflinesorting==1
            exptDataP(channel).(nameOfUnit) = Clusters{channel}.times(find(Clusters{channel}.spike_clusts==iUnit+1));
            %exptDataMUA(channel).(nameOfUnit) = allchan_spktimes{2,channel};
        else
            exptDataP(channel).(nameOfUnit) = spikes(channel).(nameOfUnit).timeStamps;
        end
    end
    if useofflinesorting==1
        if ismember(channel,iso_SUs)
        exptDataMUA{channel} = Clusters{channel}.times(find(Clusters{channel}.spike_clusts==1));
        else
        exptDataMUA{channel} = Clusters{channel}.times(find(Clusters{channel}.spike_clusts>=1));
        end
    end
end

% iSessions = sessionsStartIDX(iSessions);
plexonSyncStrobesInThisSessionTS = events.timeStamps;
plexonSyncStrobesInThisSessionStrobeID =   events.strobeNumber;
plexonSyncStrobesInThisSessionTS =  plexonSyncStrobesInThisSessionTS(plexonSyncStrobesInThisSessionStrobeID == syncStrobeID);
kofikoSyncStrobesInThisSessionTS = g_strctDAQParams.LastStrobe.TimeStamp;

kofikoSyncStrobesInThisSessionStrobeID = g_strctDAQParams.LastStrobe.Buffer;
kofikoSyncStrobesInThisSessionTS = kofikoSyncStrobesInThisSessionTS(kofikoSyncStrobesInThisSessionStrobeID == syncStrobeID);
if isempty(kofikoSyncStrobesInThisSessionTS) || isempty(plexonSyncStrobesInThisSessionTS)
    sprintf('could not process recording %s, session %i, timestamp missing.', filenameP, iSessions)
end

end

%%
thisExptEyeRaw=g_strctEyeCalib.EyeRaw.Buffer;
thisExptEyeTimes=g_strctEyeCalib.EyeRaw.TimeStamp;
N_recalibs=g_strctEyeCalib.CenterX.TimeStamp;

% figure;
% for cals=1:N_recalibs-1
%     curfixinds=find(thisExptEyeTimes>g_strctEyeCalib.CenterX.TimeStamp(cals) & thisExptEyeTimes<g_strctEyeCalib.CenterX.TimeStamp(cals+1));
%     y=thisExptEyeRaw(curfixinds,1);
%     x=thisExptEyeTimes(curfixinds);
%     X = [ones(length(x),1) x'];
%     b = X\y;
%     plot(x,y,'b'); hold on
%     plot(x,sum(y*b',2), 'r'); hold off
% %    plot(x,y-sum(y*b',2), 'r--'); hold off
% figtitle(num2str(b));
%     pause
% end
%% 
%[ExperimentRecording{:,12}] == 1);
trialsInThisSession{iSessions} = find(allStartTS> sessionStartTS(1) & ...
    allStartTS< sessionEndTS(end));
trialIter = 1;
for iFixationCheck = 1:ntrials; %trialsInThisSession{iSessions};%
    if isfield(ExperimentRecording{iFixationCheck, 1},'m_bMonkeyFixated') && ExperimentRecording{iFixationCheck, 1}.m_bMonkeyFixated
        ExperimentRecording{iFixationCheck,7} = 1;
        sessionIndex(trialIter).m_bFixated = 1;
        sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = [];
        % append eye trace information to this trial
        thisTrialFlipon = ExperimentRecording{iFixationCheck,2} ;
            eyeCenterXID = find(thisTrialFlipon > g_strctEyeCalib.CenterX.TimeStamp,1,'last'); if isempty(eyeCenterXID); eyeCenterXID=1;end
            eyeCenterYID = find(thisTrialFlipon > g_strctEyeCalib.CenterY.TimeStamp,1,'last'); if isempty(eyeCenterYID);eyeCenterYID=1;end
            eyeGainXID = find(thisTrialFlipon > g_strctEyeCalib.GainX.TimeStamp,1,'last'); if isempty(eyeGainXID);eyeGainXID=1;end
            eyeGainYID = find(thisTrialFlipon > g_strctEyeCalib.GainY.TimeStamp,1,'last'); if isempty(eyeGainYID);eyeGainYID=1;end
            
            thisTrialEyeBufferIDX = find(g_strctEyeCalib.EyeRaw.TimeStamp > thisTrialFlipon + g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(1) & ...
                g_strctEyeCalib.EyeRaw.TimeStamp < thisTrialFlipon + g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(2) );
            
            thisTrialRawEyeData = g_strctEyeCalib.EyeRaw.Buffer(thisTrialEyeBufferIDX,:);
            thisTrialRawEyeDatatimes = g_strctEyeCalib.EyeRaw.TimeStamp(1,thisTrialEyeBufferIDX);
           
            fEyeXPix = (thisTrialRawEyeData(:, 1) - g_strctEyeCalib.CenterX.Buffer(eyeCenterXID)) * g_strctEyeCalib.GainX.Buffer(eyeGainXID) + (g_strctStimulusServer.m_aiScreenSize(3)/2);
            fEyeYPix = (thisTrialRawEyeData(:, 2) - g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)) * g_strctEyeCalib.GainY.Buffer(eyeGainYID) + (g_strctStimulusServer.m_aiScreenSize(4)/2);
            ExperimentRecording{(iFixationCheck),1}.m_afEyeXPositionScreenCoordinates = fEyeXPix;
            ExperimentRecording{(iFixationCheck),1}.m_afEyeYPositionScreenCoordinates = fEyeYPix;
            ExperimentRecording{iFixationCheck,1}.m_afEyePositiontimes = thisTrialRawEyeDatatimes;
        ExperimentRecording{iFixationCheck,1}.ETthisTrialRawEyeData = thisTrialRawEyeData;
        ExperimentRecording{iFixationCheck,1}.ETCenter = [g_strctEyeCalib.CenterX.Buffer(eyeCenterXID),g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)];
        ExperimentRecording{iFixationCheck,1}.ScreenCenter = [(g_strctStimulusServer.m_aiScreenSize(3)/2),(g_strctStimulusServer.m_aiScreenSize(4)/2)];
            sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = fEyeXPix;
            sessionIndex(trialIter).m_afEyeYPositionScreenCoordinates = fEyeYPix;
    elseif ~isfield(ExperimentRecording{(iFixationCheck), 1},'m_bMonkeyFixated') || ~ExperimentRecording{(iFixationCheck), 1}.m_bMonkeyFixated
        ExperimentRecording{(iFixationCheck),7} = 0;
        sessionIndex(trialIter).m_bFixated = 0;
        sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = [];        % append eye trace information to this trial
        thisTrialFlipon = ExperimentRecording{(iFixationCheck),2} ;
        eyeCenterXID = find(thisTrialFlipon > g_strctEyeCalib.CenterX.TimeStamp,1,'last'); if isempty(eyeCenterXID); eyeCenterXID=1;end
        eyeCenterYID = find(thisTrialFlipon > g_strctEyeCalib.CenterY.TimeStamp,1,'last'); if isempty(eyeCenterYID);eyeCenterYID=1;end
        eyeGainXID = find(thisTrialFlipon > g_strctEyeCalib.GainX.TimeStamp,1,'last'); if isempty(eyeGainXID);eyeGainXID=1;end
        eyeGainYID = find(thisTrialFlipon > g_strctEyeCalib.GainY.TimeStamp,1,'last'); if isempty(eyeGainYID);eyeGainYID=1;end

        thisTrialEyeBufferIDX = g_strctEyeCalib.EyeRaw.TimeStamp > thisTrialFlipon + g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(1) & ...
            g_strctEyeCalib.EyeRaw.TimeStamp < thisTrialFlipon + g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(2) ;

        thisTrialRawEyeData = g_strctEyeCalib.EyeRaw.Buffer(thisTrialEyeBufferIDX,:);
        thisTrialRawEyeDatatimes = g_strctEyeCalib.EyeRaw.TimeStamp(1,thisTrialEyeBufferIDX);

        fEyeXPix = (thisTrialRawEyeData(:, 1) - g_strctEyeCalib.CenterX.Buffer(eyeCenterXID)) * g_strctEyeCalib.GainX.Buffer(eyeGainXID) + (g_strctStimulusServer.m_aiScreenSize(3)/2);
        fEyeYPix = (thisTrialRawEyeData(:, 2) - g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)) * g_strctEyeCalib.GainY.Buffer(eyeGainYID) + (g_strctStimulusServer.m_aiScreenSize(4)/2);
        ExperimentRecording{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates = fEyeXPix;
        ExperimentRecording{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates = fEyeYPix;
        ExperimentRecording{iFixationCheck,1}.m_afEyePositiontimes = thisTrialRawEyeDatatimes;
        ExperimentRecording{iFixationCheck,1}.ETthisTrialRawEyeData = thisTrialRawEyeData;
        ExperimentRecording{iFixationCheck,1}.ETCenter = [g_strctEyeCalib.CenterX.Buffer(eyeCenterXID),g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)];
        ExperimentRecording{iFixationCheck,1}.ScreenCenter = [(g_strctStimulusServer.m_aiScreenSize(3)/2),(g_strctStimulusServer.m_aiScreenSize(4)/2)];
        sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = fEyeXPix;
        sessionIndex(trialIter).m_afEyeYPositionScreenCoordinates = fEyeYPix;

    else
        spintf('big oof: no ET data')
    end
    trialIter = trialIter + 1;
%     figure(1);
%     plot(ExperimentRecording{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
%     plot(ExperimentRecording{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
% 
%     figure(3);
%     plot(ExperimentRecording{iFixationCheck,5}')
%     pause
end

%%
if ~BehaviorOnly

trialIter = 1;
for iTrials = 1:ntrials; %trialsInThisSession{iSessions}

    if iTrials==1
        ExperimentRecording{iTrials,1}.isCatchTrial=0;
    else
        if ExperimentRecording{iTrials-1,1}.m_iMkTurkTaskType==ExperimentRecording{iTrials,1}.m_iMkTurkTaskType;
            ExperimentRecording{iTrials,1}.isCatchTrial=0;
        else
            ExperimentRecording{iTrials,1}.isCatchTrial=1; 
        end
    end
    
    spikesInThisTrial = [];
    sessionIndex(trialIter).m_iGlobalTrialIndex = iTrials;
%    [~,trialSyncStrobeID] = min(abs(ExperimentRecording{iTrials, 2}- kofikoSyncStrobesInThisSessionTS));
    trialSyncStrobeID = find(kofikoSyncStrobesInThisSessionTS < ExperimentRecording{iTrials, 2},1,'last');

    if trialSyncStrobeID>length(kofikoSyncStrobesInThisSessionTS)
        sprintf(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(kofikoSyncStrobesInThisSessionTS))]);
        trialSyncStrobeID=length(kofikoSyncStrobesInThisSessionTS);
    end
    if trialSyncStrobeID>length(plexonSyncStrobesInThisSessionTS)
        sprintf(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(plexonSyncStrobesInThisSessionTS))]);
        trialSyncStrobeID=length(plexonSyncStrobesInThisSessionTS);
    end
    
    kofikoSyncTime = kofikoSyncStrobesInThisSessionTS(trialSyncStrobeID);
    plexonSyncTime = plexonSyncStrobesInThisSessionTS(trialSyncStrobeID);
    if ~isfield( ExperimentRecording{iTrials, 1},'m_aiStimColor') || isempty(ExperimentRecording{iTrials, 1}.m_aiStimColor)
        ExperimentRecording{iTrials, 1}.m_aiStimColor = [NaN, NaN, NaN];
    end
%    ExperimentRecording{iTrials, 3} = ExperimentRecording{iTrials, 1}.m_aiStimColor;
%        ExperimentRecording{iTrials, 4} = ExperimentRecording{iTrials, 1}.m_fRotationAngle;
%    ExperimentRecording{iTrials, 6} = strcmpi(ExperimentRecording{iTrials, 1}.m_strTrialType,'moving bar');
%        sessionIndex(trialIter).m_fRotationAngle = ExperimentRecording{iTrials, 1}.m_fRotationAngle;
%    sessionIndex(trialIter).m_bIsMovingBarTrial =  ExperimentRecording{iTrials, 6};
    ExperimentRecording{iTrials, 3} = ExperimentRecording{iTrials, 1}.m_strctTrialOutcome.m_strResult;

    plexontrialID = find(plexontrialstartTS - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) > 0, 1, 'first');
    
    ExperimentRecording{iTrials,1}.plexonSyncTime = plexonSyncTime;
    ExperimentRecording{iTrials,1}.kofikoSyncTime = kofikoSyncTime;
    ExperimentRecording{iTrials,1}.plexontrialstart = plexontrialstartTS(plexontrialID);
    ExperimentRecording{iTrials,1}.plexontrialend = plexontrialendTS(plexontrialID);

    ExperimentRecording{iTrials,1}.PreCueFixMS  = ExperimentRecording{iTrials,1}.m_strctPreCuePeriod.m_fPreCueFixationPeriodMS;
    ExperimentRecording{iTrials,1}.CueFixMS  = ExperimentRecording{iTrials,1}.m_strctCuePeriod.m_fCuePeriodMS;
    ExperimentRecording{iTrials,1}.MemoryMS  = ExperimentRecording{iTrials,1}.m_strctMemoryPeriod.m_fMemoryPeriodMS;
    ExperimentRecording{iTrials,1}.ChoiceMemoryMS  = ExperimentRecording{iTrials, 1}.m_strctMemoryPeriod.m_fMemoryChoicePeriodMS;
    ExperimentRecording{iTrials,1}.trialend = plexontrialendTS(plexontrialID)-plexontrialstartTS(plexontrialID);

%    ncols=size(ExperimentRecording,2);
    plexonLFPDataAlignedToThisTrial = LFP_ad(:, LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,4} = plexonLFPDataAlignedToThisTrial;
   
    plexonETDataAlignedToThisTrial = PlexET_ad(:, PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,5} = plexonETDataAlignedToThisTrial;

% for aligning EMs to Kofiko timestamps
if ~skip_saccades
    EMsInThisTrial =  saccade_times(saccade_times >= ExperimentRecording{iTrials, 2}-g_strctStatistics.preTrialWindow & ...
        saccade_times <= ExperimentRecording{iTrials, 2}+g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,6} = EMsInThisTrial;
end 

% for aligning EMs to plexon timestamps    
%     plexonEMsAlignedToThisTrial =  saccade_times - plexonSyncTime ;
%     EMsInThisTrial = plexonEMsAlignedToThisTrial(plexonEMsAlignedToThisTrial - ...
%         (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
%         plexonEMsAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
%         - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);

    for channel=1:nChans
        for iUnit = 1:allNumUnits(channel)
            plexonDataAlignedToThisTrial = [];
            nameOfUnit = ['unit',num2str(iUnit)];
%            nameOfUnit = ['unit',num2str(unitsInFile(iUnit))];            
            plexonDataAlignedToThisTrial =  exptDataP(channel).(nameOfUnit) - plexonSyncTime ;
            spikesInThisTrial.(nameOfUnit) = plexonDataAlignedToThisTrial(plexonDataAlignedToThisTrial - ...
                (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
                plexonDataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
                - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);
        end
        plexonMUADataAlignedToThisTrial =  exptDataMUA{channel} - plexonSyncTime ;
        MUAspikesInThisTrial = plexonMUADataAlignedToThisTrial(plexonMUADataAlignedToThisTrial - ...
            (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            plexonMUADataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
            - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);

        ExperimentRecording{iTrials,8+3*(channel-1)} = spikesInThisTrial;
        ExperimentRecording{iTrials,9+3*(channel-1)} = length(spikesInThisTrial.unit1);
        sessionIndex(trialIter).m_afSpikesInThisTrial = spikesInThisTrial;
        ExperimentRecording{iTrials,10+3*(channel-1)} = MUAspikesInThisTrial;
        
    end
    
    trialIter = trialIter + 1;
end 

% if ~isempty(experimentIndex)
%     save([strExperimentPath,filesep,'experimentIndex.mat'],'experimentIndex')
%     globalSessionTags = experimentIndex{:,3};
%     save([strExperimentPath,filesep,'experimentTags.mat'],'globalSessionTags')
% end

end
cd(currentDirectory);
fprintf('done \n')
% %
save([strExperimentPath,filesep,'ExperimentRecording_AllTrials.mat'],'ExperimentRecording','-v7.3')

% %% adjust ExperimentRecording if not done before
% ExperimentRecording{1,1}.isCatchTrial=0;
% for ii=2:length(ExperimentRecording)
%     if ExperimentRecording{ii-1,1}.m_iMkTurkTaskType==ExperimentRecording{ii,1}.m_iMkTurkTaskType
%         ExperimentRecording{ii,1}.isCatchTrial=0;
%     else
%         ExperimentRecording{ii,1}.isCatchTrial=1; 
%     end
% end
end
%%
%% get indices
% correct and incorrect
Corr_inds=[]; Incorr_inds=[]; Ttasktype=[]; Toutcome=[]; Tcatchtrial=[];
for ii=1:length(ExperimentRecording)
    Toutcome{ii}=ExperimentRecording{ii, 1}.m_strctTrialOutcome.m_strResult;
    if strcmp(ExperimentRecording{ii, 1}.m_strctTrialOutcome.m_strResult, 'Incorrect')
        Incorr_inds = [Incorr_inds,ii];
    elseif strcmp(ExperimentRecording{ii, 1}.m_strctTrialOutcome.m_strResult, 'Correct')
        Corr_inds = [Corr_inds,ii];
    end
    Ttasktype(ii) = ExperimentRecording{ii, 1}.m_iMkTurkTaskType;
    Tcatchtrial(ii) = ExperimentRecording{ii, 1}.isCatchTrial;
end


%%
targ_inds=find(Ttasktype>=1);
[Toutcomelist, ~, Toutcomeinds]= unique(Toutcome(targ_inds));
Toutcome_abort=length(find(Toutcomeinds==1))./length(targ_inds);
Toutcome_correct=length(find(Toutcomeinds==2))./length(targ_inds);
Toutcome_incorrect=length(find(Toutcomeinds==3))./length(targ_inds);
Toutcome_timeout=length(find(Toutcomeinds==4))./length(targ_inds);

% % %
% bar([Toutcome_correct, Toutcome_incorrect, Toutcome_timeout, Toutcome_abort])
% ax=gca; ax.XTickLabels={'Correct','Incorrect','Timeout','Aborted'}
% title(filenameP)


%% show STAs and LFPs for Correct and Incorrect trials
Spks_corr=[]; Spks_incorr=[]; Spks_corr_Raligned=[]; Spks_incorr_Raligned=[]; 
STA_corr=[]; STA_incorr=[]; STA_corr_Raligned=[]; STA_incorr_Raligned=[]; 
LFP_corr = []; LFP_incorr = []; LFP_corr_Raligned=[]; LFP_incorr_Raligned=[];

for trl=1:length(Corr_inds)
    for cc=1:24
    Spks_corr(cc,:,trl) = histcounts(ExperimentRecording{Corr_inds(trl), 8+(cc-1)*3}.unit1,[g_strctStatistics.preTrialWindow:.02:g_strctStatistics.postTrialWindow]);
    cur_trialend = round(ExperimentRecording{Corr_inds(trl),1}.trialend*1000)/1000;
    Spks_corr_Raligned(cc,:,trl) = histcounts(ExperimentRecording{Corr_inds(trl), 8+(cc-1)*3}.unit1,[cur_trialend-1.5:.02:cur_trialend+1]);
    end
    LFP_corr(:,:,trl) = ExperimentRecording{Corr_inds(trl),4};
    LFP_corr_Raligned(:,:,trl) = ExperimentRecording{Corr_inds(trl),4}(:,1000*cur_trialend-1500:1000*cur_trialend+1000);
end
STA_corr = squeeze(mean(Spks_corr,3)); 
STA_corr_Raligned = squeeze(mean(Spks_corr_Raligned,3)); 
LFPavg_corr = squeeze(mean(LFP_corr,3)); 
LFPavg_corr_Raligned = squeeze(mean(LFP_corr_Raligned,3)); 

for trl=1:length(Incorr_inds)
    for cc=1:24
    Spks_incorr(cc,:,trl) = histcounts(ExperimentRecording{Incorr_inds(trl), 8+(cc-1)*3}.unit1,[g_strctStatistics.preTrialWindow:.02:g_strctStatistics.postTrialWindow]);
    cur_trialend = round(ExperimentRecording{Incorr_inds(trl),1}.trialend*1000)/1000;
    Spks_incorr_Raligned(cc,:,trl) = histcounts(ExperimentRecording{Incorr_inds(trl), 8+(cc-1)*3}.unit1,[cur_trialend-1.5:.02:cur_trialend+1]);
    end
    LFP_incorr(:,:,trl) = ExperimentRecording{Incorr_inds(trl),4};
    LFP_incorr_Raligned(:,:,trl) = ExperimentRecording{Incorr_inds(trl),4}(:,1000*cur_trialend-1500:1000*cur_trialend+1000);
end
STA_incorr = squeeze(mean(Spks_incorr,3)); 
STA_incorr_Raligned = squeeze(mean(Spks_incorr_Raligned,3)); 
LFPavg_incorr = squeeze(mean(LFP_incorr,3)); 
LFPavg_incorr_Raligned = squeeze(mean(LFP_incorr_Raligned,3)); 

%% plot STAs of correct and incorrect trials
for cc=[2:24];
    figure(1)
    subplot(2,2,1); plot([-500:20:6480],STA_corr(cc,:)*50); title('Correct Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(0, 'k', {'Cue onset'})
    vline(1000, 'k',{'Cue offset'})
    vline(1500, 'k',{'Choice onset'});
    subplot(2,2,3); plot([-500:20:6480],STA_incorr(cc,:)*50); title('Incorrect Trials')
    vline(0, 'k', {'Cue onset'})
    vline(1000, 'k',{'Cue offset'})
    vline(1500, 'k',{'Choice onset'});
    
    subplot(2,2,2); plot([-1480:20:1000],STA_corr_Raligned(cc,:)*50); title('Correct Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(-500, 'k',{'Saccade'})
    vline(0, 'k',{'Choice locked'})
    subplot(2,2,4); plot([-1480:20:1000],STA_incorr_Raligned(cc,:)*50); title('Incorrect Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(-500, 'k',{'Saccade'})
    vline(0, 'k',{'Choice locked'})
    
    figtitle(['Probe '  num2str(cc) ' Unit 1'])
    pause
end

%% plot LFPs of correct and incorrect trials
figure(2)
subplot(2,2,1); imagesc(LFPavg_corr); title('Correct Trials'); xlabel('time (ms)'); ylabel('Channel')
subplot(2,2,3); imagesc(LFPavg_incorr); title('Incorrect Trials')
    
subplot(2,2,2); imagesc(LFPavg_corr_Raligned); title('Correct Trials'); xlabel('time (ms)'); ylabel('Channel')
subplot(2,2,4); imagesc(LFPavg_incorr_Raligned); title('Incorrect Trials')
figtitle(['LFP'])


%% Now compare arbitrary sets of indices
%% show STAs and LFPs for any subset of trials
% Corr_inds; Incorr_inds; Ttasktype; Toutcome; Tcatchtrial;
var1_inds = intersect(Corr_inds, find(Ttasktype==1));
var2_inds = intersect(Corr_inds, find(Ttasktype==2));

Spks_var1=[]; Spks_var2=[]; Spks_var1_Raligned=[]; Spks_var2_Raligned=[]; 
STA_var1=[]; STA_var2=[]; STA_var1_Raligned=[]; STA_var2_Raligned=[]; 
LFP_var1 = []; LFP_var2 = []; LFP_var1_Raligned=[]; LFP_var2_Raligned=[];

for trl=1:length(var1_inds)
    for cc=1:24
    Spks_var1(cc,:,trl) = histcounts(ExperimentRecording{var1_inds(trl), 8+(cc-1)*3}.unit1,[g_strctStatistics.preTrialWindow:.02:g_strctStatistics.postTrialWindow]);
    cur_trialend = round(ExperimentRecording{var1_inds(trl),1}.trialend*1000)/1000;
    Spks_var1_Raligned(cc,:,trl) = histcounts(ExperimentRecording{var1_inds(trl), 8+(cc-1)*3}.unit1,[cur_trialend-1.5:.02:cur_trialend+1]);
    end
    LFP_var1(:,:,trl) = ExperimentRecording{var1_inds(trl),4};
    LFP_var1_Raligned(:,:,trl) = ExperimentRecording{var1_inds(trl),4}(:,1000*cur_trialend-1500:1000*cur_trialend+1000);
end
STA_var1 = squeeze(mean(Spks_var1,3)); 
STA_var1_Raligned = squeeze(mean(Spks_var1_Raligned,3)); 
LFPavg_var1 = squeeze(mean(LFP_var1,3)); 
LFPavg_var1_Raligned = squeeze(mean(LFP_var1_Raligned,3)); 

for trl=1:length(var2_inds)
    for cc=1:24
    Spks_var2(cc,:,trl) = histcounts(ExperimentRecording{var2_inds(trl), 8+(cc-1)*3}.unit1,[g_strctStatistics.preTrialWindow:.02:g_strctStatistics.postTrialWindow]);
    cur_trialend = round(ExperimentRecording{var2_inds(trl),1}.trialend*1000)/1000;
    Spks_var2_Raligned(cc,:,trl) = histcounts(ExperimentRecording{var2_inds(trl), 8+(cc-1)*3}.unit1,[cur_trialend-1.5:.02:cur_trialend+1]);
    end
    LFP_var2(:,:,trl) = ExperimentRecording{var2_inds(trl),4};
    LFP_var2_Raligned(:,:,trl) = ExperimentRecording{var2_inds(trl),4}(:,1000*cur_trialend-1500:1000*cur_trialend+1000);
end
STA_var2 = squeeze(mean(Spks_var2,3)); 
STA_var2_Raligned = squeeze(mean(Spks_var2_Raligned,3)); 
LFPavg_var2 = squeeze(mean(LFP_var2,3)); 
LFPavg_var2_Raligned = squeeze(mean(LFP_var2_Raligned,3)); 

%% plot STAs of var1 and var2 trials
for cc=[1:24];% [1 3 6 9 11 13 15 20 21 22];
    figure(1)
    subplot(2,2,1); plot([-500:20:6480],STA_var1(cc,:)*50); title('Color Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(0, 'k', {'Cue onset'})
    vline(1000, 'k',{'Cue offset'})
    vline(1500, 'k',{'Choice onset'});
    subplot(2,2,3); plot([-500:20:6480],STA_var2(cc,:)*50); title('Shape Trials')
    vline(0, 'k', {'Cue onset'})
    vline(1000, 'k',{'Cue offset'})
    vline(1500, 'k',{'Choice onset'});
    
    subplot(2,2,2); plot([-1480:20:1000],STA_var1_Raligned(cc,:)*50); title('Color Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(-500, 'k',{'Saccade'})
    vline(0, 'k',{'Choice locked'})
    subplot(2,2,4); plot([-1480:20:1000],STA_var2_Raligned(cc,:)*50); title('Shape Trials'); xlabel('time (ms)'); ylabel('Firing rate')
    vline(-500, 'k',{'Saccade'})
    vline(0, 'k',{'Choice locked'})
    
    figtitle(['Probe '  num2str(cc) ' Unit 1'])
    pause
end

%% plot STAs of var1 and var2 trials on top of each other
for cc=[5:15 20:22];
    figure(1)
    subplot(2,1,1); plot([-500:20:6480],STA_var1(cc,:)*50); hold on
                    plot([-500:20:6480],STA_var2(cc,:)*50); 
    xlabel('time (ms)'); ylabel('Firing rate');
    vline(0, 'k', {'Cue onset'})
    vline(1000, 'k',{'Cue offset'})
    vline(1500, 'k',{'Choice onset'});
    legend({'Color trials', 'shape trials'}); hold off
    
    subplot(2,1,2); plot([-1480:20:1000],STA_var1_Raligned(cc,:)*50); hold on
                    plot([-1480:20:1000],STA_var2_Raligned(cc,:)*50);
    xlabel('time (ms)'); ylabel('Firing rate')
    vline(-500, 'k',{'Saccade'})
    vline(0, 'k',{'Choice locked'})
    legend({'Color trials', 'shape trials'}); hold off

    
    figtitle(['Probe '  num2str(cc) ' Unit 1'])
    pause
end
%% plot raster plot for var1 and var2
SpksRaster_var1=[]; SpksRaster_var2=[];
for trl=1:length(var1_inds)
    for cc=1:24
    SpksRaster_var1{cc}{trl} = ExperimentRecording{var1_inds(trl), 8+(cc-1)*3}.unit1;
    end
end
for trl=1:length(var2_inds)
    for cc=1:24
    SpksRaster_var2{cc}{trl} = ExperimentRecording{var2_inds(trl), 8+(cc-1)*3}.unit1;
    end
end
%%
for cc=1:24;
    subplot(2,1,1); plotSpikeRaster(SpksRaster_var1{cc});
    subplot(2,1,2); plotSpikeRaster(SpksRaster_var2{cc});
    pause
end
%% plot LFPs of var1 and var2 trials
figure(4)
subplot(2,2,1); imagesc(LFPavg_var1); title('Color Trials'); xlabel('time (ms)'); ylabel('Channel')
subplot(2,2,3); imagesc(LFPavg_var2); title('Shape Trials')
    
subplot(2,2,2); imagesc(LFPavg_var1_Raligned); title('Response Aligned'); xlabel('time (ms)'); ylabel('Color Trials')
subplot(2,2,4); imagesc(LFPavg_var2_Raligned); ylabel('Shape Trials')
figtitle(['LFP'])

%% time-frequency plot
for chan=[15, 22];
    figure;
newtimef(LFP_var1(chan,:,:), 7000,[-500 3000], 1000);
figtitle(['Color Trial time-freq for channel ' num2str(chan)])
pause
end
%%
for chan=[15, 22];%[2:6, 9, 12, 14, 16, 19];
    figure;
newtimef(LFP_var2(chan,:,:), 7000,[-500 3000], 1000);
figtitle(['Shape Trial time-freq for channel ' num2str(chan)])
pause
end
%% 
%% for reading performance across days (loop over days defined as eeday)
ToutcomeMat(1,eeday)=Toutcome_correct;
ToutcomeMat(2,eeday)=Toutcome_incorrect;
ToutcomeMat(3,eeday)=Toutcome_abort;
ToutcomeMat(4,eeday)=Toutcome_timeout;

%%
figure;
plot(ToutcomeMat');
legend({'Correct','Incorrect','Abort','Timeout'})
vline(4,'k', 'Delay implemented')
vline(11,'k', 'Feedback implemented')
vline(19,'k', 'Back from Quarantine')

%%
for eeday=1:size(ToutcomeMat,2)
TPcorrectMat(1,eeday)=ToutcomeMat(1,eeday)./[ToutcomeMat(1,eeday)+ToutcomeMat(2,eeday)];
end
%%
figure;
 plot(TPcorrectMat)
 title('% Correct of completed trials')
vline(4,'k', 'Delay implemented')
vline(11,'k', 'Feedback implemented')
