%%
clear all
%%
%filenames = {'210123_130251_Jacomo'}; iso_SUs=[5,6,8,9,10,13,16,17,21];
%filenames = {'210130_162547_Jacomo'}; iso_SUs=[2, 4:8, 11, 13, 14, 15, 19, 20, 24];
filenames = {'210202_101158_Jacomo'}; iso_SUs=[8, 10, 17, 18, 19, 20, 23, 24];

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
eye_smooth=3; eye_dt=1;
sm_avg_eyepos = ET_ad';
sm_avg_eyepos(:,1) = smooth(sm_avg_eyepos(:,1),eye_smooth);
sm_avg_eyepos(:,2) = smooth(sm_avg_eyepos(:,2),eye_smooth);
eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]./eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]./eye_dt;
all_eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
ET_n=length(all_eye_speed);
% %
%parameters
et_params.eye_fs = 120;%ET_adfreq;
sac_thresh = 5; %threshold eye speed
peri_thresh = 1; %threshold eye speed for defining saccade boundary inds
min_isi = 0.1; max_isi = Inf; %min/max inter-saccade intervals
isi_cutoff=min_isi*et_params.eye_fs;

%find local maxima of eye speed signal exceeding sac_thresh
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);
saccade_inds(saccade_inds<isi_cutoff)=[];
saccade_inds(saccade_inds>(ET_n-isi_cutoff))=[];
saccade_inds = [1; saccade_inds; ET_n];
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
sac_start_inds(1) = 1;
sac_stop_inds(end)=size(ET_times,2);

%% get rid of double-peaks
% isis = [Inf; diff(sac_start_inds)]/et_params.eye_fs;
% bad_isis = (isis < min_isi | isis > max_isi);
% bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
% saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
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

%% correct for amplitude
amp_cutoff=20; %arcmins
nsacs=length(saccade_inds);
win=[3 8];
ssi=1; cur_ssi=3;
while ssi<nsacs-2
    cur_preampX=median(sm_avg_eyepos(saccade_inds(cur_ssi)-win(2):saccade_inds(cur_ssi)-win(1),1));
    cur_postampX=median(sm_avg_eyepos(saccade_inds(cur_ssi)+win(1):saccade_inds(cur_ssi)+win(2),1));
    cur_preampY=median(sm_avg_eyepos(saccade_inds(cur_ssi)-win(2):saccade_inds(cur_ssi)-win(1),2));
    cur_postampY=median(sm_avg_eyepos(saccade_inds(cur_ssi)+win(1):saccade_inds(cur_ssi)+win(2),2));
    if [abs(cur_postampX-cur_preampX)+abs(cur_postampY-cur_preampY)] < amp_cutoff;
        saccade_inds(cur_ssi)=[];
        sac_start_inds(cur_ssi) = []; 
        sac_stop_inds(cur_ssi) = [];
        ssi=ssi+1;
    else
        cur_ssi=cur_ssi+1;
        ssi=ssi+1;
    end
end

% next_isi = [isis(2:end); Inf];
saccade_times = ET_times(saccade_inds); %saccade peak times
sac_start_times = ET_times(sac_start_inds); %saccade start times
sac_stop_times = ET_times(sac_stop_inds); %saccade end times
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

%end
