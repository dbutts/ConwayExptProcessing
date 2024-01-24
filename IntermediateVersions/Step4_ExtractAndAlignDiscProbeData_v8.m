%%
clear all

%% Extract Bevil's data

filenameP = '220722_144441_Jacomo'; ET_Eyelink=1; iso_SUs=[];

metadata.filename_plexon = filenameP;
filedate=filenameP(1:6);

useofflinesorting = 1;
skipLFP=0; 
%LFPchans=[1:24];
LFPchans=[1:24, 33,40,46,47,52,53,54,59,65,67,71,81,83,89,90,95,98,102,103,109,112,131,138,139,145,146,152,158, 161:256];
nChans=24; spks.nchans = nChans; %channel=17;

sessionsToProcess = [];
ExperimentFolder = '/media/felix/Internal_1/Data/BevilColor/';

%plxFilePath = [ExperimentFolder filenameP '.pl2'];
%plxFilePath = ['/media/felix/Seagate Portable Drive/Data_Color/' filenameP '.pl2'];
plxFilePath = ['/home/felix/NTlab_dataserver4/Conway/raw/' filenameP '.pl2'];

matFilePath = [ExperimentFolder filenameP '/'];
configFilePath = [ExperimentFolder filenameP '.mat'];

% KSstitched=0; ksFilePath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_laminar/']; arraylabel ='laminar';
% KSstitched=1; ksFilePath = ['/media/felix/Internal_1/Data/BevilColor/'...
% filenameP '/kilosorting_Utah_24chs_stitch_Ethan/']; arraylabel ='Ethan_Utah';
KSstitched=1; ksFilePath = [ExperimentFolder ...
filenameP '/kilosorting_stitch/']; arraylabel ='full';

sessionTimeOffsets = [];                                     
strExperimentPath = [matFilePath '/Analysis/'];
if ~exist(strExperimentPath,'dir');
    mkdir(strExperimentPath);
end
filenameK = matFilePath;
load(configFilePath)
%%
global strctColorValues unitsInFile PlottingVars g_strctStatistics ExperimentRecording topLevelIndex 

%persistent plexonDataAlignedToThisTrial

experimentIndex = {};

g_strctStatistics.preTrialWindow = -0.100;
g_strctStatistics.postTrialWindow = 0.300; %0.4

currentDirectory = pwd;
cd(ExperimentFolder);
allMatFiles = dir('*.mat');
[~, indices] = sort(vertcat(allMatFiles(:).datenum));
allMatFiles = allMatFiles(indices);

allPLXfiles = dir([pwd,filesep, '*.plx']);
g_strctStatistics.trackedXYZCoordinates = [];
g_strctStatistics.m_aiXYZCoorindateColorMatchingTable = [];
g_strctStatistics.bPlotRFFieldWithinStimRectangle = 1;
g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod = [-.050,4];

g_strctStatistics.numBins = 200;
g_strctStatistics.positionBinGrid = [12,20];
g_strctStatistics.bPlotRFFieldWithinStimRectangle = 1;
g_strctStatistics.RFMappingTrialIntegrationTime(:,1) = -.05:.01:.34 ;
g_strctStatistics.RFMappingTrialIntegrationTime(:,2) = -.04:.01:.35 ;
g_strctStatistics.m_bPlotRFFieldByColor = 1;

syncStrobeID = 32757;
eventChannelNumber = 257;
startRecordID = 32767;
stopRecordID = 32766;

load([ExperimentFolder filenameP '.mat'], 'g_astrctAllParadigms')
ProbeCols=g_astrctAllParadigms{1, 1}.m_cPresetProbeColors;
ProbeColIDs=g_astrctAllParadigms{1, 1}.m_cPresetProbeColorIDS;

%%
if useofflinesorting==1

    if KSstitched==1
        load([ksFilePath 'KS_stitched.mat'])
    else
        spk_times = readNPY([ksFilePath 'spike_times_seconds.npy']);
        spk_clusters = readNPY([ksFilePath 'spike_clusters.npy']);
        spk_info = tdfread([ksFilePath 'cluster_info.tsv']);
    end
%    cd('/media/felix/Internal_1/Data/BevilColor/220209_163451_Jacomo/kilosorting_laminar')
%     spk_times = readNPY('spike_times_seconds.npy');
%     spk_clusters = readNPY('spike_clusters.npy');
%     spk_labels = tdfread('cluster_KSLabel.tsv');

   spk_clustIDs = unique(spk_clusters); nclusts=length(spk_clustIDs);
%    spk_clustIDs = spk_info.cluster_id; nclusts=length(spk_clustIDs);
     
    spk_labels_SU=[]; spk_labels_MU=[];
    for cc=1:nclusts
        if strcmp(deblank(spk_info.group(cc,:)), 'good')
            spk_labels_SU = [spk_labels_SU,cc];
        elseif strcmp(spk_info.group(cc,:), 'noise')
            % % do nothing about noise
        else
             spk_labels_MU = [spk_labels_MU,cc];
        end
    end
    bad_chans_SU = find(spk_info.n_spikes(spk_labels_SU)<2000); spk_labels_SU(bad_chans_SU)=[];
    bad_chans_MU = find(spk_info.n_spikes(spk_labels_MU)<2000); spk_labels_MU(bad_chans_MU)=[];
    spk_ID_SU = (spk_clustIDs(spk_labels_SU));
    try
        spk_rating_SU = (spk_info.Rating(spk_labels_SU));
        spk_rating_MU = (spk_info.Rating(spk_labels_MU));
    catch
        spk_rating_SU = (spk_info.rating(spk_labels_SU));
        spk_rating_MU = (spk_info.rating(spk_labels_MU));
    end
    spk_channels_SU = spk_info.ch(spk_labels_SU);
    spk_ID_MU = (spk_clustIDs(spk_labels_MU));
    spk_channels_MU = spk_info.ch(spk_labels_MU);
    nSU=length(spk_ID_SU);
    nMU=length(spk_ID_MU);
end
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

try
    [events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plxFilePath, eventChannelNumber);
%    [events] = PL2EventTs(plxFilePath, eventChannelNumber);
catch
    [~,sessionName] = fileparts(plxFilePath);
    sprintf('corrupt or missing information from Plexon file, experiment %s', sessionName)
    return;
end
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

firstStrobePlexonTS = events.timeStamps(find(events.strobeNumber(events.strobeNumber == syncStrobeID),1));
plexonStrobeIDX = find(events.strobeNumber == syncStrobeID);
lastStrobePlexonTS = events.timeStamps(find(events.strobeNumber(events.strobeNumber == syncStrobeID),1,'last'));
plexonStrobeAllTS = events.timeStamps(plexonStrobeIDX);
firstStrobeKofikoTS = g_strctDAQParams.LastStrobe.TimeStamp(find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID, 1,'first'));
lastStrobeKofikoTS = g_strctDAQParams.LastStrobe.TimeStamp(find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID, 1,'last'));
kofikoStrobeIDX = find(g_strctDAQParams.LastStrobe.Buffer == syncStrobeID);

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
sessionsStartIDX = find(g_strctDAQParams.LastStrobe.Buffer == startRecordID);
sessionStartTS = g_strctDAQParams.LastStrobe.TimeStamp(sessionsStartIDX);% - firstStrobeKofikoTS;

sessionsEndIDX = find(g_strctDAQParams.LastStrobe.Buffer == stopRecordID);
sessionEndTS = g_strctDAQParams.LastStrobe.TimeStamp(sessionsEndIDX);% - firstStrobeKofikoTS;
numSessions = numel(sessionsStartIDX);

if strfind(filenameP, '190409_134012_Jacomo')
numSessions=8;
end
sessionStartPlexonIDX = find(events.strobeNumber == startRecordID);
sessionEndPlexonIDX = find(events.strobeNumber == stopRecordID);
sessionsStartPlexonTS = events.timeStamps(sessionStartPlexonIDX);
sessionEndPlexonTS = events.timeStamps(sessionEndPlexonIDX);
sessionsStartID = 1;
%%
for iTrials = 1:size(ExperimentRecording,1)
    ExperimentRecording{iTrials, 2} = ExperimentRecording{iTrials, 1}.m_fImageFlipON_TS_Kofiko;
    trialIter(iTrials) = ExperimentRecording{iTrials, 1}.m_iTrialNumber;
    ExperimentRecording{iTrials, 1}.SessionID = find(ExperimentRecording{iTrials, 2}>sessionStartTS,1,'last');
end
%% if picking a particular subset of trials
targ_trials=[];
for tt=1:length(ExperimentRecording)
    if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Disc Probe') && ...
                ExperimentRecording{tt, 1}.m_bMonkeyFixated==1;
    targ_trials=[targ_trials,tt];
    end
end
%%
ExperimentRecording=ExperimentRecording(targ_trials,:);


%% if going by trial
% for iTrials = 1:size(ExperimentRecording,1)
%     ExperimentRecording{iTrials, 2} = 
% end

ntrials=length(ExperimentRecording);

topLevelIndex = [];
sessionIndex  = [];
sessionTags   = {};

cd ..
allPLXfiles = vertcat(dir([pwd,filesep, '*.plx']),dir([pwd,filesep, '*.pl2']));

%disp(sprintf('No individual PLX file found for %s session %i, using original plexon file', filenameP ,sessionsToProcess(itSessions)))

%allPLXfiles = dir([pwd,filesep,'*.plx']);
% [allDirs, allFiles] = subdir(pwd);
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
%% align spiking data

numUnitsInSession = 0;
[tscounts, wfcounts, evcounts, contcounts] = plx_info([thisSessionFile], false);
%    for channel=1:32
%        numUnits = find(sum(wfcounts(2:end,:),2));
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
%% get ET data
load([ExperimentFolder, filenameP, '.mat'], 'g_strctEyeCalib');
load([ExperimentFolder, filenameP, '.mat'], 'g_strctStimulusServer');

%%
% figure; 
% [~, ~, ~, ~, PlexET_trigs] = plx_ad_v(thisSessionFile, 'AI02');
% plot(PlexET_trigs(1e5:2e5))
%%
%/{
[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI05');
[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI06');
[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
%ET_ad(1,:)=ET_ad(1,:)-median(ET_ad(1,:));
%ET_ad(2,:)=ET_ad(2,:)-median(ET_ad(2,:));
PlexET_times=[1:ET_n]/ET_adfreq;
%}
%%
ET_ad = g_strctEyeCalib.EyeRaw.Buffer(:,1:2)';
ET_times = g_strctEyeCalib.EyeRaw.TimeStamp;
%%
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
%%
%{
figure;
oL = length(ET_ad);
downsampled = interp1(1:oL, ET_ad(1,:), linspace(1,oL,round(oL/10)));
for tt=1000:1000:10e6
    plot(downsampled(1,tt:tt+999)'-mean(downsampled(1,tt:tt+999)')); hold on
    plot(diff(downsampled(1,tt:tt+999))); hold off
    pause
end
%}
%%
ET_times_up = ET_times(1):.001:ET_times(end);
ET_ad_up(1,:) = interp1(ET_times, ET_ad(1,:), ET_times_up, 'spline');
ET_ad_up(2,:) = interp1(ET_times, ET_ad(2,:), ET_times_up, 'spline');

ET_ad_up_s(1,:) = smooth(ET_ad_up(1,:), 'sgolay', 3);
ET_ad_up_s(2,:) = smooth(ET_ad_up(2,:), 'sgolay', 3);

%% SACCADE DETECTION (from detect_saccades_v2
fprintf('Detecting saccades\n');

sac_thresh = 4; %threshold eye speed
peri_thresh = 1; %threshold eye speed for defining saccade boundary inds
min_isi = 0.1; max_isi = Inf; %min/max inter-saccade intervals

eye_smooth=3; eye_dt=1;
sm_avg_eyepos = ET_ad_up_s';

%for amplitude correction
amp_cutoff=6; %arcmins
win=[3 8];

%sm_avg_eyepos(:,1) = smooth(sm_avg_eyepos(:,1),eye_smooth);
%sm_avg_eyepos(:,2) = smooth(sm_avg_eyepos(:,2),eye_smooth);

eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]./eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]./eye_dt;
all_eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
ET_n=length(all_eye_speed);
% %
%parameters
et_params.eye_fs = 1000;%ET_adfreq;
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
sac_stop_inds(end)=ET_n;

%% identify artifacts (saccade stayed above peri_thresh for too long)
bad_sac_inds=find(diff(sac_start_inds)==0)+1;
sac_start_inds(bad_sac_inds)=[];
sac_stop_inds(bad_sac_inds)=[];
saccade_inds(bad_sac_inds)=[];

bad_sac_inds=find(saccade_inds>[ET_n-win(2)-1]);
sac_start_inds(bad_sac_inds)=[];
sac_stop_inds(bad_sac_inds)=[];
saccade_inds(bad_sac_inds)=[];
%{
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
nsacs=length(saccade_inds);
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
%}
%%
sac_stop_inds(end)=ET_n;
% next_isi = [isis(2:end); Inf];
saccade_times = ET_times_up(saccade_inds); %saccade peak times
sac_start_times = ET_times_up(sac_start_inds); %saccade start times
sac_stop_times = ET_times_up(sac_stop_inds); %saccade end times
%sac_durs = sac_stop_times - sac_start_times; %saccade durations 
%sac_peakvel = all_eye_speed(saccade_inds); %peak eye vels

sprintf('detected %d saccades: %d per second \n',length(saccade_inds), length(saccade_inds)./(length(ET_ad_up_s)./et_params.eye_fs))
%% % 
%{
figure;
winpre=200; winpost=800;
for sact=3500:length(saccade_inds)
    subplot(2,1,1)
    ETwintimes = ET_times_up(saccade_inds(sact)-winpre:saccade_inds(sact)+winpost)-ET_times_up(saccade_inds(sact));
    plot(ETwintimes,ET_ad_up(:,saccade_inds(sact)-winpre:saccade_inds(sact)+winpost)'); hold on
    plot(ETwintimes,all_eye_speed(saccade_inds(sact)-winpre:saccade_inds(sact)+winpost),'k');
vline(0); hold off; %axis tight
ylim([-60 60]); %xlim([-.5 1])
ylabel('Arcmin'); xlabel('Time (ms)'); 

subplot(2,1,2)
plot(ETwintimes,sm_avg_eyepos(saccade_inds(sact)-winpre:saccade_inds(sact)+winpost,:)); hold on
plot(ETwintimes,all_eye_speed(saccade_inds(sact)-winpre:saccade_inds(sact)+winpost),'k');
vline(0); hold off; %axis tight
ylim([-60 60]); %xlim([-.5 1])
ylabel('Arcmin'); xlabel('Time (ms)');

pause
end
%}
%%
sprintf('done with saccades \n')
%%
LFPcc=1;
if ~skipLFP
    for chan=LFPchans
        [LFP_adfreq, LFP_n, LFP_ts, LFP_fn, LFP_ad(:,LFPcc)] = plx_ad_v(thisSessionFile, ['FP' num2str(chan,'%03.f')]);
        LFPcc=LFPcc+1;
    end
    LFP_ad=LFP_ad';
    LFP_times=[1:LFP_n]/LFP_adfreq;
end
%% previous location of reading in kilosort outputs
%%
exptDataP = []; exptDataP2=[]; exptDataMUA=[];

if useofflinesorting==1
    for iUnit=1:length(spk_labels_SU)
        exptDataP(iUnit).spkID = double(spk_ID_SU(iUnit));
        exptDataP(iUnit).spkCh = spk_channels_SU(iUnit);
        exptDataP(iUnit).rating = spk_rating_SU(iUnit);
        exptDataP(iUnit).unit1 = spk_times(find(spk_clusters==spk_ID_SU(iUnit)));
    end

    for iUnit=1:length(spk_labels_MU)
        exptDataMUA(iUnit).spkID = double(spk_ID_MU(iUnit));
        exptDataMUA(iUnit).spkCh = spk_channels_MU(iUnit);
        exptDataMUA(iUnit).rating = spk_rating_MU(iUnit);
        exptDataMUA(iUnit).unit1 = spk_times(find(spk_clusters==spk_ID_MU(iUnit)));
    end

else
    for channel=1:nChans    
        for iUnit = 1:allNumUnits(channel)
            nameOfUnit = ['unit',num2str(iUnit)];
            if useofflinesorting==1
                %{
                exptDataP(channel).(nameOfUnit) = Clusters{channel}.times(find(Clusters{channel}.spike_clusts==iUnit+1));
                %exptDataP(channel).(nameOfUnit) = allchan_spktimes{1,channel};
                %exptDataMUA(channel).(nameOfUnit) = allchan_spktimes{2,channel};
                %}
            else
                try
                exptDataP(channel).(nameOfUnit) = spikes(channel).(nameOfUnit).timeStamps;
                catch
                exptDataP(channel).(nameOfUnit) = 0;   
                end
            end
        end

    end

end
iSessions = 1;%
allStartTS = [ExperimentRecording{:,2}];
%%
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
%{
for iFixationCheck = 1:ntrials; %trialsInThisSession{iSessions};%
figure(1); 
plot(ExperimentRecording{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
plot(ExperimentRecording{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
figure(2); 
plot(ExperimentRecording{iFixationCheck,1}.ETthisTrialRawEyeData);
pause
end
%}
%%
trialIter = 1;
disp(['found ' num2str(ntrials) ' trials'])
for iTrials = 1:ntrials;%trialsInThisSession{iSessions}
%%
    spikesInThisTrial = [];
    sessionIndex(trialIter).m_iGlobalTrialIndex = iTrials;
    [~,trialSyncStrobeID] = min(abs(ExperimentRecording{iTrials, 2} - kofikoSyncStrobesInThisSessionTS));

    if trialSyncStrobeID>length(kofikoSyncStrobesInThisSessionTS)
        error(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(kofikoSyncStrobesInThisSessionTS))]);
        %trialSyncStrobeID=length(kofikoSyncStrobesInThisSessionTS);
    end
    if trialSyncStrobeID>length(plexonSyncStrobesInThisSessionTS)
        error(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(plexonSyncStrobesInThisSessionTS))]);
        %trialSyncStrobeID=length(plexonSyncStrobesInThisSessionTS);
    end
    
    kofikoSyncTime = kofikoSyncStrobesInThisSessionTS(trialSyncStrobeID);
    plexonSyncTime = plexonSyncStrobesInThisSessionTS(trialSyncStrobeID);
% check if offsets need work?
%     if (ExperimentRecording{iTrials, 2} - kofikoSyncTime)>0
%         kofikoSyncTime = kofikoSyncStrobesInThisSessionTS(trialSyncStrobeID+1);
%         plexonSyncTime = plexonSyncStrobesInThisSessionTS(trialSyncStrobeID+1);
%     end
    ExperimentRecording{iTrials,1}.kofikoSyncTime = kofikoSyncTime; 
    ExperimentRecording{iTrials,1}.PlexonSyncTime = plexonSyncTime;
    ExperimentRecording{iTrials,1}.PlexonOnsetTime = plexonSyncTime + (ExperimentRecording{iTrials, 2} - kofikoSyncTime);
%    [ExperimentRecording{iTrials,1}.PlexonOnsetTime, (ExperimentRecording{iTrials, 2} - kofikoSyncTime)]
%%

    if ~isfield( ExperimentRecording{iTrials, 1},'m_aiStimColor') || isempty(ExperimentRecording{iTrials, 1}.m_aiStimColor)
        ExperimentRecording{iTrials, 1}.m_aiStimColor = [NaN, NaN, NaN];
    end
%    ExperimentRecording{iTrials, 3} = ExperimentRecording{iTrials, 1}.m_aiStimColor;
%        ExperimentRecording{iTrials, 4} = ExperimentRecording{iTrials, 1}.m_fRotationAngle;
%    ExperimentRecording{iTrials, 6} = strcmpi(ExperimentRecording{iTrials, 1}.m_strTrialType,'moving bar');
%        sessionIndex(trialIter).m_fRotationAngle = ExperimentRecording{iTrials, 1}.m_fRotationAngle;
%    sessionIndex(trialIter).m_bIsMovingBarTrial =  ExperimentRecording{iTrials, 6};
    ExperimentRecording{iTrials, 3} = ExperimentRecording{iTrials, 1}.m_strTrialType;

    sessionIndex(trialIter).m_aiStimulusColor = ExperimentRecording{iTrials, 1}.m_aiStimColor;
    sessionIndex(trialIter).m_strTrialType = ExperimentRecording{iTrials, 1}.m_strTrialType;


    if useofflinesorting==1
        for iUnit=1:nSU
            plexonDataAlignedToThisTrial = [];         
            plexonDataAlignedToThisTrial =  exptDataP(iUnit).unit1 - plexonSyncTime;
            spikesInThisTrial.unit1 = plexonDataAlignedToThisTrial(plexonDataAlignedToThisTrial - ...
                (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
                plexonDataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow) ...
                - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);
            spikesInThisTrial.spkID = exptDataP(iUnit).spkID;
            spikesInThisTrial.rating = exptDataP(iUnit).rating;
            ExperimentRecording{iTrials,10+3*(iUnit-1)} = spikesInThisTrial;
            ExperimentRecording{iTrials,11+3*(iUnit-1)} = length(spikesInThisTrial.unit1);
            %sessionIndex(trialIter).m_afSpikesInThisTrial = spikesInThisTrial;
        end

        for iUnit=1:nMU
            plexonMUADataAlignedToThisTrial=[];
            plexonMUADataAlignedToThisTrial =  exptDataMUA(iUnit).unit1 - plexonSyncTime ;
            MUAspikesInThisTrial.unit1 = plexonMUADataAlignedToThisTrial(plexonMUADataAlignedToThisTrial - ...
                (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
                plexonMUADataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
                - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);

            MUAspikesInThisTrial.spkID = exptDataMUA(iUnit).spkID;
            ExperimentRecording{iTrials,10+3*(iUnit-1+nSU)} = MUAspikesInThisTrial.unit1;
            ExperimentRecording{iTrials,11+3*(iUnit-1+nSU)} = length(MUAspikesInThisTrial.unit1);
        end

    else
        for channel=find(allNumUnits);%1:nChans
            for iUnit = 1:allNumUnits(channel)
                plexonDataAlignedToThisTrial = [];
                nameOfUnit = ['unit',num2str(iUnit)];
    %            nameOfUnit = ['unit',num2str(unitsInFile(iUnit))];            
                plexonDataAlignedToThisTrial =  exptDataP(channel).(nameOfUnit) - plexonSyncTime ;
                spikesInThisTrial.(nameOfUnit) = plexonDataAlignedToThisTrial(plexonDataAlignedToThisTrial - ...
                    (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
                    plexonDataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
                    - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);;
            end
    %         plexonMUADataAlignedToThisTrial =  exptDataMUA{channel} - plexonSyncTime ;
    %         MUAspikesInThisTrial = plexonMUADataAlignedToThisTrial(plexonMUADataAlignedToThisTrial - ...
    %             (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
    %             plexonMUADataAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
    %             - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);
    
    
            ExperimentRecording{iTrials,10+3*(channel-1)} = spikesInThisTrial;
            ExperimentRecording{iTrials,11+3*(channel-1)} = length(spikesInThisTrial.unit1);
            sessionIndex(trialIter).m_afSpikesInThisTrial = spikesInThisTrial;
    %        ExperimentRecording{iTrials,12+3*(channel-1)} = MUAspikesInThisTrial;
            
        end
    end
        
    if ~skipLFP
%    ncols=size(ExperimentRecording,2);
    plexonLFPDataAlignedToThisTrial = LFP_ad(:, LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,4} = plexonLFPDataAlignedToThisTrial;
    end


    curColID=find(sum(ProbeCols==ExperimentRecording{iTrials, 1}.DiscprobeColor,2)==3,1,'first');
    ExperimentRecording{iTrials,8}=ProbeColIDs(curColID,1:3);
    ExperimentRecording{iTrials,9}=ProbeColIDs(curColID,[6 4 5]);

% % for aligning EMs to plexon timestamps    
%     plexonEMsAlignedToThisTrial =  saccade_times - plexonSyncTime ;
%     EMsInThisTrial = plexonEMsAlignedToThisTrial(plexonEMsAlignedToThisTrial - ...
%         (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
%         plexonEMsAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
%         - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);

% for aligning EMs to Kofiko timestamps
    EMsInThisTrial.saccades =  saccade_times(saccade_times >= ExperimentRecording{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
        saccade_times <= ExperimentRecording{iTrials, 2}+g_strctStatistics.postTrialWindow);
    EMsInThisTrial.saccade_start =  sac_start_times(saccade_times >= ExperimentRecording{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
        sac_start_times <= ExperimentRecording{iTrials, 2}+g_strctStatistics.postTrialWindow);
    EMsInThisTrial.saccade_stop =  sac_stop_times(saccade_times >= ExperimentRecording{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
        sac_stop_times <= ExperimentRecording{iTrials, 2}+g_strctStatistics.postTrialWindow);
    EMsInThisTrial.ET_times =  ET_times(ET_times >= ExperimentRecording{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
        ET_times <= ExperimentRecording{iTrials, 2}+g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,6} = EMsInThisTrial;

%     plexonETDataAlignedToThisTrial.ET_trace = PlexET_ad(:, PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
%             PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
%     plexonETDataAlignedToThisTrial.ET_times = PlexET_times(:, PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
%             PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    plexonETDataAlignedToThisTrial.ET_trace = PlexET_ad(:, PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    plexonETDataAlignedToThisTrial.ET_times = PlexET_times(:, PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
            PlexET_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
    ExperimentRecording{iTrials,7} = plexonETDataAlignedToThisTrial;

    trialIter = trialIter + 1;
    
    if mod(iTrials,50)==0
        fprintf('finished trial %d of %d \n',iTrials, ntrials)
    end
    
end 

cd(currentDirectory);
disp('done. saving...')

save([strExperimentPath,filesep,'FullDiscProbeExperimentRecordingwLFP_ks' num2str(useofflinesorting) '_' arraylabel '_v08.mat'],'ExperimentRecording','-v7.3')
disp('saved!')

%%
trlColIDs_polar=cell2mat(ExperimentRecording(:,8));
trlspikes=cell2mat(ExperimentRecording(:,[11:3:size(ExperimentRecording,2)]));

spike_ID=double([spk_ID_SU; spk_ID_MU]);
probe_ID=[spk_channels_SU; spk_channels_MU];
%% check acorss saturation
cols_satur=ProbeCols(1:45,:); 
cols_satur2=reshape(cols_satur, 3,15,3);

saturIDs=[0.33 0.66 1];
trlinds_satur=[];
for theta=1:15
    for satur = 1:3
        trlinds_satur{theta,satur}=intersect(find(trlColIDs_polar(:,1)==theta), find(trlColIDs_polar(:,2)==saturIDs(satur)));
    end
end
%%
cur_tuning_satur=[];

probetuningfig=figure;

for cc=1:size(trlspikes,2)
    for theta=1:15
        for satur = 1:3
            cur_tuning_satur(theta,satur) = mean(trlspikes(trlinds_satur{theta,satur},cc));
        end
    end
    subplot(2,1,1);
    imagesc(cur_tuning_satur'); title('Mean response')
    h=colorbar;title(h, 'Firing rate'); 
    set(gca, 'XtickLabels',[]); set(gca, 'YtickLabels',[]); ylabel('Saturation');

    subplot(2,1,2);
    image(cols_satur2./255); title('Corresponding screen colors')
    xlabel('Theta'); set(gca, 'YtickLabels',[]); ylabel('Color'); colorbar;

    figtitle(['Rec ' filedate  ' Spike ID: ' num2str(spike_ID(cc)) ' - Probe ' num2str(probe_ID(cc))])
    saveas(probetuningfig,[strExperimentPath 'Spike ID ' num2str(spike_ID(cc)) ' - Probe ' num2str(probe_ID(cc)) '_DiscProbe_Saturation.png'])
    %pause
    all_cells_disctuning(cc).probe = probe_ID(cc);
    all_cells_disctuning(cc).unitID = spike_ID(cc);
    all_cells_disctuning(cc).tuning_satur=cur_tuning_satur;

end

%% check acorss elevation
elevIDs = [-.8 -.5 -.2 .2 .5 .8];
trlinds_elev=[];
for theta=1:15
    for elev = 1:6
        trlinds_elev{theta,elev}=intersect(find(trlColIDs_polar(:,1)==theta), find(trlColIDs_polar(:,3)==elevIDs(elev)));
    end
end
cols_elev=ProbeCols(46:135,:); 
cols_elev2=reshape(cols_elev, 6,15,3);
%%
cur_tuning_elev=[];

probetuningfig=figure;

for cc=1:size(trlspikes,2)
    for theta=1:15
        for elev = 1:6
            cur_tuning_elev(theta,elev) = mean(trlspikes(trlinds_elev{theta,elev},cc));
        end
    end
    subplot(2,1,1);
    imagesc(cur_tuning_elev'); h=colorbar;title(h, 'firing rate'); set(gca, 'XtickLabels',[]);
    ylabel('luminance'); set(gca, 'YtickLabels',[]); 

    subplot(2,1,2);
    image(cols_elev2./255); colorbar;
    xlabel('theta'); ylabel('Color'); set(gca, 'YtickLabels',[]); 
    figtitle(['Rec ' filedate ' Spike ID: ' num2str(spike_ID(cc)) ' - Probe ' num2str(probe_ID(cc))])
    saveas(probetuningfig,[strExperimentPath 'Spike ID ' num2str(spike_ID(cc)) ' - Probe ' num2str(probe_ID(cc)) '_DiscProbe_Elevation.png'])

    all_cells_disctuning(cc).tuning_elev=cur_tuning_elev;
    %pause
end
%%
cur_filename=[strExperimentPath 'Jocamo_' filenameP(1:6) '_' arraylabel '_discprobetuning.mat'];
save(cur_filename, 'all_cells_disctuning');
%%
disp('DONE')