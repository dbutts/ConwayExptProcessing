%%
clear all

% DAN NEEDS PATH -- all dependencies subdirectories
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/CSD_functions/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilosort2/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Plexon-Matlab Offline Files SDK/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')

%% Extract Bevil's data

%filenameP = '220616_130047_Jacomo'; ET_Eyelink=2; spk_offset=0; % set to 2 for 0616 to 0709 for eyelink with monoc-only traces
filenameP = '230510_141725_Jacomo'; spk_offset=0; ET_Eyelink=3; % set to 3 for any dDPI tracked files
iso_SUs=[]; 

useofflinesorting = 1; % set to 1 in order to use kilosort outputs, otherwise 0

skipLFP=0; 
LFPchans=[1:24];
LFP_chanmap=LFPchans;
nChans=24; %for CSD only now

sessionsToProcess = [];
ExperimentFolder = '/Users/dbutts/Data/Conway/';  % EXPERIMENT FOLDER
plxFilePath = [ExperimentFolder filenameP '.pl2'];
matFilePath = [ExperimentFolder filenameP '/']; % NEED SLASHES THE CORRECT DIRECTION
configFilePath = [ExperimentFolder filenameP '.mat'];

% ALSO FOLLOWING LINE BACKWARDS /
KSstitched=0; ksFilePath = [ExperimentFolder filenameP '/kilosorting_laminar/']; arraylabel ='lam';
%KSstitched=0; ksFilePath = [ExperimentFolder filenameP '\kilosorting_laminar\']; arraylabel ='lam';

sessionTimeOffsets = [];
%strExperimentPath = [matFilePath '\Analysis\'];
strExperimentPath = [matFilePath '/Analysis/'];
if ~exist(strExperimentPath,'dir')
	mkdir(strExperimentPath);
end
filenameK = matFilePath;
load(configFilePath)

%% test for misalignment
%/{
[SPKC_adfreq, SPKC_n, SPKC_ts, SPKC_fn, SPKC_ad_test] = plx_ad_v(plxFilePath, 'SPKC003');
[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad_test] = plx_ad_v(plxFilePath, 'AI07');
[LFP_adfreq, LFP_n, LFP_ts, LFP_fn, LFP_ad_test] = plx_ad_v(plxFilePath, ['FP' num2str(1,'%03.f')]);
ET_n./LFP_n
ET_n-LFP_n
if abs(ET_n-LFP_n)>100
	warning('Danger - Plexon timing bug dectected')
else
	disp('Plexon timing check passed')
end
%}

%%
global strctColorValues unitsInFile PlottingVars g_strctStatistics ExperimentRecording topLevelIndex 

%persistent plexonDataAlignedToThisTrial
experimentIndex = {};

g_strctStatistics.preTrialWindow = 0;
g_strctStatistics.postTrialWindow = 4; %0.4

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

%%
if useofflinesorting==1

	if KSstitched==1
		load([ksFilePath 'KS_stitched.mat'])
	else
		spk_times = readNPY([ksFilePath 'spike_times_seconds.npy']) +spk_offset;
		spk_clusters = readNPY([ksFilePath 'spike_clusters.npy']);
		spk_info = tdfread([ksFilePath 'cluster_info.tsv']);
	end
	
	spk_clustIDs = unique(spk_clusters); nclusts=length(spk_clustIDs);
	spk_labels_SU=[]; spk_labels_MU=[];
	for cc=1:nclusts
		if strcmp(deblank(spk_info.group(cc,:)), 'good')
			spk_labels_SU = [spk_labels_SU,cc];
		elseif strcmp(spk_info.group(cc,:), 'noise')
			% do nothing about noise
		else
			spk_labels_MU = [spk_labels_MU,cc];
		end
	end
	bad_chans_SU = find(spk_info.n_spikes(spk_labels_SU)<2000); spk_labels_SU(bad_chans_SU)=[];
	bad_chans_MU = find(spk_info.n_spikes(spk_labels_MU)<2000); spk_labels_MU(bad_chans_MU)=[];
	spk_ID_SU = (spk_clustIDs(spk_labels_SU));
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
	% [events] = PL2EventTs(plxFilePath, eventChannelNumber);
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
	sprintf('strobe ID mismatch, different number of sync timestamps detected')
	%sprintf('Will attempt to reconstruct sessions use sync timestamps from within recording sessions')
end

colorsInThisExperiment = zeros(1,3);
anglesInThisExperiment = zeros(1,1);

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
	if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Dual Stim');
	% if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Fivedot');
		targ_trials=[targ_trials,tt];
	end
end

%%
ExperimentRecording=ExperimentRecording(targ_trials,:);
%% ensure trials are in correct order
ntrials=length(ExperimentRecording);

for tt=1:ntrials
	trialstart_raw(tt)=ExperimentRecording{tt, 2}; 
end
[~,order] = sort(trialstart_raw);
ExperimentRecording = ExperimentRecording(order, :);  

%%
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
			% [~, ~, tempTS, ~] = PL2Waves([thisSessionFile], channel, iUnit);
		catch
			fprintf('big oof \n')
		end
		if tempTS > 0
			numUnitsInSession = numUnitsInSession + 1
			unitsInFile = [unitsInFile, iUnit];
			[spikes(channel).(nameOfUnit).count, spikes(channel).(nameOfUnit).numWaves, spikes(channel).(nameOfUnit).timeStamps, spikes(channel).(nameOfUnit).Waves] = ...
				plx_waves_v([thisSessionFile], channel, iUnit);
			% [spikes.(nameOfUnit).count, spikes.(nameOfUnit).numWaves, spikes.(nameOfUnit).timeStamps, spikes.(nameOfUnit).Waves] = ...
			%		PL2Waves([thisSessionFile], channel, iUnit);
		end
	end
	allNumUnits=[allNumUnits,numUnitsInSession];
end

numUnitsInSession=max(allNumUnits);
%% get ET data
load([ExperimentFolder, filenameP, '.mat'], 'g_strctEyeCalib');
load([ExperimentFolder, filenameP, '.mat'], 'g_strctStimulusServer');

%/{
if ET_Eyelink == 1
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI05');
	[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI06');
	[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);
	PlexET_ad_calib(3,:) = (PlexET_ad_calib(3,:)-median(PlexET_ad_calib(3,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(4,:) = (PlexET_ad_calib(4,:)-median(PlexET_ad_calib(4,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 2
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 3
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI03');
	[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI04');
	[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(3,:) = (PlexET_ad_calib(3,:)-median(PlexET_ad_calib(3,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(4,:) = (PlexET_ad_calib(4,:)-median(PlexET_ad_calib(4,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 0
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad(1,:)=PlexET_ad(1,:)-median(PlexET_ad(1,:));
	PlexET_ad(2,:)=PlexET_ad(2,:)-median(PlexET_ad(2,:));
end
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

% ET_ad_up_s(1,:) = smooth(ET_ad_up(1,:), 'sgolay', 3);
% ET_ad_up_s(2,:) = smooth(ET_ad_up(2,:), 'sgolay', 3);

%% SACCADE DETECTION (from detect_saccades_v2
fprintf('Detecting saccades\n');

sac_thresh = 3; %threshold eye speed Eyelink0 means it should be set to [1.2] 
peri_thresh = 1.5; %threshold eye speed for defining saccade boundary inds
min_isi = 0.1; max_isi = Inf; %min/max inter-saccade intervals

eye_smooth=3; eye_dt=1;
sm_avg_eyepos = ET_ad_up';

%for amplitude correction
amp_cutoff=8; %arcmins
win=[3 8];

%sm_avg_eyepos(:,1) = smooth(sm_avg_eyepos(:,1),eye_smooth);
%sm_avg_eyepos(:,2) = smooth(sm_avg_eyepos(:,2),eye_smooth);

eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]./eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]./eye_dt;
all_eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
ET_n=length(all_eye_speed);

% parameters
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

sprintf('detected %d saccades: %0.1f per second \n',length(saccade_inds), length(saccade_inds)./(length(ET_ad_up)./et_params.eye_fs))
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
if spk_offset ~=0
	LFP_times=LFP_times+spk_offset;
end

%% previous location of reading in kilosort outputs

%%
exptDataP = []; exptDataP2=[]; exptDataMUA=[];
if useofflinesorting==1
	for iUnit=1:length(spk_labels_SU)
		exptDataP(iUnit).spkID = double(spk_ID_SU(iUnit));
		exptDataP(iUnit).spkCh = spk_channels_SU(iUnit);
		% exptDataP(iUnit).rating = spk_rating_SU(iUnit);  
		exptDataP(iUnit).unit1 = spk_times(find(spk_clusters==spk_ID_SU(iUnit)));
	end

	for iUnit=1:length(spk_labels_MU)
		exptDataMUA(iUnit).spkID = double(spk_ID_MU(iUnit));
		exptDataMUA(iUnit).spkCh = spk_channels_MU(iUnit);
		% exptDataMUA(iUnit).rating = spk_rating_MU(iUnit);
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

		%         if useofflinesorting==1
		%             if ismember(channel,iso_SUs)
		%                 exptDataMUA{channel} = Clusters{channel}.times(find(Clusters{channel}.spike_clusts==1));
		%             else
		%                 exptDataMUA{channel} = Clusters{channel}.times(find(Clusters{channel}.spike_clusts>=1));
		%             end
		%         end  
	end
end

iSessions = 1;
allStartTS = [ExperimentRecording{:,2}];

%%
% iSessions = sessionsStartIDX(iSessions);
plexonSyncStrobesInThisSessionTS = events.timeStamps;
plexonSyncStrobesInThisSessionStrobeID = events.strobeNumber;
plexonSyncStrobesInThisSessionTS = plexonSyncStrobesInThisSessionTS(plexonSyncStrobesInThisSessionStrobeID == syncStrobeID);

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
			g_strctEyeCalib.EyeRaw.TimeStamp < thisTrialFlipon + g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(2);
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
	%    [ExperimentRecording{iTrials,1}.PlexonOnsetTime = (ExperimentRecording{iTrials, 2} - kofikoSyncTime)]

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
			%  spikesInThisTrial.rating = exptDataP(iUnit).rating;
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
				% nameOfUnit = ['unit',num2str(unitsInFile(iUnit))];            
         
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
		plexonLFPDataAlignedToThisTrial = LFP_ad(1:24, LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
			LFP_times - plexonSyncTime - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
		ExperimentRecording{iTrials,4} = plexonLFPDataAlignedToThisTrial;  
	end


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

	%     for iUnit2 = 1:allNumUnits2
	%         SortedSpikesAlignedToThisTrial = [];
	%         nameOfUnit = ['unit',num2str(iUnit2)];
	% %            nameOfUnit = ['unit',num2str(unitsInFile(iUnit))];            
	%         SortedSpikesAlignedToThisTrial =  exptDataP2.(nameOfUnit);
	%         SortedSpikesInThisTrial.(nameOfUnit) = SortedSpikesAlignedToThisTrial(SortedSpikesAlignedToThisTrial - ...
	%             (ExperimentRecording{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
	%             SortedSpikesAlignedToThisTrial - (ExperimentRecording{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
	%             - (ExperimentRecording{iTrials, 2} - kofikoSyncTime);
	%         ExperimentRecording{iTrials,9}(iUnit2) = length(SortedSpikesInThisTrial);
	%   
	%     end
	%     ExperimentRecording{iTrials,8} = SortedSpikesInThisTrial;
	trialIter = trialIter + 1;  
	if mod(iTrials,50)==0
		fprintf('finished trial %d of %d \n',iTrials, ntrials)
	end
end 
%%
% if ~isempty(experimentIndex)
%     save([strExperimentPath,filesep,'experimentIndex.mat'],'experimentIndex')
%     globalSessionTags = experimentIndex{:,3};
%     save([strExperimentPath,filesep,'experimentTags.mat'],'globalSessionTags')
% end
%%
cd(currentDirectory);
disp('done. saving...')

save([strExperimentPath,filesep,'FullExperimentRecordingwLFP_ks' num2str(useofflinesorting) '_' arraylabel '_v08.mat'],'ExperimentRecording','-v7.3')
disp('saved!')

%% check eye traces
% for iFixationCheck = 1:ntrials; %trialsInThisSession{iSessions};%
% figure(1); 
% plot(ExperimentRecording{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
% plot(ExperimentRecording{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
% title(['Processed ET signal - fixation target @ ' num2str(ExperimentRecording{iFixationCheck, 1}.m_pt2iFixationSpot)]);
% figure(2); 
% plot(ExperimentRecording{iFixationCheck,1}.ETthisTrialRawEyeData); title('Raw ET signal')
% figure(3);
% plot(ExperimentRecording{iFixationCheck,5}'); title('Plexon ET signal')
% pause
% end

%% Now we process LFPs
targ_trials=[];
for tt=1:size(ExperimentRecording,1)
	if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Dual Stim') && ...
			ExperimentRecording{tt, 1}.m_bMonkeyFixated==1
		targ_trials=[targ_trials,tt];    
	end
end
%targ_trials=targ_trials(end-1000:end); %for only picking the first 1000 trials to compare over time

%% Get CSD
%all_lfps=zeros(nChans,800,size(ExperimentRecording,1)-2);
all_lfps=zeros(nChans,800,length(targ_trials));

secCutoff=2.1; %seconds duration in between trials
iii=1;
%for tt=2:size(all_lfps,3)-1
for tt=targ_trials
	%if ExperimentRecording{tt, 1}.CSDtrigframe==1
	if ~isempty(ExperimentRecording{tt,4})   % BANDAID to get rid of bad first trial in ExperimentRecord
		all_lfps(:,:,iii)=ExperimentRecording{tt,4}(:,1:800);
		all_lfptimes(iii)=ExperimentRecording{tt,2};
		iii=iii+1;        
	end
end
iii

%%
cur_lfp_inds=[1:iii-1];
%cur_lfp_inds=[160:iii-1];
twin = all_lfptimes([cur_lfp_inds(1) cur_lfp_inds(end)])/60

%%
vars.Fs = 1000;
vars.BrainBound = 1;
vars.ChanSep = 0.05;
vars.diam = 2; %0.5
CSD = PettersenCSD(all_lfps(:,:,cur_lfp_inds),'spline',vars);

%%
expt_CSDs = mean(CSD,3);
CSDfig=figure;
imagesc(expt_CSDs)
title([filenameP ' iCSD'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 600]); colorbar
%ax=gca; ax.XTickLabel={'-100','0','100','200','300','400','500'};
ax=gca; ax.XTickLabel={'0','100','200','300','400','500','600'};
%caxis([-4e5 1e5])
%xlim([0 500]);
%ax=gca; ax.XTickLabel={'-50','50','150','250','350','450'};
%colormap(hotcold); caxis([-max(expt_CSDs(:))/3 max(expt_CSDs(:))/3]);
saveas(CSDfig,[strExperimentPath 'iCSD.pdf'])
saveas(CSDfig,[strExperimentPath 'iCSD.png'])

CSDfig2=figure;
imagesc(expt_CSDs)
title([filenameP ' iCSD'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 250]); colorbar
ax=gca; ax.XTickLabel={'0','50','100','150','200','250'};
saveas(CSDfig2,[strExperimentPath 'iCSD_zoom.pdf'])
saveas(CSDfig2,[strExperimentPath 'iCSD_zoom.png'])


%%
lfps2=squeeze(nanmean(all_lfps(:,:,cur_lfp_inds), 3));
LFPfig=figure;
imagesc(lfps2)
title([filenameP ' LFP'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 600]); colorbar
%ax=gca; ax.XTickLabel={'-100','0','100','200','300','400','500'};
ax=gca; ax.XTickLabel={'0','100','200','300','400','500','600'};
%colormap(hotcold); caxis([-max(lfps2(:))/3 max(lfps2(:))/3]);
saveas(LFPfig,[strExperimentPath 'LFP.pdf'])
saveas(LFPfig,[strExperimentPath 'LFP.png'])

%% 
% %% time-frequency plot
% % [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% %                  newtimef(data, frames, epochlim, srate, cycles,...
% %                       'key1',value1, 'key2',value2, ... );
% figure;
% chan=1;
% [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
%                   newtimef(all_lfps(chan,:,:), 800, [0 800], 1000, 0, 'freqs', [2 100]);
% %%
% for chan=[1:24];
%     [erspmat(chan,:,:),itc,powbase,times,freqs,erspboot,itcboot] = ...
% newtimef(all_lfps(chan,:,:), 800,[0 800], 1000, 0, 'freqs', [2 100]);
% figtitle(['Color Trial time-freq for channel ' num2str(chan)])
% %pause
% end
% % %%
% % figure; 
% % ersp_mean=squeeze(mean(erspmat(:,:,10:40),3));
% % imagesc(ersp_mean)
% % freqvec={num2str(freqs(5)), num2str(freqs(10)), num2str(freqs(15)), num2str(freqs(20))};
% % ax=gca;
% % ax.XTickLabel=(freqvec);
% 
% %%
% TFfig=figure; 
% ersp_mean=squeeze(erspmat(:,:,30));
% imagesc(ersp_mean)
% freqvec={num2str(freqs(5)), num2str(freqs(10)), num2str(freqs(15)), num2str(freqs(20))};
% ax=gca;
% ax.XTickLabel=(freqvec);
% ylabel('Channel'); xlabel('frequency')
% 
% saveas(TFfig,[strExperimentPath 'LFP_timefreq.png'])
% 
% disp('alignment - done!')
