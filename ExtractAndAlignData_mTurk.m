function [ExptTrials, ExptInfo] = ExtractAndAlignData_mTurk( exptname, dirpath, which_computer )
%
% Usage: ExtractAndAlignData( exptname, dirpath, <which_computer> )
%
% exptname of the form '230510_141725_Jacomo'
% dirpath is directory where these files are
%
% All these things directories need to be in path (or add subdirctory for Felix's dependencies folder)
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/CSD_functions/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilosort2/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Plexon-Matlab Offline Files SDK/')
% addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')

if nargin < 3
	% This can be used to set default directories
	% Dan's laptop = 0
	% Bevil office desktop = 1
    % LSR 2A58 sorting rig = 2
	which_computer = 2; % default value
end

if nargin < 2
	switch(which_computer)
		case 0, dirpath = '/Users/dbutts/Data/Conway/';
		case 1, dirpath = 'C:\SpkSort2023\AAActiveData\';
        case 2, dirpath = '/home/conwaylab/Data/';
		otherwise
			disp('which_computer is not specified')
	end
else
	if dirpath(end) ~= filesep
		dirpath = [dirpath filesep];
	end
end

%ExperimentFolder = dirpath;  % EXPERIMENT FOLDER
%filenameP = exptname;

%% Extract Bevil's data
spk_offset=0; 
ET_Eyelink=3; % set to 3 for any dDPI tracked files

useofflinesorting = 1; % set to 1 in order to use kilosort outputs, otherwise 0

skipLFP=0; 
nChans=24; %for CSD only now
LFPchans=1:nChans;
%LFP_chanmap=LFPchans;

%sessionsToProcess = [];
plxFilePath = [dirpath exptname '.pl2'];
matFilePath = [dirpath exptname filesep]; 
KSstitched=0; 
ksFilePath = [dirpath exptname filesep 'kilosorting_laminar' filesep]; 
arraylabel ='lam';
%strExperimentPath = [matFilePath 'Analysis' filesep];
output_directory = [matFilePath 'Analysis' filesep];

configFilePath = [dirpath exptname '.mat'];

%sessionTimeOffsets = [];
if ~exist(output_directory,'dir')
	mkdir(output_directory);
end
%filenameK = matFilePath;
load(configFilePath)

%% test for misalignment
%/{
[SPKC_adfreq, SPKC_n, SPKC_ts, SPKC_fn, SPKC_ad_test] = plx_ad_v(plxFilePath, 'SPKC003');
[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad_test] = plx_ad_v(plxFilePath, 'AI07');
[LFP_adfreq, LFP_n, LFP_ts, LFP_fn, LFP_ad_test] = plx_ad_v(plxFilePath, ['FP' num2str(1,'%03.f')]);
ET_n./LFP_n;
ET_n-LFP_n;
if abs(ET_n-LFP_n)>100
	warning('Danger - Plexon timing bug dectected')
else
	disp('Plexon timing check passed')
end
%}

%%
global strctColorValues unitsInFile PlottingVars g_strctStatistics ExptTrials topLevelIndex 

%persistent plexonDataAlignedToThisTrial
experimentIndex = {};

g_strctStatistics.preTrialWindow = 0;
g_strctStatistics.postTrialWindow = 6; %0.4

currentDirectory = pwd;
cd(dirpath);
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

%% Prepare Kilosort information
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

%% Load Kofiko files
allMatFiles = dir([matFilePath, '*.mat']);
[~, indices] = sort(vertcat(allMatFiles(:).datenum));
allMatFiles = allMatFiles(indices);
g_strctStatistics.ExptTrials = {};
for iFiles = 1:size(allMatFiles,1)
	if allMatFiles(iFiles).isdir
		continue
	end
	load([matFilePath,allMatFiles(iFiles).name])
	fprintf('-> loading kofiko file %s\n', allMatFiles(iFiles).name);
	if exist('g_strctLocalExperimentRecording') && size(g_strctLocalExperimentRecording{1},1) > 1
		warning('incorrect format detected in save structure in file  %s, skipping', allMatFiles(iFiles).name)
		%tmp = g_strctLocalExperimentRecording{1};
		continue
	end

	% details = whos([output_directory,allMatFiles(iFiles).name])
	%g_strctLocalg_strctStatistics.ExptTrials(cellfun('isempty',g_strctLocalg_strctStatistics.ExptTrials)) = [];
	g_strctStatistics.ExptTrials(cellfun('isempty',g_strctStatistics.ExptTrials)) = [];
	%g_strctLocalg_strctStatistics.ExptTrials = vertcat(g_strctLocalg_strctStatistics.ExptTrials,g_strctDynamicStimLog.TrialLog);
	%g_strctLocalg_strctStatistics.ExptTrials = vertcat(g_strctLocalg_strctStatistics.ExptTrials{cellfun(@isempty,g_strctLocalg_strctStatistics.ExptTrials),dataToSave');
	if (exist('g_strctLocalExperimentRecording') == 1)
		% disp('case 1')
		g_strctStatistics.ExptTrials = vertcat(g_strctStatistics.ExptTrials,vertcat({g_strctLocalExperimentRecording{find(~cellfun(@isempty,g_strctLocalExperimentRecording))}})');
		clear g_strctLocalExperimentRecording
	elseif (exist('strctLocalExperimentRecording') == 1)
		%disp('case 2')
		g_strctStatistics.ExptTrials = vertcat(g_strctStatistics.ExptTrials,vertcat({strctLocalExperimentRecording{find(~cellfun(@isempty,strctLocalExperimentRecording))}})');
		clear('strctLocalExperimentRecording')
	else
		%disp('case 3')
		try
			g_strctStatistics.ExptTrials = vertcat(g_strctStatistics.ExptTrials,vertcat({dataToSave{find(~cellfun(@isempty,dataToSave))}})');
		catch
			warning(sprintf('file %s did not contain appropriately formatted trials', allMatFiles(iFiles).name)); 
		end
	end
end

ExptTrials = g_strctStatistics.ExptTrials;
clear g_strctStatistics.ExptTrials;

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
cd(dirpath);

try
load(matFilePath, 'g_strctDAQParams');
catch
	try
		load([dirpath, exptname, '.mat'], 'g_strctDAQParams');
	catch
		[~, filename] = fileparts(matFilePath); 
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

if strfind(exptname, '190409_134012_Jacomo')
	numSessions=8;
end
sessionStartPlexonIDX = find(events.strobeNumber == startRecordID);
sessionEndPlexonIDX = find(events.strobeNumber == stopRecordID);
sessionsStartPlexonTS = events.timeStamps(sessionStartPlexonIDX);
sessionEndPlexonTS = events.timeStamps(sessionEndPlexonIDX);
sessionsStartID = 1;

%% Establish trial times based on Kofiko-time stamp of image-flip
for iTrials = 1:size(ExptTrials,1)
	ExptTrials{iTrials, 2} = ExptTrials{iTrials, 1}.m_fImageFlipON_TS_Kofiko;
	% DAN COMMENTED OUT: UNUSED trialIter(iTrials) = ExptTrials{iTrials, 1}.m_iTrialNumber;

	% Add session information (SessionID) which is a manual way in Kofiko to annotate
	% Note this is not strict and not every trial will have one
	% Every "start of recording" (a button) in Kofiko is a different session time 
	ExptTrials{iTrials, 1}.SessionID = find(ExptTrials{iTrials, 2} > sessionStartTS,1,'last');
end

%% To picking a particular subset of trials -- including checking timing validity
% Use timing of beginning and end of Kofiko time stamps to remove erroneous trial information  
RECstart = kofikoStrobeAllTS(1);
RECend = kofikoStrobeAllTS(end);

% targ_trials=[];
% for tt=1:length(ExptTrials)
% 	if strcmp(ExptTrials{tt, 1}.m_strTrialType, 'Dual Stim')
% 		% if strcmp(ExptTrials{tt, 1}.m_strTrialType, 'Fivedot');
% 		if (ExptTrials{tt, 2} >= RECstart) && (ExptTrials{tt, 2} <= RECend) 
% 			targ_trials = [targ_trials, tt];
% 		else
% 			fprintf('   Invalid Dual-Stim trial detected (#%d) -- eliminating.\n', tt)
% 		end
% 	end
% end
% 
% %% Reduce trials to valid subset
% ExptTrials=ExptTrials(targ_trials,:);

%% Ensure trials are in correct order
ntrials=length(ExptTrials);

for tt=1:ntrials
	trialstart_raw(tt)=ExptTrials{tt, 2}; 
end
[~,order] = sort(trialstart_raw);
ExptTrials = ExptTrials(order, :);  

%%
topLevelIndex = [];
sessionIndex  = [];
sessionTags   = {};

cd ..
allPLXfiles = vertcat(dir([pwd,filesep, '*.plx']),dir([pwd,filesep, '*.pl2']));

thisSessionFile = plxFilePath;
if ~exist(plxFilePath) && isempty(allPLXfiles) && ~any(~arrayfun(@isempty,strfind({allPLXfiles(:).name}, [exptname,'.plx']))) && ~any(~arrayfun(@isempty,strfind({allPLXfiles(:).name}, [exptname,'.pl2'])))
	warning('could not find PLX file for this experiment')
else
	%  thisSessionFile = [pwd,'\',exptname];
end

cd(dirpath)
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

%% This goes through online sorting -- commenting out
%{
for channel=1:nChans
	numUnitsInSession = 0;
	for iUnit = 1:numel(numUnits)
		nameOfUnit = ['unit',num2str(numUnits(iUnit))];
		tempTS = 0;
		try
			%%%% [~, ~, tempTS, ~] = plx_waves_v([thisSessionFile], channel, iUnit);  % ONLY FOR ONLINE SORTING
			% [~, ~, tempTS, ~] = PL2Waves([thisSessionFile], channel, iUnit);
		catch
			fprintf('big oof \n')
		end
		if tempTS > 0
			numUnitsInSession = numUnitsInSession + 1;
			unitsInFile = [unitsInFile, iUnit];
			[spikes(channel).(nameOfUnit).count, spikes(channel).(nameOfUnit).numWaves, spikes(channel).(nameOfUnit).timeStamps, spikes(channel).(nameOfUnit).Waves] = ...
				plx_waves_v([thisSessionFile], channel, iUnit);
			% [spikes.(nameOfUnit).count, spikes.(nameOfUnit).numWaves, spikes.(nameOfUnit).timeStamps, spikes.(nameOfUnit).Waves] = ...
			%		PL2Waves([thisSessionFile], channel, iUnit);
		end
	end
	allNumUnits=[allNumUnits,numUnitsInSession];
end
%}
%numUnitsInSession=max(allNumUnits);  % not used

%% Get ET data
load([dirpath, exptname, '.mat'], 'g_strctEyeCalib');
load([dirpath, exptname, '.mat'], 'g_strctStimulusServer');

if ET_Eyelink == 1 %binocular Eyelink
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI05');
	[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI06');
	[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);
	PlexET_ad_calib(3,:) = (PlexET_ad_calib(3,:)-median(PlexET_ad_calib(3,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(4,:) = (PlexET_ad_calib(4,:)-median(PlexET_ad_calib(4,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 2 %monocular Eyelink
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 3 % DDPI
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI03'); %contains the DDPI trigger information we need for alignment
	[~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI04');
	[~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad_calib=PlexET_ad;
	PlexET_ad_calib(3,:) = (PlexET_ad_calib(3,:)-median(PlexET_ad_calib(3,:)))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	PlexET_ad_calib(4,:) = (PlexET_ad_calib(4,:)-median(PlexET_ad_calib(4,:)))*(g_strctEyeCalib.GainY.Buffer(end)./1000);

elseif ET_Eyelink == 0 % Iscan
	[~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
	[ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
	PlexET_ad(1,:)=PlexET_ad(1,:)-median(PlexET_ad(1,:));
	PlexET_ad(2,:)=PlexET_ad(2,:)-median(PlexET_ad(2,:));
end
PlexET_times=[1:ET_n]/ET_adfreq;


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
	ET_ad(1,CalibEyeBufferIDX) = (ET_ad(1,CalibEyeBufferIDX) - g_strctEyeCalib.CenterX.Buffer(eyeCenterXID)) * g_strctEyeCalib.GainX.Buffer(eyeGainXID);% + (g_strctStimulusServer.m_aiScreenSize(3)/2)-ExptTrials{1, 1}.m_pt2iFixationSpot(1);
	ET_ad(2,CalibEyeBufferIDX) = (ET_ad(2,CalibEyeBufferIDX) - g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)) * g_strctEyeCalib.GainY.Buffer(eyeGainYID);% + (g_strctStimulusServer.m_aiScreenSize(4)/2)-ExptTrials{1, 1}.m_pt2iFixationSpot(2);
end

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

fprintf('-> detected %d saccades: %0.1f per second\n',length(saccade_inds), ...
	length(saccade_inds)./(length(ET_ad_up)./et_params.eye_fs))

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
fprintf('Done with rudimentary microsaccade detect\n\n')

%% PROCESS LFPs
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

%% Organize spike time metadata
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
    %this code is how we would extract onlines-roted spikes
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
allStartTS = [ExptTrials{:,2}];

%%
% iSessions = sessionsStartIDX(iSessions);
plexonSyncStrobesInThisSessionTS = events.timeStamps;
plexonSyncStrobesInThisSessionStrobeID =   events.strobeNumber;
plexonSyncStrobesInThisSessionTS =  plexonSyncStrobesInThisSessionTS(plexonSyncStrobesInThisSessionStrobeID == syncStrobeID);

kofikoSyncStrobesInThisSessionTS = g_strctDAQParams.LastStrobe.TimeStamp;
kofikoSyncStrobesInThisSessionStrobeID = g_strctDAQParams.LastStrobe.Buffer;
kofikoSyncStrobesInThisSessionTS = kofikoSyncStrobesInThisSessionTS(kofikoSyncStrobesInThisSessionStrobeID == syncStrobeID);
if isempty(kofikoSyncStrobesInThisSessionTS) || isempty(plexonSyncStrobesInThisSessionTS)
	sprintf('could not process recording %s, session %i, timestamp missing.', exptname, iSessions)
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
%[ExptTrials{:,12}] == 1);
trialsInThisSession{iSessions} = find(allStartTS > sessionStartTS(1) & allStartTS < sessionEndTS(end));

trialIter = 1;

for iFixationCheck = 1:ntrials  %trialsInThisSession{iSessions};%
	if isfield(ExptTrials{iFixationCheck, 1},'m_bMonkeyFixated') && ExptTrials{iFixationCheck, 1}.m_bMonkeyFixated
		ExptTrials{iFixationCheck,7} = 1;
		sessionIndex(trialIter).m_bFixated = 1;
		sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = [];

		% append eye trace information to this trial
		thisTrialFlipon = ExptTrials{iFixationCheck,2};
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
		ExptTrials{(iFixationCheck),1}.m_afEyeXPositionScreenCoordinates = fEyeXPix;
		ExptTrials{(iFixationCheck),1}.m_afEyeYPositionScreenCoordinates = fEyeYPix;
		ExptTrials{iFixationCheck,1}.m_afEyePositiontimes = thisTrialRawEyeDatatimes;
		ExptTrials{iFixationCheck,1}.ETthisTrialRawEyeData = thisTrialRawEyeData;
		ExptTrials{iFixationCheck,1}.ETCenter = [g_strctEyeCalib.CenterX.Buffer(eyeCenterXID),g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)];
		ExptTrials{iFixationCheck,1}.ScreenCenter = [(g_strctStimulusServer.m_aiScreenSize(3)/2),(g_strctStimulusServer.m_aiScreenSize(4)/2)];
		sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = fEyeXPix;
		sessionIndex(trialIter).m_afEyeYPositionScreenCoordinates = fEyeYPix;
    
	elseif ~isfield(ExptTrials{(iFixationCheck), 1},'m_bMonkeyFixated') || ~ExptTrials{(iFixationCheck), 1}.m_bMonkeyFixated
		ExptTrials{(iFixationCheck),7} = 0;
		sessionIndex(trialIter).m_bFixated = 0;
		sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = [];        % append eye trace information to this trial
		thisTrialFlipon = ExptTrials{(iFixationCheck),2} ;
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
		ExptTrials{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates = fEyeXPix;
		ExptTrials{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates = fEyeYPix;
		ExptTrials{iFixationCheck,1}.m_afEyePositiontimes = thisTrialRawEyeDatatimes;
		ExptTrials{iFixationCheck,1}.ETthisTrialRawEyeData = thisTrialRawEyeData;
		ExptTrials{iFixationCheck,1}.ETCenter = [g_strctEyeCalib.CenterX.Buffer(eyeCenterXID),g_strctEyeCalib.CenterY.Buffer(eyeCenterYID)];
		ExptTrials{iFixationCheck,1}.ScreenCenter = [(g_strctStimulusServer.m_aiScreenSize(3)/2),(g_strctStimulusServer.m_aiScreenSize(4)/2)];
		sessionIndex(trialIter).m_afEyeXPositionScreenCoordinates = fEyeXPix;
		sessionIndex(trialIter).m_afEyeYPositionScreenCoordinates = fEyeYPix;
	else
		spintf('big oof: no ET data')
	end
	trialIter = trialIter + 1;
%     figure(1);
%     plot(ExptTrials{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
%     plot(ExptTrials{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
% 
%     figure(3);
%     plot(ExptTrials{iFixationCheck,5}')
%     pause
end

%%
%{
for iFixationCheck = 1:ntrials; %trialsInThisSession{iSessions};%
figure(1); 
plot(ExptTrials{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
plot(ExptTrials{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
figure(2); 
plot(ExptTrials{iFixationCheck,1}.ETthisTrialRawEyeData);
pause
end
%}

%% Align all time signals
trialIter = 1;
disp(['General time-alignment: ' num2str(ntrials) ' trials'])
for iTrials = 1:ntrials   %trialsInThisSession{iSessions}
	spikesInThisTrial = [];
	[~,trialSyncStrobeID] = min(abs(ExptTrials{iTrials, 2} - kofikoSyncStrobesInThisSessionTS));

	if trialSyncStrobeID>length(kofikoSyncStrobesInThisSessionTS)
		error(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(kofikoSyncStrobesInThisSessionTS))]);
		%trialSyncStrobeID=length(kofikoSyncStrobesInThisSessionTS);
	end
  
	if trialSyncStrobeID>length(plexonSyncStrobesInThisSessionTS)
		error(['Error with sync strobes? off by ' num2str(trialSyncStrobeID-length(plexonSyncStrobesInThisSessionTS))]);
		%trialSyncStrobeID=length(plexonSyncStrobesInThisSessionTS);
	end

    %append synchronization info
	kofikoSyncTime = kofikoSyncStrobesInThisSessionTS(trialSyncStrobeID);
	plexonSyncTime = plexonSyncStrobesInThisSessionTS(trialSyncStrobeID);
	ExptTrials{iTrials,1}.kofikoSyncTime    = kofikoSyncTime; 
	ExptTrials{iTrials,1}.PlexonSyncTime    = plexonSyncTime;
	ExptTrials{iTrials,1}.PlexonOnsetTime   = plexonSyncTime + (ExptTrials{iTrials, 2} - kofikoSyncTime);
	ExptTrials{iTrials,3}                   = plexonSyncTime + (ExptTrials{iTrials, 2} - kofikoSyncTime); %this gets its own field to make it easier to read out trial onset information on the plexon clock 

	if ~isfield( ExptTrials{iTrials, 1},'m_aiStimColor') || isempty(ExptTrials{iTrials, 1}.m_aiStimColor)
		ExptTrials{iTrials, 1}.m_aiStimColor = [NaN, NaN, NaN];
	end

%	ExptTrials{iTrials, 4} = ExptTrials{iTrials, 1}.m_strTrialType;

	if useofflinesorting==1
		for iUnit=1:nSU
			plexonDataAlignedToThisTrial = [];         
			plexonDataAlignedToThisTrial =  exptDataP(iUnit).unit1 - plexonSyncTime;
      
			spikesInThisTrial.unit1 = plexonDataAlignedToThisTrial(plexonDataAlignedToThisTrial - ...
				(ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
				plexonDataAlignedToThisTrial - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow) ...
				- (ExptTrials{iTrials, 2} - kofikoSyncTime);
           
			spikesInThisTrial.spkID = exptDataP(iUnit).spkID;
			%  spikesInThisTrial.rating = exptDataP(iUnit).rating;
			ExptTrials{iTrials,10+2*(iUnit-1)} = spikesInThisTrial;
			ExptTrials{iTrials,11+2*(iUnit-1)} = length(spikesInThisTrial.unit1);
		end

		for iUnit=1:nMU
			plexonMUADataAlignedToThisTrial=[];
			plexonMUADataAlignedToThisTrial =  exptDataMUA(iUnit).unit1 - plexonSyncTime ;          
			MUAspikesInThisTrial.unit1 = plexonMUADataAlignedToThisTrial(plexonMUADataAlignedToThisTrial - ...
				(ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
				plexonMUADataAlignedToThisTrial - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
				- (ExptTrials{iTrials, 2} - kofikoSyncTime);

			MUAspikesInThisTrial.spkID = exptDataMUA(iUnit).spkID;
			ExptTrials{iTrials,10+2*(iUnit-1+nSU)} = MUAspikesInThisTrial.unit1;
			ExptTrials{iTrials,11+2*(iUnit-1+nSU)} = length(MUAspikesInThisTrial.unit1);	
		end
	else		
		for channel = find(allNumUnits) %1:nChans
			for iUnit = 1:allNumUnits(channel)
				plexonDataAlignedToThisTrial = []; 
				nameOfUnit = ['unit',num2str(iUnit)];
				% nameOfUnit = ['unit',num2str(unitsInFile(iUnit))];            
         
				plexonDataAlignedToThisTrial =  exptDataP(channel).(nameOfUnit) - plexonSyncTime ;
				spikesInThisTrial.(nameOfUnit) = plexonDataAlignedToThisTrial(plexonDataAlignedToThisTrial - ...
					(ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
					plexonDataAlignedToThisTrial - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
					- (ExptTrials{iTrials, 2} - kofikoSyncTime);;
			end
			%         plexonMUADataAlignedToThisTrial =  exptDataMUA{channel} - plexonSyncTime ;
			%         MUAspikesInThisTrial = plexonMUADataAlignedToThisTrial(plexonMUADataAlignedToThisTrial - ...
			%             (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
			%             plexonMUADataAlignedToThisTrial - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
			%             - (ExptTrials{iTrials, 2} - kofikoSyncTime);
			ExptTrials{iTrials,10+2*(channel-1)} = spikesInThisTrial;
			ExptTrials{iTrials,11+2*(channel-1)} = length(spikesInThisTrial.unit1);
			%        ExptTrials{iTrials,12+3*(channel-1)} = MUAspikesInThisTrial;
		end 
	end

	if ~skipLFP
		%    ncols=size(ExptTrials,2);
		plexonLFPDataAlignedToThisTrial = LFP_ad(1:24, LFP_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
			LFP_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
		ExptTrials{iTrials,5} = plexonLFPDataAlignedToThisTrial;  
	end


	% % for aligning EMs to plexon timestamps    
	%     plexonEMsAlignedToThisTrial =  saccade_times - plexonSyncTime ;
	%     EMsInThisTrial = plexonEMsAlignedToThisTrial(plexonEMsAlignedToThisTrial - ...
	%         (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
	%         plexonEMsAlignedToThisTrial - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow)...
	%         - (ExptTrials{iTrials, 2} - kofikoSyncTime);

	% for aligning EMs to Kofiko timestamps
	EMsInThisTrial.saccades =  saccade_times(saccade_times >= ExptTrials{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
		saccade_times <= ExptTrials{iTrials, 2}+g_strctStatistics.postTrialWindow);
	EMsInThisTrial.saccade_start =  sac_start_times(saccade_times >= ExptTrials{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
		sac_start_times <= ExptTrials{iTrials, 2}+g_strctStatistics.postTrialWindow);
	EMsInThisTrial.saccade_stop =  sac_stop_times(saccade_times >= ExptTrials{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
		sac_stop_times <= ExptTrials{iTrials, 2}+g_strctStatistics.postTrialWindow);
	EMsInThisTrial.ET_times =  ET_times(ET_times >= ExptTrials{iTrials, 2}+g_strctStatistics.preTrialWindow & ...
		ET_times <= ExptTrials{iTrials, 2}+g_strctStatistics.postTrialWindow);
	ExptTrials{iTrials,6} = EMsInThisTrial;

	%     plexonETDataAlignedToThisTrial.ET_trace = PlexET_ad(:, PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
	%             PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
	%     plexonETDataAlignedToThisTrial.ET_times = PlexET_times(:, PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
	%             PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);

	plexonETDataAlignedToThisTrial.ET_trace = PlexET_ad(:, PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
		PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);
	plexonETDataAlignedToThisTrial.ET_times = PlexET_times(:, PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) >= g_strctStatistics.preTrialWindow & ...
		PlexET_times - plexonSyncTime - (ExptTrials{iTrials, 2} - kofikoSyncTime) <=  g_strctStatistics.postTrialWindow);  
	ExptTrials{iTrials,7} = plexonETDataAlignedToThisTrial;

	trialIter = trialIter + 1;  
	if mod(iTrials,100)==0
		fprintf('-> trial %4d of %d\n',iTrials, ntrials)
	end
end 
%%
% if ~isempty(experimentIndex)
%     save([output_directory,'experimentIndex.mat'],'experimentIndex')
%     globalSessionTags = experimentIndex{:,3};
%     save([output_directory,'experimentTags.mat'],'globalSessionTags')
% end

fprintf('Done combining information.\n')

%%
cd(currentDirectory);

outfile = sprintf( '%s_FullExpt_ks%d_%s_v09.mat', exptname, useofflinesorting, arraylabel );

% Save other variables to continue the process
ExptInfo.exptname = exptname;
ExptInfo.array_label ='lam';
ExptInfo.nSU = nSU;
ExptInfo.nMU = nMU;
ExptInfo.spk_ID_SU = spk_ID_SU; 
ExptInfo.spk_ID_MU = spk_ID_MU;
ExptInfo.spk_channels_SU = spk_channels_SU;
ExptInfo.spk_channels_MU = spk_channels_MU;
ExptInfo.ET_Eyelink = ET_Eyelink;
ExptInfo.use_offline_sorting = useofflinesorting;
ExptInfo.expt_folder = dirpath;

% Add experiment metadata from base file (moved from cloud stim 2: general expt parameters here)
load([dirpath exptname '.mat'], 'g_astrctAllParadigms')
ExptInfo.g_astrctAllParadigms = g_astrctAllParadigms{1};  % this is only one cell array -- what is it?
ExptInfo.g_strctEyeCalib = g_strctEyeCalib;
ExptInfo.g_strctDAQParams = g_strctDAQParams;  % don't know if needed, but just in case

fprintf('Saving %s in %s....\n', outfile, output_directory)
save([output_directory outfile], 'ExptTrials', 'ExptInfo', '-v7.3')  % HDF5 output

% Output RAW plexon ET trace in separate file (for time-alignment etc)
EToutfile = sprintf( '%s_FullExpt_ET.mat', exptname );

fprintf('Saving %s in %s....\n', EToutfile, output_directory)

PlexET_ad_calib = PlexET_ad_calib';
PlexET_times = PlexET_times';
save([output_directory EToutfile], 'PlexET_ad_calib', 'PlexET_times', '-v7.3')  % HDF5 output
disp('Saved!')

%% check eye traces
% for iFixationCheck = 1:ntrials; %trialsInThisSession{iSessions};%
% figure(1); 
% plot(ExptTrials{iFixationCheck,1}.m_afEyeXPositionScreenCoordinates); hold on
% plot(ExptTrials{iFixationCheck,1}.m_afEyeYPositionScreenCoordinates); hold off
% title(['Processed ET signal - fixation target @ ' num2str(ExptTrials{iFixationCheck, 1}.m_pt2iFixationSpot)]);
% figure(2); 
% plot(ExptTrials{iFixationCheck,1}.ETthisTrialRawEyeData); title('Raw ET signal')
% figure(3);
% plot(ExptTrials{iFixationCheck,5}'); title('Plexon ET signal')
% pause
% end

%% Processing the LFPs/CSDs
fprintf('\nProcessing CSDs (all trials):\n')
CSDprocess( ExptTrials, [], 1, exptname, output_directory, 1 );
% Usage: [csd, lfp] = CSDprocess( ExptRecording, trial_select, save_images, exptname, savedir, verbose )
return 
