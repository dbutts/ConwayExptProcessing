function data = PackageTaskData_v9( exptdata, metadata_struct, task_type, targ_ETstimtype, stimFilePath, output_path, skipLFP, which_computer )
%
% Usage: data = PackageTaskData_v9( exptdata, metadata_struct, <task_type>, <targ_ETstimtype>, <stimFilePath>, output_path, skipLFP, which_computer )
%
% Can check experiment composition
% experiment_composition( exptdata );

% Process depending on task type
% 1 = old mTurk

% for ETstimtype: we use only 1, 7 (will work with both of these)
% 7 = Color cloud in both windows
% 1 = 1-D bars alternating each frame in both windows
% 2 = 1-D horizontal in 1, vertical in other
% 3 = Alternate to 2...

if (nargin < 3) || isempty(task_type)
	task_type = 1;    % this selects which stimulus to process: 8 looks for cloud stims
end
if (nargin < 3) || isempty(targ_ETstimtype)
	targ_ETstimtype = 0;    % zero means that there is no eye-tracking stim (or at least don't load
end

skipLP = 0; %set to 1 if you want to skip laminar probe data and only analyzing the ET stims
%skipET = 1; %set to 1 if you want to skip the ETstim (such as if you moved them all the way out of the way); set to on on 5/12/23 on bevil's request to not have to worry about this
% Just set targ_ETstimtype to zero (default if you dont want this
if targ_ETstimtype == 0
	skipET = 1;
end

plot_intermediate=0;

%% reload data for analysis
if nargin < 8
	% This can be used to set default directories
	% Dan's laptop = 0
	% Bevil office desktop = 1
	which_computer = 0; % default value
end

if (nargin < 7) || isempty(skipLTP)
	skipLFP = 1;
end

if (nargin < 6) 
	output_path = [];
end
if isempty(output_path)
	output_path = metadata_struct.expt_folder;
end

if (nargin < 5) || isempty(stimFilePath)
	switch(which_computer)
		case 0, stimFilePath = '/Users/dbutts/Data/Conway/Cloudstims_calib_01_2022/';
		case 1, stimFilePath = 'C:\SpkSort2023\Cloudstims_calib_01_2022\';
		otherwise
			disp('which_computer not defined')
	end
end

% Pull necessary variables from meta-data
exptname = metadata_struct.exptname;
arraylabel = metadata_struct.array_label;
useofflinesorting = metadata_struct.use_offline_sorting;
nSU = metadata_struct.nSU;
nMU = metadata_struct.nMU;
g_astrctAllParadigms = metadata_struct.g_astrctAllParadigms;  % expt configuration information
g_strctEyeCalib = metadata_struct.g_strctEyeCalib;
output_directory = [output_path filesep];

% Set up data-struct for output
data.exptname = exptname;
data.stim = [];   % these are the big variables that not worth representing twice
data.stimET = [];

%% params
if useofflinesorting == 1
	nChans=metadata_struct.nSU+metadata_struct.nMU;
	spk_ID=[metadata_struct.spk_ID_SU; metadata_struct.spk_ID_MU];
	spk_ch=[metadata_struct.spk_channels_SU+1; metadata_struct.spk_channels_MU+1];
else
	nChans=24; %targchans=1:24;
end

if any([strfind(exptname,'220203'), strfind(exptname,'220205'), strfind(exptname,'220207')])
	trlbins=160; dt=.024; repframes=3;
else
	trlbins=240; dt=.016; repframes=2;
end
%fixvarcutoff=6;
%blink_thresh = 60;
cur_BlockID=1; 

%%
numSpks = zeros(nChans);
for cc=1:nChans
	numSpks(cc)=sum(cell2mat(exptdata(:, 11+2*(cc-1))));  % this just seems to have total number of spikes
end

%targchans=find(numSpks>2000);

if useofflinesorting == 1
	targchans=1:nChans;  % pass all channels through
end

%% Assign list of trials to teach type of data from the experiment
trlonset_diffs = [4; diff(cell2mat(exptdata(:, 2)))];
exptdata_mod = exptdata;
%for switch_stimtype=0:8

	%targ_trials=[];  
	%for tt=1:length(exptdata)
		%    if strcmp(exptdata{tt, 1}.m_strTrialType, 'Dense Noise');
	%	if strcmp(exptdata{tt, 1}.m_strTrialType, 'Dual Stim') && ...
	%		(exptdata{tt, 1}.m_bMonkeyFixated == 1) && (trlonset_diffs(tt) > 0) && ...
	%		(exptdata{tt, 1}.DualstimPrimaryuseRGBCloud == switch_stimtype) % &&...
			% exptdata{tt, 1}.DualstimSecondaryUseCloud==targ_ETstimtype; % && exptdata{tt, 1}.m_aiStimulusRect(1)==975;
	%		targ_trials=[targ_trials,tt];
	%	end
	%end
	% targ_trials = 1:size(exptdata,1);
	%switch switch_stimtype
	%	case 0; curstimstype='GT';  exptdata_GT=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_GT','-v7.3')
	%	case 3; curstimstype='HL';  exptdata_HartleyLum=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_HartleyLum','-v7.3')
	%	case 6; curstimstype='HC';  exptdata_HartleyCol=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_HartleyCol','-v7.3')
	%	case 8; curstimstype='CC';    exptdata_ColCloud=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_ColCloud','-v7.3')
	%end

%end

%%
%{
for tt=1:length(exptdata)
figure(1);
plot(exptdata{tt, 1}.m_afEyeXPositionScreenCoordinates-(exptdata{tt, 1}.m_pt2iFixationSpot(1))); hold on
plot(exptdata{tt, 1}.m_afEyeYPositionScreenCoordinates-(exptdata{tt, 1}.m_pt2iFixationSpot(2))); hold off

title(num2str(exptdata{tt, 1}.m_bMonkeyFixated));
%exptdata{4353, 1}.m_afEyeXPositionScreenCoordinates  

figure(2);
plot(exptdata{tt, 5}');

pause
end
%}

% DAN commented out: This now comes in the experiment info metadata
%load([metadata_struct.expt_folder exptname '.mat'], 'g_astrctAllParadigms')
%DualstimETbars = int8(squeeze((metadata_struct.g_astrctAllParadigms.DualstimETbars-128)/127)');

%% Stimulus selection
%switch targ_stimtype
%	case 3; curstimstype='HL';	exptdata_mod = exptdata_HartleyLum;
%	case 6; curstimstype='HC';	exptdata_mod = exptdata_HartleyCol;
%	case 7; curstimstype='LC';	exptdata_mod = exptdata_LumCloud;
%	case 8; curstimstype='CC';	exptdata_mod = exptdata_ColCloud;
%end
switch targ_ETstimtype
	case 1; curETstimtype='1D';
	case 7; curETstimtype='CC';
end

%% Random assortment of setup things -- before going into trial-by-trial processing
%ntrls = min([size(exptdata_mod,1),500]); % only using the first 500
ntrls = size(exptdata_mod,1);
NT=trlbins*ntrls;

if any([strfind(exptname,'220203'), strfind(exptname,'220205'), strfind(exptname,'220207')])
	trlsecs=2.67;
	edges = linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);
	cur_trial_ETsamples = [0:.001:2.669];
else
	trlsecs = 6;
	edges = linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);  
	cur_trial_ETsamples = [0:.001:3.999];
end

if metadata_struct.ET_Eyelink == 3
	numETtraces = 2;
	ET_trace_raw_1khz = zeros(ntrls*trlsecs*1000, 3);
else
	numETtraces = size(exptdata_mod{end, 7}.ET_trace, 1);
	ET_trace_raw_1khz = zeros(ntrls*trlsecs*1000, numETtraces);
end
%ET_ad_up=[];

%data.stim = int8(zeros(NT,60,60,3)); 
stimtype = zeros(NT,1);  

if ~skipET
	if targ_ETstimtype == 1
		data.stimET = int8(zeros(NT,60)); 
		stimtypeET = zeros(NT,1); 
		stimETori = zeros(NT,1); 
	else
		data.stimET = int8(zeros(NT,60,60,3)); 
		stimtypeET = zeros(NT,1); 
	end
end

trialstart_plx=zeros(ntrls,1);

hartleystim_metas=zeros(NT,4); 
binned_SU1 = int8(zeros(NT,nChans)); binned_MUA=[]; ET_trace=zeros(NT,2); ET_trace_raw=zeros(NT,2);
%binned_SU2=int8(zeros(NT,nChans)); binned_SU3=int8(zeros(NT,nChans)); %binned_SU4=zeros(NT,nChans); 
bad_inds_fix = [1:10]; bad_inds_sac = []; sacc_inds=[];
Block_onsetinds=1; Block_offsetinds=[]; Block_onsetinds_blink=1; Block_offsetinds_blink=[];
BlockID=[]; TrialID=[]; useLeye=zeros(NT,1); useReye=zeros(NT,1);
%stimloc=[]; % DAB: old variable 
spk_times_all=cell([nSU+nMU,1]);

cloud_scale=zeros(NT,1); 
cloud_binary=zeros(NT,1);     

%cd(stimFilePath)
%cur_scale=g_astrctAllParadigms{1, 1}.DualstimScale.Buffer(find(g_astrctAllParadigms{1, 1}.DualstimScale.TimeStamp<exptdata_mod{2,2},1,'last'));
%cur_scale=g_astrctAllParadigms.DualstimScale.Buffer(find(g_astrctAllParadigms.DualstimScale.TimeStamp<exptdata_mod{2,2},1,'last'));
%load(sprintf('Cloudstims_Chrom_size60_scale%d_%02d.mat', cur_scale, cur_BlockID));
%load(sprintf([stimFilePath 'Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID)) 
%DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));

%% Load external eye-tracking file info -- has two variables (at 1 kHZ)
% NOT NEEDED -- plexon information already went into trial-level information
% load(sprintf('%s%s_FullExpt_ET.mat', metadata_struct.expt_folder, exptname ))
% variables: PlexET_ad_calib, PlexET_times 

%% CREATE TASK VARIABLES
% Classify each trial as correct or incorrect -- 
data.task_correct = zeros(ntrls);
data.task_type = zeros(ntrls);
data.task_outcome = {};
%data.task_catchtrial = zeros(ntrls);
data.cueID = zeros(ntrls);
data.choiceIDs = {}; %zeros{ntrls};

%% PROCESS EACH TRIAL -- LOOP
for tt = 1:ntrls
	cur_trlinds = [1:trlbins]+trlbins*(tt-1);
	%pixelscaf = round(exptdata_mod{tt, 1}.m_aiStimulusArea/60);
	%stimloc(tt,:) = exptdata_mod{tt, 1}.m_aiStimulusRect;  
	trialstart_plx(tt) = exptdata_mod{tt, 1}.PlexonOnsetTime;

	if ~isfield(exptdata_mod{tt,1}, 'usebinary')
		exptdata_mod{tt, 1}.usebinary=0;
	end

	if ~isfield(exptdata_mod{tt,1}, 'UseLeye')
		exptdata_mod{tt, 1}.UseLeye=1;
		exptdata_mod{tt, 1}.UseReye=1;
	end
	useLeye(cur_trlinds)=exptdata_mod{tt,1}.UseLeye;
	useReye(cur_trlinds)=exptdata_mod{tt,1}.UseReye;

	%% Get task information
	data.task_outcome{tt} = exptdata_mod{tt, 1}.m_strctTrialOutcome.m_strResult; % self-explanatory
	if strcmp(exptdata_mod{tt, 1}.m_strctTrialOutcome.m_strResult, 'Incorrect')
		data.task_correct(tt) = -1; 
	elseif strcmp(exptdata_mod{tt, 1}.m_strctTrialOutcome.m_strResult, 'Correct')
		data.task_correct(tt) = 1;
	end

	data.task_type(tt) = exptdata_mod{tt, 1}.m_iMkTurkTaskType; % 1 means a color trial, 2 means a shape trial
	%data.task_catchtrial(tt) = exptdata_mod{tt, 1}.isCatchTrial; % 1 means the task type has just switched
	data.cueID(tt) = exptdata_mod{tt, 1}.m_iMkTurkTargetID; % the correct color/shape this trial
	data.choiceIDs{tt} = exptdata_mod{tt, 1}.m_iMkTurkChoiceIDs; % the possible color/shape IDs this trial (Inote that the correct one is always in the first index)

    
  %% Pull spike times/binned
	for channel = 1:nChans
		if ~isempty(exptdata_mod{tt, 10+2*(channel-1)})
			try
				cur_spks = exptdata_mod{tt,10+2*(channel-1)}.unit1;
				cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
				binned_SU1(cur_trlinds,channel) = histcounts(exptdata_mod{tt,10+2*(channel-1)}.unit1,edges_hist);
				spk_times_all{channel} = [spk_times_all{channel}; cur_spks+(trlsecs*(tt-1))];

			catch
				cur_spks = exptdata_mod{tt,10+2*(channel-1)};
				cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
				binned_SU1(cur_trlinds,channel) = histcounts(exptdata_mod{tt,10+2*(channel-1)},edges_hist);
				spk_times_all{channel,1} = [spk_times_all{channel,1}; cur_spks+(trlsecs*(tt-1))];
			end
		end   
%         try
%             binned_SU2(cur_trlinds,channel)=histcounts(exptdata_mod{tt,10+3*(channel-1)}.unit2,edges_hist);
% %         catch
% %             binned_SU2(cur_trlinds,targchans(channel))=0;
%         end 
%         try
%             binned_SU3(cur_trlinds,channel)=histcounts(exptdata_mod{tt,10+3*(channel-1)}.unit3,edges_hist);
%         end
%         try
%             binned_SU4(cur_trlinds,targchans(channel))=histcounts(exptdata_mod{tt,10+3*(channel-1)}.unit4,edges_hist);
%         end
    %binned_MUA(cur_trlinds,channel)=histc(exptdata_mod{tt,12+3*(channel-1)},edges);
    %LFP_trace(channel,cur_trlinds)=splitapply(@mean, exptdata_mod{tt,4}(channel,:), LFP_labels);
	end

%%
%     cur_trial_ET_trace_full=[]; cur_trial_ET_trace=[];
%     cur_trial_ETsamples=exptdata_mod{tt, 1}.m_afEyePositiontimes-exptdata_mod{tt, 1}.m_afEyePositiontimes(1)+g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(1);
%     cur_trial_ET_trace_full(:,1)=([exptdata_mod{tt,1}.m_afEyeXPositionScreenCoordinates-exptdata_mod{tt,1}.m_pt2iFixationSpot(1)])./pixelscaf;
%     cur_trial_ET_trace_full(:,2)=([exptdata_mod{tt,1}.m_afEyeYPositionScreenCoordinates-exptdata_mod{tt,1}.m_pt2iFixationSpot(2)])./pixelscaf;

	cur_trial_ET_trace_full=[]; cur_trial_ET_trace = [];
	%    cur_trial_ETsamples=[0:.001:3.999];
	cur_trial_ET_trace_full(:,1)=(exptdata_mod{tt, 7}.ET_trace(1,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
	cur_trial_ET_trace_full(:,2)=(exptdata_mod{tt, 7}.ET_trace(2,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
	
	if numETtraces > 2
    	cur_trial_ET_trace_full(:,3)=(exptdata_mod{tt, 7}.ET_trace(3,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
    	cur_trial_ET_trace_full(:,4)=(exptdata_mod{tt, 7}.ET_trace(4,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
	end
	if metadata_struct.ET_Eyelink == 3
		cur_trial_ET_trace_full(:,1)=(exptdata_mod{tt, 7}.ET_trace(3,:)' - median(exptdata_mod{tt, 7}.ET_trace(3,:)'))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
		cur_trial_ET_trace_full(:,2)=(exptdata_mod{tt, 7}.ET_trace(4,:)' - median(exptdata_mod{tt, 7}.ET_trace(4,:)'))*(g_strctEyeCalib.GainY.Buffer(end)./1000);    
	end 

	% % if needing to interpolate time positions
	%         cur_trial_ET_trace(:,1)=interp1(cur_trial_ETsamples,cur_trial_ET_trace_full(:,1),edges,'linear');
	%         cur_trial_ET_trace(:,2)=interp1(cur_trial_ETsamples,cur_trial_ET_trace_full(:,2),edges,'linear');
	
	%%

	if metadata_struct.ET_Eyelink >= 1
		eye_smooth = [15, 3, 3]; sgolay_deg = [2,2,2];
		sac_thresh = 9; %3.5; %threshold eye speed % default 6 for 0314
		peri_thresh = 3; %threshold eye speed for defining saccade boundary inds % default 2.5 for 0314  
	else
		cur_trial_ET_trace_full(:,1) = (exptdata_mod{tt, 7}.ET_trace(1,:)' - median(exptdata_mod{tt, 7}.ET_trace(1,:)'))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
		cur_trial_ET_trace_full(:,2) = (exptdata_mod{tt, 7}.ET_trace(2,:)' - median(exptdata_mod{tt, 7}.ET_trace(2,:)'))*(g_strctEyeCalib.GainY.Buffer(end)./1000);    
	
		eye_smooth = [181, 3, 5]; sgolay_deg = [3,2,4];
		sac_thresh = 2.6; %3.5; %threshold eye speed %2.5 with 1.2
		peri_thresh = 1.6; %threshold eye speed for defining saccade boundary inds
	end

	if plot_intermediate==1
		plot_flag=3;
	else
		plot_flag=0; 
	end

	% if any([strfind(exptname,'220203'), strfind(exptname,'220205'), strfind(exptname,'220207')])
	%     [~, ~, ~, saccade_times, sac_start_times, sac_stop_times] = detect_saccades_bevfel(...
	%         cur_trial_ET_trace_full', cur_trial_ETsamples, 3, eye_smooth, sgolay_deg, sac_thresh, peri_thresh, plot_flag);
	% else
	%     [~, ~, ~, saccade_times, sac_start_times, sac_stop_times] = detect_saccades_bevfel(...
	%         cur_trial_ET_trace_full', cur_trial_ETsamples, 3, eye_smooth, sgolay_deg, sac_thresh, peri_thresh, plot_flag);
	% end
    	
	%%
	if plot_intermediate==1
		pause
	end
			
	%% Processing of eye-tracking signals
	% do same processing of ET data as done for saccade detection
	if metadata_struct.ET_Eyelink == 3
		ET_trace_raw_1khz([1:trlsecs*1000]+trlsecs*1000*(tt-1),1:2)=cur_trial_ET_trace_full;
		ET_trace_raw_1khz([1:trlsecs*1000]+trlsecs*1000*(tt-1),3) = exptdata_mod{tt, 7}.ET_trace(1,:)';
	
		cur_trial_ETsamples_kofiko = exptdata_mod{tt, 1}.m_afEyePositiontimes - exptdata_mod{tt, 1}.m_afEyePositiontimes(1);
		%cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeXPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(1), linspace(0,4,240))';
		%cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeYPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(2), linspace(0,4,240))';    
		cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeXPositionScreenCoordinates'-960, linspace(0,4,240))';
		cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeYPositionScreenCoordinates'-540, linspace(0,4,240))';    
	else
		ET_trace_raw_1khz([1:trlsecs*1000]+trlsecs*1000*(tt-1),:)=cur_trial_ET_trace_full;
		% ET_ad_up(1,:) = smooth(cur_trial_ET_trace_full(:,1)', eye_smooth(1), 'sgolay', sgolay_deg(1));
		% ET_ad_up(2,:) = smooth(cur_trial_ET_trace_full(:,2)', eye_smooth(1), 'sgolay', sgolay_deg(1));
		% 
		% et_params.eye_fs = 60; %ET_adfreq;
		% ET_times_up = cur_trial_ETsamples(1):(1/et_params.eye_fs):cur_trial_ETsamples(end);
		% cur_trial_ET_trace(1,:) = interp1(cur_trial_ETsamples, ET_ad_up(1,:), ET_times_up, 'spline');
		% cur_trial_ET_trace(2,:) = interp1(cur_trial_ETsamples, ET_ad_up(2,:), ET_times_up, 'spline'); 
		% 
		% eye_smooth2=3;
		% cur_trial_ET_trace(1,:) = smooth(cur_trial_ET_trace(1,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
		% cur_trial_ET_trace(2,:) = smooth(cur_trial_ET_trace(2,:), eye_smooth(2), 'sgolay', sgolay_deg(2));		
	end

	% %%
%     cur_sac_start_inds=[]; cur_sac_stop_inds=[];
%     for sacs=1:length(saccade_times);
%         cur_sac_start_inds(sacs)=find(edges>sac_start_times(sacs),1,'first');
%         if isempty(find(edges>sac_stop_times(sacs),1,'first'))
%             cur_sac_stop_inds(sacs)=length(edges);
%         else
%         cur_sac_stop_inds(sacs)=find(edges>sac_stop_times(sacs),1,'first');
%         end
%     end

%     ET_trace_raw(cur_trlinds,1)=cur_trial_ET_trace_full(:,1);
%     ET_trace_raw(cur_trlinds,2)=cur_trial_ET_trace_full(:,2);
%    ET_trace(cur_trlinds,1)=smooth(ET_trace_raw(cur_trlinds,1),5);
%    ET_trace(cur_trlinds,2)=smooth(ET_trace_raw(cur_trlinds,2),5);

%    ET_trace_raw(cur_trlinds,1) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(:,1), ET_times_up, 'spline');
%    ET_trace_raw(cur_trlinds,2) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(:,2), ET_times_up, 'spline'); 
	
	ET_trace(cur_trlinds,1)=cur_trial_ET_trace(1,:)';
	ET_trace(cur_trlinds,2)=cur_trial_ET_trace(2,:)';

 %   sacc_inds=[sacc_inds; cur_trlinds([cur_sac_start_inds' cur_sac_stop_inds'])];
 %   bad_inds_sac=[bad_inds_sac, cur_trlinds(cur_sac_start_inds(sacs)):cur_trlinds(cur_sac_stop_inds(sacs))];

    %{
tt=3;
figure;
subplot(2,1,1); plot(exptdata_ColCloud{tt, 7}.ET_trace'); title('Plexon ET traces')
subplot(2,1,2); plot(exptdata_ColCloud{tt, 1}.ETthisTrialRawEyeData); title('Kofiko ET traces'); axis tight

    %}
% previous ET integration code    
%{
    try
    cur_trial_ETsamples=exptdata_mod{tt, 1}.m_afEyePositiontimes-exptdata_mod{tt, 1}.m_afEyePositiontimes(1)+g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(1);
    catch
%        cur_trial_ETsamples=2*[1:length(exptdata_mod{tt,1}.m_afEyeXPositionScreenCoordinates)]./length(exptdata_mod{tt,1}.m_afEyeXPositionScreenCoordinates);
        sprintf('ET data missing? \n')
    end
    cur_trial_ET_trace(:,1)=([interp1(cur_trial_ETsamples,exptdata_mod{tt,1}.m_afEyeXPositionScreenCoordinates,edges,'linear')-exptdata_mod{tt,1}.m_pt2iFixationSpot(1)])./pixelscaf;
    cur_trial_ET_trace(:,2)=([interp1(cur_trial_ETsamples,exptdata_mod{tt,1}.m_afEyeYPositionScreenCoordinates,edges,'linear')-exptdata_mod{tt,1}.m_pt2iFixationSpot(2)])./pixelscaf;


    ET_trace_raw(cur_trlinds,1)=cur_trial_ET_trace(:,1);
    ET_trace_raw(cur_trlinds,2)=cur_trial_ET_trace(:,2);

    ET_trace(cur_trlinds,1)=ET_trace_raw(cur_trlinds,1);
    ET_trace(cur_trlinds,2)=ET_trace_raw(cur_trlinds,2);
    
    if ~isempty(exptdata_mod{tt,6}.saccades);
%         cur_sac_inds=[1];
%         sacc_inds=[sacc_inds; [cur_trlinds(1)-1, cur_trlinds(1)]];
%         sac_edges=[0,1];
        cur_sac_inds=[]; cur_sac_start_inds=[]; cur_sac_stop_inds=[];
        
        if length(exptdata_mod{tt,6}.saccade_stop)<length(exptdata_mod{tt,6}.saccade_start)
            exptdata_mod{tt,6}.saccade_start=exptdata_mod{tt,6}.saccade_start(1:length(exptdata_mod{tt,6}.saccade_stop));
            exptdata_mod{tt,6}.saccades=exptdata_mod{tt,6}.saccades(1:length(exptdata_mod{tt,6}.saccade_stop));
        end
        if any(exptdata_mod{tt,6}.saccades-exptdata_mod{tt,2}>max(edges))
            exptdata_mod{tt,6}.saccades(exptdata_mod{tt,6}.saccades>max(edges))=[];
        end
        
        for sacs=1:length(exptdata_mod{tt,6}.saccades);
            curtrl_sactimes = exptdata_mod{tt,6}.saccades-exptdata_mod{tt,2};
            curtrl_sactimes_start = exptdata_mod{tt,6}.saccade_start-exptdata_mod{tt,2};
            curtrl_sactimes_stop = exptdata_mod{tt,6}.saccade_stop-exptdata_mod{tt,2};
            cur_sac_inds(sacs)=find(edges>curtrl_sactimes(sacs),1,'first');
            cur_sac_start_inds(sacs)=find(edges>curtrl_sactimes_start(sacs),1,'first');
            if isempty(find(edges>curtrl_sactimes_stop(sacs),1,'first'))
                cur_sac_stop_inds(sacs)=length(edges);
            else
            cur_sac_stop_inds(sacs)=find(edges>curtrl_sactimes_stop(sacs),1,'first');
            end
        end
            sacc_inds=[sacc_inds; cur_trlinds([cur_sac_start_inds' cur_sac_stop_inds'])];
            %sac_edges=[sac_edges;[cur_sac_start_inds' cur_sac_stop_inds']];
            bad_inds_sac=[bad_inds_sac, cur_trlinds(cur_sac_start_inds(sacs)):cur_trlinds(cur_sac_stop_inds(sacs))];
        
        %if cur_sac_inds(1)>1; sac_edges=[1,cur_sac_inds]; else; sac_edges=[cur_sac_inds];end
        %if cur_sac_stop_inds(end)<trlbins; sacc_inds=[sacc_inds; [trlbins-1 trlbins]]; sac_edges=[sac_edges,trlbins];end
    
        % detect saccades separately

    else
%        sacc_inds=[sacc_inds; [cur_trlinds(1), cur_trlinds(1); cur_trlinds(trlbins)-1, cur_trlinds(trlbins)]];
%        sac_edges=[0,1; trlbins, trlbins+1];
        bad_inds_sac=[bad_inds_sac cur_trlinds(1:6)];
    end
%}

%     for fixs=1:length(sac_edges)-1
%         cur_fixs_adj_inds=[sac_edges(fixs,2):sac_edges(fixs+1,1)];
%     %             ET_trace(cur_trlinds(cur_fixs_adj_inds),1)=mean(exptdata_mod{tt,3}(1,cur_fixs_adj_inds));
%     %             ET_trace(cur_trlinds(cur_fixs_adj_inds),2)=mean(exptdata_mod{tt,3}(2,cur_fixs_adj_inds));
%         ET_trace(cur_trlinds(cur_fixs_adj_inds),1)=mean(ET_trace(cur_trlinds(cur_fixs_adj_inds),1));
%         ET_trace(cur_trlinds(cur_fixs_adj_inds),2)=mean(ET_trace(cur_trlinds(cur_fixs_adj_inds),2));
%     end
    
   
	if tt >= 2
		Block_offsetinds = [Block_offsetinds,cur_trlinds(1)-1]; %end of last block
		Block_onsetinds = [Block_onsetinds,cur_trlinds(1)]; %start of this block block  
	end
%{    
    figure(1);
%     plot([ET_trace(cur_trlinds,1)-offsets(1)]./24,'r'); hold on
%     plot([ET_trace(cur_trlinds,2)-offsets(2)]./24,'b'); 
%     plot([ET_trace_raw(cur_trlinds,1)-offsets(1)]./24,'--r');
%     plot([ET_trace_raw(cur_trlinds,2)-offsets(2)]./24,'--b'); 
%     hold off
    plot([ET_trace(cur_trlinds,1)],'r'); hold on
    plot([ET_trace(cur_trlinds,2)],'b'); 
    plot([ET_trace_raw(cur_trlinds,1)],'--r');
    plot([ET_trace_raw(cur_trlinds,2)],'--b'); 
    if ~isempty(intersect(bad_inds_fix,cur_trlinds))
    vline(intersect(bad_inds_fix,cur_trlinds)-cur_trlinds(1),'b');
    end
    if ~isempty(intersect(bad_inds_artifact,cur_trlinds))
    vline(intersect(bad_inds_artifact,cur_trlinds)-cur_trlinds(1),'b');
    end
    title([num2str(var(ET_trace_raw(cur_trlinds,1))) '   ' num2str(var(ET_trace_raw(cur_trlinds,2)))])
    hold off    
%    vline(intersect(cur_trlinds,sac_inds)-cur_trlinds(1))
    vline(cur_sac_inds)
%    ylim([-10 10])
    pause
%}  
%/{
% remove to only extract ET info   
        
	if ~skipET
		switch exptdata_mod{tt, 1}.DualstimSecondaryUseCloud
			case 1
				data.stimET(cur_trlinds,:) = DualstimETbars(exptdata_mod{tt, 1}.stimseq_ET_bars(1:repframes:end),:);
				%stimulus(cur_trlinds_stim,:,:,:)=permute(exptdata_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
				stimETori(cur_trlinds,:) = exptdata_mod{tt, 1}.stimseq_ET_baroris(1:repframes:end);
				stimtypeET(cur_trlinds) = 1;
        
			case 7
				data.stimET(cur_trlinds,:,:,:) = permute( DensenoiseChromcloud_DKlspace(:,:,exptdata_mod{tt,1}.stimseq_ET_Cclouds(1:repframes:repframes*trlbins),:),[3 1 2 4]);
				%stimulus(cur_trlinds_stim,:,:,:)=permute(exptdata_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
				stimtypeET(cur_trlinds)=7;            
		end
	end
%}     

	if mod(tt,100)==0
		fprintf( '  Finished trial %4d of %4d\n', tt, ntrls )
	end

end

%%
tvec=[1:length(binned_SU1)];
Block_offsetinds=[Block_offsetinds,length(binned_SU1)];
%Block_offsetinds_blink=[Block_offsetinds_blink,length(binned_SU1)];
%Block_inds=[unique([Block_onsetinds,Block_onsetinds_blink]);unique([Block_offsetinds,Block_offsetinds_blink])];
Block_inds=[Block_onsetinds; Block_offsetinds];

bad_inds_fix=sort([Block_inds(1,:), Block_inds(1,:)+1, Block_inds(1,:)+2,Block_inds(1,:)+3,Block_inds(1,:)+4,Block_inds(1,:)+5,Block_inds(1,:)+6],1);
use_inds_fix=setdiff(tvec,unique(bad_inds_fix));
%{
use_inds_fix=setdiff(tvec,unique(bad_inds_fix));

use_inds_microsacc=setdiff(tvec,unique(bad_inds_sac));

use_inds_artifact=setdiff(tvec,unique(bad_inds_artifact));
use_inds_artifact(abs(ET_trace(use_inds_artifact,1))>25)=[];
use_inds_artifact(abs(ET_trace(use_inds_artifact,2))>25)=[];

use_inds=setdiff(tvec,unique([bad_inds_fix,bad_inds_sac,bad_inds_artifact]));
%}
%%

%% DAB: This stuff below was half-commented out and probably useless (%%%%)
%figure;
%subplot(3,1,1); imagesc(binned_SU1)
%subplot(3,1,2); imagesc(binned_SU2)
%subplot(3,1,3); imagesc(binned_SU3)
%figure; 
%imagesc(binned_SUSort)
%%
%%%% spksums=sum(binned_SU1);
%figure; histogram(spksums)
%%%% goodspkinds=find(spksums>3000);
%%
%%%% binned_SU=double(binned_SU1);

%{
%% quick STA test
cur_use_inds=find(stimtype>=2);
cur_use_inds(end-10:end)=[];
curlag=3;

figure;
stim2=reshape(stimulus,length(stimulus),3*25*25);
for cc=1:25;%[1 2 3 4 5 10 11 13 16:18 20 24];%
    cur_STA=binned_SU(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
    cur_StA2=reshape(cur_STA,3,25,25);
    curlim=max(abs(cur_StA2(:)));
    subplot(3,1,1)
    imagesc(squeeze(cur_StA2(1,:,:))); caxis([-curlim curlim]); title('Luminance axis'); pbaspect([1 1 1])
    subplot(3,1,2)
    imagesc(squeeze(cur_StA2(2,:,:))); caxis([-curlim curlim]); title('L-M axis'); pbaspect([1 1 1])
    subplot(3,1,3)
    imagesc(squeeze(cur_StA2(3,:,:))); caxis([-curlim curlim]); title('S axis'); pbaspect([1 1 1])
    figtitle(['unit ' num2str(cc) '   lag ' num2str(curlag)])
    colormap(gray)
    pause
end

trololo
%}
%{
%% quick STA to check (does not have eye correction)
cur_use_inds=intersect(find(stimtype==8),use_inds_fix);
cur_use_inds(end-10:end)=[];

stim2=(single(reshape(stim,size(stim,1),3*60*60))./127);
%stim3=reshape(stimET,size(stim,1),3*60*60);

%%
for cc=targchans
    tic
    figure;

    for curlag=1:6
        cur_STA1(curlag,:)=binned_SU(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
    end
    
    for curlag=1:6
    cur_StA2=reshape(cur_STA1(curlag,:),60,180);
    %curlim=150;%
    curlim=max(abs(cur_STA1(:)'));

    cur_StA2(:,[60 120])=curlim;
    subplot(6,1,curlag)
    imagesc(cur_StA2); clim([-curlim curlim]); pbaspect([3 1 1])

    colormap(gray); xlabel(['Lag ' num2str(curlag)]);
    if curlag==1; ylabel('S          L-M          Lum'); title(['Unit ' num2str(cc) ]); end
    end
%    figtitle(['Probe ' num2str(Robs_probe_ID(cc)) '   Unit ' num2str(cc) ])

    toc
%    pause
end  
%}
%% for bar stimulus STAs
%{
%binned_SU1=double(Robs');
binned_SU=double(binned_SU1);
stim3=double(stimET');
%stim3b=double(stimET(2:2:end,:));
%% STA across several lags
cur_use_inds=intersect(2:2:length(stim3),use_inds_fix);
%cur_use_inds=intersect(2:2:78420,valid_data);
cur_use_inds(end-10:end)=[];

for cc=1:134
    tic
    for curlag=1:10
%    cur_StA2=reshape(cur_STA{curlag},3,25,25);
%    cur_StA3=reshape(permute(cur_StA2,[3 2 1]),25,75);
%   cur_StA2=reshape(cur_STA{curlag},60,180);
%    cur_StA2(curlag,:)=double(binned_SU(cur_use_inds+curlag,cc))'*stim3(cur_use_inds,:);
    cur_StA2(curlag,:)=binned_SU(cur_use_inds+curlag,cc)'*stim3(cur_use_inds,:);
    
    %curlim=max(abs(cell2mat(cur_STA(:)')));
    curlim=max(abs(cur_StA2(:)'));
    imagesc(cur_StA2); caxis([-curlim curlim]); %pbaspect([1 3 1])
    colormap(gray); xlabel(['Lag ' num2str(curlag)]);
    if curlag==1; ylabel('S          L-M          Lum'); end
    end
    title(['SU Unit ' num2str(cc)])
    toc
    pause
end  
%%
%use_SUs=[1 2 4 10 13 17 20 24];
%binned_sua=binned_SU(:,use_SUs);
% good_SU=[1 3 5 6 7];
% ok_SU=[4 8:10 15];
binned_SU=binned_SU(:,iso_SUs);
figure; imagesc(binned_SU)

%}

%%


%% Now we package the data for later modeling
disp('Converting...') 
cd(output_directory)

%% Determine stimulus windows
%nofix = 1;  % what does this mean?
pixel_size = 1;

exptname = metadata_struct.exptname;
exptdate = metadata_struct.exptname(1:6);
curstimstype = 'mTurk';

cur_filename=['Jocamo_' exptname(1:6) '_' arraylabel '_mTurk_v09.mat'];
cur_filename_targchans=['Jocamo_' exptname(1:6) '_' arraylabel '_CC_ETCC_v09_cloudSUinds.mat'];

%dt=.016;
electrode_info = [];

% Determine stim locations and whether they move
% stim_location = exptdata_mod{end, 1}.m_aiStimulusRect;

if isfield(exptdata_mod{end,1}, 'm_aiTiledStimulusRect')
	num_stim_locs = size(exptdata_mod{end, 1}.m_aiTiledStimulusRect, 1);
else
	num_stim_locs = 1;
end


%if nofix
%	fixdot = [];
%else
%	overlap_d1 = any(fixlocs1 > stim_location(1) & fixlocs1 < stim_location(3));
%	dist_d1 = fix_location(1)-stim_location(1);   
%	overlap_d2=any(fixlocs2 > stim_location(2) & fixlocs2 < stim_location(4));
%	dist_d2=fix_location(2)-stim_location(2);

%	if overlap_d1 && overlap_d2

%		if stimscale > 1
			%todo: get this for new fixation calc
%			if mod(fix_size,2)
%				if mod(dist_d1,2)
%					fix_d1=ceil((dist_d1-fix_size-1)/stimscale):floor((dist_d1+fix_size-1)/stimscale);
%				else
%					fix_d1=ceil((dist_d1-fix_size)/stimscale):floor((dist_d1+fix_size)/stimscale);  
%				end

%				if mod(dist_d2,2)
%					fix_d2=ceil((dist_d2-fix_size-1)/stimscale):floor((dist_d2+fix_size-1)/stimscale);
%				else
%					fix_d2=ceil((dist_d2-fix_size)/stimscale):floor((dist_d2+fix_size)/stimscale);
%				end

%			else

%				if mod(dist_d1,2)
%					fix_d1=ceil((dist_d1-fix_size)/stimscale):floor((dist_d1+fix_size)/stimscale);
%				else
%					fix_d1=ceil((dist_d1-fix_size-1)/stimscale):floor((dist_d1+fix_size-1)/stimscale);
%				end

%				if mod(dist_d2,2)
%					fix_d2=ceil((dist_d2-fix_size)/stimscale):floor((dist_d2+fix_size)/stimscale);
%				else
%					fix_d2=ceil((dist_d2-fix_size-1)/stimscale):floor((dist_d2+fix_size-1)/stimscale);
%				end    
%			end
%			fixdot=[fix_d1;fix_d2];
%		else
			%         fix_d1=(dist_d1-fix_size):(dist_d1+fix_size);
			%         fix_d2=(dist_d2-fix_size):(dist_d2+fix_size);    
%			fix_d1=(fixlocs1-stim_location(1)); fix_d1=fix_d1(fix_d1>0); fix_d1=fix_d1(fix_d1<61);
%			fix_d2=(fixlocs2-stim_location(2)); fix_d2=fix_d2(fix_d2>0); fix_d2=fix_d2(fix_d2<61);
%			fixdot=[{fix_d2;fix_d1}];     
%		end
%		stim(fix_d2,fix_d1,1,:)=0;
%		stim(fix_d2,fix_d1,2,:)=0;
%		stim(fix_d2,fix_d1,3,:)=0;

%	else
%		fixdot=[0,0;0,0];
%	end
%end

%% Strange case of fixdotET
if ~skipET
	if targ_ETstimtype == 1
		data.stimET = data.stimET';
		stimETori = stimETori';
	else
		data.stimET = permute(data.stimET,[2 3 4 1]);
    
		%if nofix
		%	fixdotET=[];
		%else
		%	overlapET1_d1=any(fixlocs1>ETstim_location(1,1) & fixlocs1<ETstim_location(1,3));
		%	overlapET1_d2=any(fixlocs2>ETstim_location(1,2) & fixlocs2<ETstim_location(1,4));%
		%	overlapET2_d1=any(fixlocs1>ETstim_location(2,1) & fixlocs1<ETstim_location(2,3));
		%	overlapET2_d2=any(fixlocs2>ETstim_location(2,2) & fixlocs2<ETstim_location(2,4));

		%	if any([(overlapET1_d1 && overlapET1_d2), (overlapET2_d1 && overlapET2_d2)])
		%		fixET_d1=unique([fixlocs1-ETstim_location(1,1),fixlocs1-ETstim_location(2,1)]); fixET_d1=fixET_d1(fixET_d1>0); fixET_d1=fixET_d1(fixET_d1<61);
		%		fixET_d2=unique([fixlocs2-ETstim_location(1,2),fixlocs2-ETstim_location(2,2)]); fixET_d2=fixET_d2(fixET_d2>0); fixET_d2=fixET_d2(fixET_d2<61);
		%		fixdotET={fixET_d2;fixET_d1};  
		%		data.stimET(fixET_d2,fixET_d1,1,:) = 0;
		%		data.stimET(fixET_d2,fixET_d1,2,:) = 0;
		%		data.stimET(fixET_d2,fixET_d1,3,:) = 0;
        
		%	else
		%		fixdotET=[];
		%	end
		%end
	end
	stimtypeET = stimtypeET';
else
	curETstimtype = 'none';
end

%%
valid_data = use_inds_fix;
block_inds = Block_inds;

trial_start_ts = trialstart_plx';
trial_start_inds = round(trialstart_plx'.*1000);

sacc_inds=sacc_inds';
ETtrace_raw=ET_trace_raw_1khz';
ETtrace=ET_trace';
%ETtrace_plex_calib = PlexET_ad_calib';
ETgains = [g_strctEyeCalib.GainX.Buffer(:)' g_strctEyeCalib.GainY.Buffer(:)'];
useLeye=useLeye';
useReye=useReye';

targchans_SU = targchans(targchans <= nSU);
%targchans=1:nSU;

SU_clusters=[];
Robs = binned_SU1(:,targchans_SU)';
Robs_probe_ID = metadata_struct.spk_channels_SU(targchans_SU');
Robs_rating=[];  %spk_rating_SU(targchans_SU');

datafilts = ones(size(Robs));
%Robs_sort_inds=[1, 3, 4, 6, 7, 52, 67, 77, 81, 109, 112, 139, 152, 161, 163,171,176,188,190,191,194,194,205,214,215,217,218,222,224,225,226,235,236,240,241,244,249,252];
%targchansMU=[4 5 10 28 47 72 80 83 84 86 95 97 98 107 114 119]; %for 220209
%targchansMU=find(sum(binned_SU2)>500);
%RobsMU=binned_SU2(:,targchansMU)';

%targchansMU=[1:nMU]+nSU;
targchansMU = targchans(targchans > nSU);
RobsMU = binned_SU1(:,targchansMU)';
RobsMU_probe_ID = metadata_struct.spk_channels_MU(targchansMU'-nSU);
datafiltsMU = ones(size(RobsMU));

%% Accumulate spike-times and spikeIDs
spk_times=[]; spk_IDs=[];
for cc=1:length(spk_times_all)
	spk_times = [spk_times, spk_times_all{cc}'];
	spk_IDs = [spk_IDs, ones(1, length(spk_times_all{cc}))*cc];
end

%% Build data structure and save
data.exptname = exptname;
data.exptdate = exptdate;
data.dt = dt;
data.electrode_info = electrode_info;
%data.stim_location = stim_location;
%data.stim_location_deltas = delta_stimlocs;
%data.ETstim_location = ETstim_location;
data.pixel_size = pixel_size;

%data.fix_location = fix_location;
%data.fix_size = fix_size;
%data.fixdot = fixdot; % now handled in NDN directly
%data.stimloc = stimloc;  % this will be empty for some
%data.stimscale = stimscale;
%data.fixdotET = [];

data.ETtrace = ETtrace;
data.ETtrace_raw = ETtrace_raw;
data.ETgains = ETgains;
data.useLeye = useLeye;
data.useReye = useReye;
data.sacc_inds = sacc_inds;

%data.stim = stim;
data.stimtype = 0;
data.stimET = [];

data.trial_start_ts = trial_start_ts;
data.block_inds = block_inds;
data.valid_data = valid_data;
%data.blockID = blockID;  % only relevant to some conditions (?) but otherwise blank
%data.trialID = trialID;  % likewise

data.Robs = Robs;
data.RobsMU = RobsMU;
data.RobsMU_probe_ID = RobsMU_probe_ID;
data.Robs_probe_ID = Robs_probe_ID;
data.RobsMU_probe_ID = RobsMU_probe_ID;
data.Robs_rating = Robs_rating;
data.RobsMU_rating = []; % this is never have information  % RobsMU_rating = []; %spk_rating_MU(targchansMU'-nSU);
data.datafilts = datafilts;
data.datafiltsMU = datafiltsMU;
data.spike_ts = spk_times;
data.spikeIDs = spk_IDs;


% ET conditionals
if ~skipET
	data.stimET = stimET;
	data.stimtypeET = stimtypeET;
	%data.fixdotET = fixdotET;  % this will usually be blank
	if targ_ETstimtype == 1  % 1D-bar ET stimuli
		data.stimETori = stimETori;
	end
end

% SAVE
disp('Saving...') 
save( cur_filename, '-struct', 'data', '-v7.3')

%%%%%%%%%%% ORIGINAL SAVING FIASCO %%%%%%%%%%%%%%%%%%%%
%switch targ_stimtype
%    case 8
%        if ~skipET
%            if targ_ETstimtype==1;
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
%                'stim', 'stimET', 'stimETori', 'stimtype', 'stimtypeET', 'stimloc', 'stimscale','valid_data', 'cloud_scale', 'cloud_binary', 'block_inds', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%            else
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot','fixdotET',...
%                'stim', 'stimET', 'stimtype', 'stimtypeET','stimloc',  'stimscale','valid_data', 'cloud_scale', 'cloud_binary', 'block_inds', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%            end
%        else
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
%                'stim', 'stimtype', 'stimscale','stimloc', 'valid_data', 'block_inds', 'cloud_scale', 'cloud_binary', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%        end
%        save(cur_filename_targchans, 'targchans')
%    case 6
%        if ~skipET
%            if targ_ETstimtype==1;
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
%                'stim', 'stimET', 'stimETori', 'stimtype', 'hartley_metas', 'stimtypeET', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%            else
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
%                'stim', 'stimET', 'stimtype', 'hartley_metas', 'stimtypeET', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw',  'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%            end
%        else
%            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
%                'stim', 'stimtype', 'hartley_metas', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
%                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
%        end
%end

disp('Done saving output for modeling.') 

%% SAVE LFPs
if ~skipLFP
	disp('Saving LFPs...') 
	switch targ_stimtype
		case 8
			cur_filename_LFP=['Jocamo_' metadata_struct.exptname(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_v09_LFP.mat'];
			save(cur_filename_LFP,'LFP_ad', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
      
		case 6
			cur_filename_LFP=['Jocamo_' metadata_struct.exptname(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_v09_LFP.mat'];
			save(cur_filename_LFP, 'trial_start_ts', 'trial_start_inds', '-v7.3' )
	end
	disp('Done with LFPs') 
end

% cur_filename_spikes=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.mat'];
% save(cur_filename_spikes, 'spk_times_all');
% cur_filename_spikes_dat=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.dat'];
% writecell(spk_times_all, cur_filename_spikes_dat);
% cur_filename_spikes_py=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.npy'];
% writeNPY(spk_times_all, cur_filename_spikes_py);

