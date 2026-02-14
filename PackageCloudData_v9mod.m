function data = PackageCloudData_v9mod( exptdata, metadata_struct, targ_stimtype, targ_ETstimtype, stimFilePath, output_path, skipLFP, which_computer, LFP_ad )
%
% Usage: data = PackageCloudData_v9( exptdata, metadata_struct, <targ_stimtype>, <targ_ETstimtype>, <stimFilePath>, output_path, skipLFP, which_computer )
%
% Can check experiment composition
% experiment_composition( exptdata );

% For stimtype (should work with 6,7,8)
% 8 = Color cloud
% 7 = luminance cloud
% 6 = Color hartleys
% 0 = ground truth

% for ETstimtype: we use only 1, 7 (will work with both of these)
% 7 = Color cloud in both windows
% 1 = 1-D bars alternating each frame in both windows
% 2 = 1-D horizontal in 1, vertical in other
% 3 = Alternate to 2...

if (nargin < 3) || isempty(targ_stimtype)
	targ_stimtype = 8;    %this selects which stimulus to process: 8 looks for cloud stims
end
if (nargin < 3) || isempty(targ_ETstimtype)
	targ_ETstimtype = 0;    % this selects which stimulus to process: 8 looks for cloud stims
end
if nargin < 9
    skipLFP = 1;
end
%skipLFP = 0; %set to 1 if you want to skip laminar probe data and only analyzing the ET stims
%skipET = 1; %set to 1 if you want to skip the ETstim (such as if you moved them all the way out of the way); set to on on 5/12/23 on bevil's request to not have to worry about this
% Just set targ_ETstimtype to zero (default if you dont want this
if targ_ETstimtype == 0
	skipET = 1;
    curETstimtype = 'NA';
end
plot_intermediate=0;

%% reload data for analysis
if nargin < 8
	% This can be used to set default directories
	% Dan's laptop = 0
	% Bevil office desktop = 1
    % LSR 2A58 sorting rig = 2
	which_computer = 2; % default value
end

if (nargin < 7) || isempty(skipLFP)
	skipLFP = 0;
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
        case 2, stimFilePath = '/home/conwaylab/Processing/Cloudstims_calib_01_2022';
		otherwise
			disp('which_computer not defined')
	end
end

% Pull necessary variables from meta-data
exptname = metadata_struct.exptname;
arraylabel = metadata_struct.array_label;

if numel(arraylabel) > 1
    arraylabel_filepart = [cellfun(@(x) [x '_'], arraylabel(1:end-1), 'UniformOutput', false) arraylabel(end)];
    arraylabel_filepart = [arraylabel_filepart{:}];
else
    arraylabel_filepart = arraylabel;
end


useofflinesorting = metadata_struct.use_offline_sorting;
nSU = metadata_struct.nSU;
nMU = metadata_struct.nMU;
g_astrctAllParadigms = metadata_struct.g_astrctAllParadigms;  % expt configuration information
g_strctEyeCalib = metadata_struct.g_strctEyeCalib;
% output_directory = [output_path 'Analysis' filesep];
output_directory = [output_path];

% Set up data-struct for output
data.exptname = exptname;
data.stim = [];   % these are the big variables that not worth representing twice
data.stimET = [];

%% params
%if useofflinesorting == 1
	nChans=metadata_struct.nSU+metadata_struct.nMU;
	spk_ID=[metadata_struct.spk_ID_SU; metadata_struct.spk_ID_MU];
	spk_ch=[metadata_struct.spk_channels_SU+1; metadata_struct.spk_channels_MU+1];
% else
% 	nChans=24;
% end

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
	numSpks(cc)=length(cell2mat(exptdata(:, 10+(cc-1))));
end
%targchans=find(numSpks>2000);

if useofflinesorting == 1
	targchans=1:nChans;  % pass all channels through
end

%% Assign list of trials to teach type of data from the experiment
trlonset_diffs = [4; diff(cell2mat(exptdata(:, 2)))];

trialdur = metadata_struct.trialdur;

for switch_stimtype=0:8

	targ_trials=[];  
	for tt=1:length(exptdata) % tt = number of trials
		%    if strcmp(exptdata{tt, 1}.m_strTrialType, 'Dense Noise');
        if ~isfield(exptdata{tt,1},'m_bMonkeyFixatedOverride'); exptdata{tt, 1}.m_bMonkeyFixatedOverride=0; end
		if strcmp(exptdata{tt, 1}.m_strTrialType, 'Dual Stim') && ...
			(exptdata{tt, 1}.m_bMonkeyFixated == 1 | exptdata{tt, 1}.m_bMonkeyFixatedOverride == 1) && ...
            (trlonset_diffs(tt) > trialdur) && ...
			(exptdata{tt, 1}.DualstimPrimaryuseRGBCloud == switch_stimtype) &&...
            (exptdata{tt, 1}.m_aiStimulusArea > 0)
			% exptdata{tt, 1}.DualstimSecondaryUseCloud==targ_ETstimtype; % && exptdata{tt, 1}.m_aiStimulusRect(1)==975;
			targ_trials=[targ_trials,tt];
		end
	end

	switch switch_stimtype
		case 0; curstimstype='GT';  exptdata_GT=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_GT','-v7.3')
		case 3; curstimstype='HL';  exptdata_HartleyLum=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_HartleyLum','-v7.3')
		case 6; curstimstype='HC';  exptdata_HartleyCol=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_HartleyCol','-v7.3')
		case 8; curstimstype='CC';    exptdata_ColCloud=exptdata(targ_trials,:);
			%save([output_directory,filesep,'exptdata4modeling' curstimstype '.mat'],'exptdata_ColCloud','-v7.3')
	end

end

try
    DualstimETbars = int8(squeeze((metadata_struct.g_astrctAllParadigms.DualstimETbars-128)/127)');
catch
    warning('No ET bars found - so there is no no 1D noise data. not a problem for experiments using all clouds.')
end
%% Stimulus selection
switch targ_stimtype
	case 3; curstimstype='HL';	exptdata_mod = exptdata_HartleyLum;
	case 6; curstimstype='HC';	exptdata_mod = exptdata_HartleyCol;
	case 7; curstimstype='LC';	exptdata_mod = exptdata_LumCloud;
	case 8; curstimstype='CC';	exptdata_mod = exptdata_ColCloud;
end
switch targ_ETstimtype
	case 1; curETstimtype='1D';
	case 7; curETstimtype='CC';
end

%% Random assortment of setup things -- before going into trial-by-trial processing
%ntrls = min([size(exptdata_mod,1),500]); % only using the first 500
ntrls = size(exptdata_mod,1); % find number of trials
NT=trlbins*ntrls;

if any([strfind(exptname,'220203'), strfind(exptname,'220205'), strfind(exptname,'220207')])
	trlsecs=2.67; % number of seconds in each trial
	edges = linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);
	cur_trial_ETsamples = [0:.001:2.669];
else
	trlsecs = 4; 
	edges = linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);  
	cur_trial_ETsamples = [0:.001:3.999];
end

if metadata_struct.ET_Eyelink == 3
	numETtraces = 2; % number of eye traces
	ET_trace_raw_1khz = zeros(ntrls*trlsecs*1000, 3);
elseif metadata_struct.ET_Eyelink == 4 %%% Added by [Ramon Bartolo, 20250522] to handle binocular dDPI
    numETtraces = 4; %number of eye traces for binocular dDPI
	ET_trace_raw_1khz = zeros(ntrls*trlsecs*1000, numETtraces);
else
	numETtraces = size(exptdata_mod{end, 7}.ET_trace, 1);
	ET_trace_raw_1khz = zeros(ntrls*trlsecs*1000, numETtraces);
end
%ET_ad_up=[];

data.stim = int8(zeros(NT,60,60,3)); 
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
cloud_area = zeros(NT,1);  

%% Load Hartleys if necessary
if (targ_stimtype == 3) || (targ_stimtype == 6)
	%load('C:\SpkSort2023\Cloudstims_calib_01_2022\hartleys_60.mat')
	load(sprintf('%shartleys_60.mat', stimFilePath))
	hartleys = hartleys60_DKL;
	hartleys_metas = hartleys60_meta;
end


cd(stimFilePath)
%cur_scale=g_astrctAllParadigms{1, 1}.DualstimScale.Buffer(find(g_astrctAllParadigms{1, 1}.DualstimScale.TimeStamp<exptdata_mod{2,2},1,'last'));
try
    cur_scale=g_astrctAllParadigms.DualstimScale.Buffer(find(g_astrctAllParadigms.DualstimScale.TimeStamp<exptdata_mod{2,2},1,'last'));
catch
    warning('Something is wrong with stim scale - how long did this experiment run?')
    cur_scale=g_astrctAllParadigms.m_fInitial_DualstimPrimaryCloudScale;
end
load(sprintf('Cloudstims_Chrom_size60_scale%d_%02d.mat', cur_scale, cur_BlockID));
%load(sprintf([stimFilePath 'Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID)) 
%DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));

%% Load external eye-tracking file info -- has two variables (at 1 kHZ)
% NOT NEEDED -- plexon information already went into trial-level information
% load(sprintf('%s%s_FullExpt_ET.mat', metadata_struct.expt_folder, exptname ))
% variables: PlexET_ad_calib, PlexET_times 


%% PROCESS EACH TRIAL -- LOOP
for tt = 1:ntrls
	cur_trlinds = [1:trlbins]+trlbins*(tt-1);
	pixelscaf = round(exptdata_mod{tt, 1}.m_aiStimulusArea/60);
	%stimloc(tt,:) = exptdata_mod{tt, 1}.m_aiStimulusRect;  
	trialstart_plx(tt) = exptdata_mod{tt, 1}.PlexonOnsetTime;

	if ~isfield(exptdata_mod{tt,1}, 'usebinary')
		exptdata_mod{tt, 1}.usebinary=0;
	end

	if targ_stimtype == 8
		if exptdata_mod{tt,1}.BlockID ~= cur_BlockID
			cur_BlockID = exptdata_mod{tt,1}.BlockID;
			%cur_scale=g_astrctAllParadigms{1, 1}.DualstimScale.Buffer(find(g_astrctAllParadigms{1, 1}.DualstimScale.TimeStamp<exptdata_mod{tt,2},1,'last'));
            try
    			cur_scale=g_astrctAllParadigms.DualstimScale.Buffer(find(g_astrctAllParadigms.DualstimScale.TimeStamp<exptdata_mod{tt,2},1,'last'));
            catch
                cur_scale=g_astrctAllParadigms.m_fInitial_DualstimPrimaryCloudScale; %if Kofiko crashed and this stuff wasn't saved, we'll just have to hope that nobody changed the default value here.
            end

			if exptdata_mod{tt, 1}.usebinary==1
				load(sprintf(['Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat'], cur_scale, cur_BlockID));            
            elseif exptdata_mod{tt, 1}.usebinary==2
				load(sprintf(['Cloudstims_ContrastMatched_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID));      
            else
				load(sprintf(['Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID));
			end
			DensenoiseChromcloud_DKlspace=int8(10*127*(DensenoiseChromcloud_DKlspace));
		end
    
		cur_TrialID = exptdata_mod{tt,1}.TrialID;
		BlockID(cur_trlinds)=cur_BlockID;
		TrialID(cur_trlinds)=cur_TrialID;
    
		cloud_scale(cur_trlinds) = cur_scale;
		cloud_binary(cur_trlinds) = exptdata_mod{tt, 1}.usebinary;  
        cloud_area(cur_trlinds) = exptdata_mod{tt,1}.m_aiStimulusArea;
	end

	if ~isfield(exptdata_mod{tt,1}, 'UseLeye')
		exptdata_mod{tt, 1}.UseLeye=1;
		exptdata_mod{tt, 1}.UseReye=1;
	end
	useLeye(cur_trlinds)=exptdata_mod{tt,1}.UseLeye;
	useReye(cur_trlinds)=exptdata_mod{tt,1}.UseReye;
    
  % Pull spike times/binned
	for channel = 1:nChans
		if ~isempty(exptdata_mod{tt, 10+(channel-1)})
			try
				cur_spks = exptdata_mod{tt,10+(channel-1)}.unit1;
				cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
				binned_SU1(cur_trlinds,channel) = histcounts(exptdata_mod{tt,10+(channel-1)}.unit1,edges_hist);
				spk_times_all{channel} = [spk_times_all{channel}; cur_spks+(trlsecs*(tt-1))];

			catch
				cur_spks = exptdata_mod{tt,10+(channel-1)};
				cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
				binned_SU1(cur_trlinds,channel) = histcounts(exptdata_mod{tt,10+(channel-1)},edges_hist);
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
    
    %%% Added by [Ramon Bartolo, 20250522] to handle binocular dDPI
    if metadata_struct.ET_Eyelink == 4 && size(exptdata_mod{tt, 7}.ET_trace,1)==8
        %%% on exptdata_mod.ET_trace from binocular dDPI data:
        %%% row 1: strobe
        %%% row 2: null (strobe artifact)
        %%% row 3: L pupil size
        %%% row 4: R pupil size
        %%% rows 5-6: XY (R)
        %%% rows 7-8: XY (L)
        cur_trial_ET_trace_full=[]; cur_trial_ET_trace = [];
        cur_trial_ET_trace_full(:,1)=(exptdata_mod{tt, 7}.ET_trace(5,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000); %X right eye
	    cur_trial_ET_trace_full(:,2)=(exptdata_mod{tt, 7}.ET_trace(6,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000); %Y right eye
		cur_trial_ET_trace_full(:,3)=(exptdata_mod{tt, 7}.ET_trace(7,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000); %X left eye
    	cur_trial_ET_trace_full(:,4)=(exptdata_mod{tt, 7}.ET_trace(8,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000); %Y left eye
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
		cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeXPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(1), linspace(0,4,240))';
		cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeYPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(2), linspace(0,4,240))';

    elseif metadata_struct.ET_Eyelink == 4 %%% Added by [Ramon Bartolo, 20250522] to handle binocular dDPI
		ET_trace_raw_1khz([1:trlsecs*1000]+trlsecs*1000*(tt-1),1:4)=cur_trial_ET_trace_full;
		ET_trace_raw_1khz([1:trlsecs*1000]+trlsecs*1000*(tt-1),5:6) = exptdata_mod{tt, 7}.ET_trace(3:4,:)'; %adding pupil sizes
	    
		cur_trial_ETsamples_kofiko = exptdata_mod{tt, 1}.m_afEyePositiontimes - exptdata_mod{tt, 1}.m_afEyePositiontimes(1);
		cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeXPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(1), linspace(0,4,240))';
		cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, exptdata_mod{tt, 1}.m_afEyeYPositionScreenCoordinates' - exptdata_mod{tt, 1}.m_pt2iFixationSpot(2), linspace(0,4,240))';

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
% comment out this section if you want to only extract ET info   
    
	switch exptdata_mod{tt,1}.DualstimPrimaryuseRGBCloud
		case 8
			data.stim(cur_trlinds,:,:,:) = permute( DensenoiseChromcloud_DKlspace(:,:,exptdata_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:),[3 1 2 4]);
			%stimulus(cur_trlinds_stim,:,:,:)=permute(exptdata_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
			stimtype(cur_trlinds)=8;

		case 7
			data.stim(cur_trlinds,:,:,:) = permute( DensenoiseAchromcloud_binned(:,:,exptdata_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:),[3 1 2 4]);
			%stimulus(cur_trlinds_stim,:,:,:)=permute(exptdata_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
			stimtype(cur_trlinds) = 8;

		case 6    
			data.stim(cur_trlinds,:,:,:) = hartleys60_DKL(exptdata_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:,:,:);
			stimtype(cur_trlinds)=6;
			hartleystim_metas(cur_trlinds,:) = hartleys_metas(exptdata_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:);

		case 3
			data.stim(cur_trlinds_stim,:,:,:) = hartleys60_DKL(exptdata_mod{tt,1}.stimseq(1:repframes:end),:,:,:);
			stimtype(cur_trlinds) = 3;
			hartleystim_metas(cur_trlinds_stim,:) = hartleys_metas(exptdata_mod{tt,1}.stimseq(1:2:end),:);
	end
    

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

bad_inds_block=sort([Block_inds(1,:), Block_inds(1,:)+1, Block_inds(1,:)+2,Block_inds(1,:)+3,Block_inds(1,:)+4,Block_inds(1,:)+5,Block_inds(1,:)+6],1);

ETdist_thresh=40; 
bad_inds_fix = unique([find(abs(ET_trace(:,1))>ETdist_thresh);find(abs(ET_trace(:,2))>ETdist_thresh)])';
bad_inds_all = unique([bad_inds_block,bad_inds_fix,bad_inds_fix-1, bad_inds_fix-2, bad_inds_fix+1, bad_inds_fix+2]); % remove indices immediately preceding and following eye movement artifacts
use_inds_fix=setdiff(tvec,bad_inds_all);

% find sequences less than 10 due to eye movement removal, and exclude them
% to avoid clogging up the modeling pipeline with tiny snippets
diffs=diff([1,use_inds_fix]);
[~,X]=find(diff(diffs)<10);  
for k= X, use_inds_fix(diffs(k):diffs(k+1)-1)=0; end
use_inds_fix(use_inds_fix==0)=[];

%{
use_inds_fix=setdiff(tvec,unique(bad_inds_fix));

use_inds_microsacc=setdiff(tvec,unique(bad_inds_sac));

use_inds_artifact=setdiff(tvec,unique(bad_inds_artifact));
use_inds_artifact(abs(ET_trace(use_inds_artifact,1))>25)=[];
use_inds_artifact(abs(ET_trace(use_inds_artifact,2))>25)=[];

use_inds=setdiff(tvec,unique([bad_inds_fix,bad_inds_sac,bad_inds_artifact]));
%}

%%


%% Now we package the data for later modeling
disp('Converting...') 
try
    cd(output_directory)
catch
    mkdir(output_directory)
    cd(output_directory)
end
%% Determine stimulus windows
%nofix = 1;  % what does this mean?
pixel_size = 1;

exptname = metadata_struct.exptname;
exptdate = metadata_struct.exptname(1:6);

cur_filename=[metadata_struct.monkey_name '_' exptname(1:6) '_' arraylabel_filepart '_' curstimstype '_ET' curETstimtype '_v09.mat'];
cur_filename_targchans=[metadata_struct.monkey_name '_' exptname(1:6) '_' arraylabel_filepart '_CC_ETCC_v09_cloudSUinds.mat'];

%dt=.016;
electrode_info = [];

% Determine stim locations and whether they move
% stim_location = exptdata_mod{end, 1}.m_aiStimulusRect;

if isfield(exptdata_mod{end,1}, 'm_aiTiledStimulusRect')
	num_stim_locs = size(exptdata_mod{end, 1}.m_aiTiledStimulusRect, 1);
else
	num_stim_locs = 1;
end
stim_locs = zeros(size(exptdata_mod,1), num_stim_locs, 4);
stim_area = zeros(size(exptdata_mod,1), 1);

for nn = 1:size(exptdata_mod,1)
	if isfield(exptdata_mod{end,1}, 'm_aiTiledStimulusRect')
		stim_locs(nn, :, :) = exptdata_mod{nn, 1}.m_aiTiledStimulusRect;
	else
		stim_locs(nn, :, :) = exptdata_mod{nn, 1}.m_aiStimulusRect;
    end
    stim_area(nn) = exptdata_mod{nn, 1}.m_aiStimulusArea;
end

if unique(stim_area)>1;
    disp('Caution: multiple stimulus areas detected. Check stim_area for consistency!')
end

L = mode(stim_area); % identify the stim area with the most trials. hopefully this will always be 60. if not, we can code in targeting for multiple stim sizes.

% try % old way to calculate stim area introduced 
%     L = stim_locs(1,3,1)-stim_locs(1,1,1);
% catch
%     L = stim_locs(1,1,3)-stim_locs(1,1,1);
% end

% Determine median stim location and use that -- but then record shifts relative to that
tcx = median(stim_locs(:,:,1));
tcy = median(stim_locs(:,:,2));
stim_location = [tcx' tcy' tcx'+L tcy'+L];
delta_stimlocs = [stim_locs(:,1,1)-tcx(1) stim_locs(:,1,2)-tcy(1)];
if sum(abs(delta_stimlocs), "all") > 0
	disp('  Stimulus window shifts detected: be sure to use stim_location_deltas')
end

ETstim_location = [exptdata_mod{end, 1}.secondarystim_bar_rect;   
exptdata_mod{end, 1}.tertiarystim_bar_rect];
fix_location = exptdata_mod{end, 1}.m_pt2iFixationSpot;
fix_size = exptdata_mod{end, 1}.m_fFixationSizePix-1;

data.stim = permute(data.stim, [2 3 4 1]);
stimtype = stimtype';
stimscale = (stim_location(3)-stim_location(1))/60;

%% Fixation dots
fixlocs1 = [fix_location(1)-fix_size:fix_location(1)+fix_size];
fixlocs2 = [fix_location(2)-fix_size:fix_location(2)+fix_size];

%{
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
%}
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
try
    ETgains = [g_strctEyeCalib.GainX.Buffer, g_strctEyeCalib.GainY.Buffer];
catch
    ETgains = [g_strctEyeCalib.GainX.Buffer(end), g_strctEyeCalib.GainY.Buffer(end)]; %kludge, but should work for msot cases: usually we set the gains at the beginning and then leave them there
end
useLeye=useLeye';
useReye=useReye';

switch targ_stimtype
  
	case 8
		blockID = BlockID;
		trialID = TrialID;
		cloud_scale = cloud_scale';
		cloud_binary = cloud_binary';
        cloud_area = cloud_area';
%		targchans=find(sum(binned_SU1)>2000);
    
	case 6
		blockID = BlockID;
		trialID = TrialID;
        hartley_metas = hartleystim_metas';
%		load(cur_filename_targchans)
end

%targchans_SU = targchans(targchans <= nSU);
targchans_SU=1:nSU;

SU_clusters=[];
Robs = binned_SU1(:,targchans_SU)';
Robs_probe_ID = metadata_struct.spk_channels_SU(targchans_SU');
Robs_rating=[];  %spk_rating_SU(targchans_SU');

datafilts = ones(size(Robs));
%Robs_sort_inds=[1, 3, 4, 6, 7, 52, 67, 77, 81, 109, 112, 139, 152, 161, 163,171,176,188,190,191,194,194,205,214,215,217,218,222,224,225,226,235,236,240,241,244,249,252];
%targchansMU=[4 5 10 28 47 72 80 83 84 86 95 97 98 107 114 119]; %for 220209
%targchansMU=find(sum(binned_SU2)>500);
%RobsMU=binned_SU2(:,targchansMU)';

%targchansMU = targchans(targchans > nSU);
targchansMU=[1:nMU]+nSU;
RobsMU = binned_SU1(:,targchansMU)';
RobsMU_probe_ID = metadata_struct.spk_channels_MU(targchansMU'-nSU);
datafiltsMU = ones(size(RobsMU));

%% Accumulate spike-times and KEEP original Phy/KS cluster IDs
% commented out by MJG 2/13/26
% spk_times = [];
% spk_IDs   = [];
% for cc = 1:length(spk_times_all)
%     if ~isempty(spk_times_all{cc})
%         spk_times = [spk_times, spk_times_all{cc}'];
%         spk_IDs   = [spk_IDs, repmat(spk_ID(cc), 1, numel(spk_times_all{cc}))];  % use true cluster ID
%     end
% end
% data.spike_ts = spk_times;
% data.spikeIDs = spk_IDs;   % now these match Phy/Kilosort cluster IDs


spk_times=[]; spk_IDs=[];
for cc=1:length(spk_times_all)
	spk_times = [spk_times,spk_times_all{cc}'];
	spk_IDs = [spk_IDs, ones(1,length(spk_times_all{cc}))*cc];
end
data.spike_ts = spk_times;
data.spikeIDs = spk_IDs; 

%% Build data structure and save
data.exptname = exptname;
data.exptdate = exptdate;
data.dt = dt;
data.electrode_info = electrode_info;
data.stim_location = stim_location;
data.stim_location_deltas = delta_stimlocs;
data.stim_area = stim_area;
data.ETstim_location = ETstim_location;
data.pixel_size = pixel_size;

data.fix_location = fix_location;
data.fix_size = fix_size;
%data.fixdot = fixdot; % now handled in NDN directly
%data.stimloc = stimloc;  % this will be empty for some
data.stimscale = stimscale;
%data.fixdotET = [];

data.ETtrace = ETtrace;
data.ETtrace_raw = ETtrace_raw;
data.ETgains = ETgains;
data.useLeye = useLeye;
data.useReye = useReye;
data.sacc_inds = sacc_inds;

data.trial_start_ts = trial_start_ts;
data.block_inds = block_inds;
data.valid_data = valid_data;
data.blockID = blockID;  % only relevant to some conditions (?) but otherwise blank
data.trialID = trialID;  % likewise

data.stimtype = stimtype;
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

% Add hartley-specific or cloud-specific information 
switch targ_stimtype
	case 8
		data.cloud_scale = cloud_scale;
		data.cloud_binary = cloud_binary;
        data.cloud_area = cloud_area;
	case 6
		data.hartley_metas = hartley_metas;
end

% ET conditionals
if ~skipET
	%data.stimET = stimET;
	data.stimtypeET = stimtypeET;
	%data.fixdotET = fixdotET;  % this will usually be blank
	if targ_ETstimtype == 1  % 1D-bar ET stimuli
		data.stimETori = stimETori;
	end
end

% SAVE
disp('Saving...') 
save( fullfile(output_path, cur_filename), '-struct', 'data', '-v7.3')

disp('Done saving output for modeling.') 

%% SAVE LFPs
if ~skipLFP
	disp('Saving LFPs...') 
	switch targ_stimtype
		case 8
			cur_filename_LFP=[metadata_struct.monkey_name '_' metadata_struct.exptname(1:6) '_' arraylabel_filepart '_' curstimstype '_ET' curETstimtype '_v09_LFP.mat'];
			try
            save(fullfile(output_path,cur_filename_LFP),'LFP_ad', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
            catch
                keyboard;
            end
		case 6
			cur_filename_LFP=[metadata_struct.monkey_name '_' metadata_struct.exptname(1:6) '_' arraylabel_filepart '_' curstimstype '_ET' curETstimtype '_v09_LFP.mat'];
			save(fullfile(output_path, cur_filename_LFP), 'trial_start_ts', 'trial_start_inds', '-v7.3' )
	end
	disp('Done with LFPs') 
end

% cur_filename_spikes=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.mat'];
% save(cur_filename_spikes, 'spk_times_all');
% cur_filename_spikes_dat=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.dat'];
% writecell(spk_times_all, cur_filename_spikes_dat);
% cur_filename_spikes_py=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.npy'];
% writeNPY(spk_times_all, cur_filename_spikes_py);

