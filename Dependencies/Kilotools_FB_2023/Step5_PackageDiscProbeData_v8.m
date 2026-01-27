%%
%load('/media/felix/Internal_1/Data/BevilColor/220727_134235_Jacomo/Analysis/FullDiscProbeExperimentRecordingwLFP_ks1_full_v08.mat')
%load('/media/felix/Internal_1/Data/BevilColor/220727_134235_Jacomo/Analysis/Jocamo_220727_full_CC_ETCC_nofix_v08_cloudSUinds.mat', 'targchans')
targchans
%%
%ET_Eyelink=1;
nChans=(size(ExperimentRecording,2)-8)/3;
%%
ntrls = size(ExperimentRecording,1);
trlbins=13;
NT=trlbins*ntrls;

trlsecs = 0.2; trlsecs_ms = trlsecs*1000;
edges=linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);

plot_intermediate=0;

if ET_Eyelink==3
    numETtraces=2;
    ET_trace_raw_1khz=zeros(ntrls*trlsecs*1000,3);
else
    numETtraces=size(ExperimentRecording{end, 7}.ET_trace,1);
    ET_trace_raw_1khz=zeros(ntrls*trlsecs*1000,numETtraces);
end
ET_ad_up=[];
ET_trialsec=101:300; % in case ET traces include pre/post trial period

% if any([strfind(filenameP,'220203'), strfind(filenameP,'220205'), strfind(filenameP,'220207')])
% edges=linspace(0,2.67,trlbins); edges_hist=linspace(0,2.67,trlbins+1);
% else
% edges=linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);
% end

stim=int8(zeros(NT,3)); 
stim_polar=double(zeros(NT,3)); 

trialstart_plx=zeros(ntrls,1);

binned_SU1=int8(zeros(NT,nChans)); binned_MUA=[]; ET_trace=zeros(NT,2); ET_trace_raw=zeros(NT,2);
%binned_SU2=int8(zeros(NT,nChans)); binned_SU3=int8(zeros(NT,nChans)); %binned_SU4=zeros(NT,nChans); 
bad_inds_fix=[1:10]; bad_inds_sac = []; sacc_inds=[];
Block_onsetinds=1; Block_offsetinds=[]; Block_onsetinds_blink=1; Block_offsetinds_blink=[];
BlockID=[]; TrialID=[]; useLeye=zeros(NT,1); useReye=zeros(NT,1);
spk_times_all=cell([nChans,1]);

%%
for tt=1:ntrls
    cur_trlinds=[1:trlbins]+trlbins*(tt-1);
    trialstart_plx(tt) = ExperimentRecording{tt, 1}.PlexonOnsetTime;

    % commented out: eye information

%/{ 
% comment out to only get ET info    
    for channel=1:nChans
        if ~isempty(ExperimentRecording{tt,10+3*(channel-1)})
            try
                cur_spks=ExperimentRecording{tt,10+3*(channel-1)}.unit1;
                cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
                binned_SU1(cur_trlinds,channel)=histcounts(ExperimentRecording{tt,10+3*(channel-1)}.unit1,edges_hist);
                spk_times_all{channel}=[spk_times_all{channel}; cur_spks+(trlsecs*(tt-1))];
            catch
                cur_spks=ExperimentRecording{tt,10+3*(channel-1)};
                cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
                binned_SU1(cur_trlinds,channel)=histcounts(ExperimentRecording{tt,10+3*(channel-1)},edges_hist);
                spk_times_all{channel,1}=[spk_times_all{channel,1}; cur_spks+(trlsecs*(tt-1))];
            end
        end   
    end
%}

%%
    cur_trial_ET_trace_full=[]; cur_trial_ET_trace=[];
%    cur_trial_ETsamples=[0:.001:0.199];
    cur_trial_ETsamples=[g_strctStatistics.preTrialWindow:.001:g_strctStatistics.postTrialWindow-.001];
    cur_trial_ET_trace_full(:,1)=(ExperimentRecording{tt, 7}.ET_trace(1,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
    cur_trial_ET_trace_full(:,2)=(ExperimentRecording{tt, 7}.ET_trace(2,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
if numETtraces>2
    cur_trial_ET_trace_full(:,3)=(ExperimentRecording{tt, 7}.ET_trace(3,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
    cur_trial_ET_trace_full(:,4)=(ExperimentRecording{tt, 7}.ET_trace(4,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
end
if ET_Eyelink==3
    cur_trial_ET_trace_full(:,1)=(ExperimentRecording{tt, 7}.ET_trace(3,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
    cur_trial_ET_trace_full(:,2)=(ExperimentRecording{tt, 7}.ET_trace(4,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
end 
ET_trace_raw_1khz([1:trlsecs_ms]+trlsecs_ms*(tt-1),:)=cur_trial_ET_trace_full(ET_trialsec,:);

%%
    if ET_Eyelink>=1
        eye_smooth = [15, 3, 3]; sgolay_deg = [2,2,2];
        sac_thresh = 9; %3.5; %threshold eye speed % default 6 for 0314
        peri_thresh = 3; %threshold eye speed for defining saccade boundary inds % default 2.5 for 0314
    else
         cur_trial_ET_trace_full(:,1)=(ExperimentRecording{tt, 7}.ET_trace(1,:)' - median(ExperimentRecording{tt, 7}.ET_trace(1,:)'))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
         cur_trial_ET_trace_full(:,2)=(ExperimentRecording{tt, 7}.ET_trace(2,:)' - median(ExperimentRecording{tt, 7}.ET_trace(2,:)'))*(g_strctEyeCalib.GainY.Buffer(end)./1000);    
    
         eye_smooth = [181, 3, 5]; sgolay_deg = [3,2,4];
         sac_thresh = 2.6; %3.5; %threshold eye speed %2.5 with 1.2
         peri_thresh = 1.6; %threshold eye speed for defining saccade boundary inds
    end

    if plot_intermediate==1
        plot_flag=3;
    else
       plot_flag=0; 
    end


    [~, ~, ~, saccade_times, sac_start_times, sac_stop_times] = detect_saccades_bevfel(...
        cur_trial_ET_trace_full', cur_trial_ETsamples, 3, eye_smooth, sgolay_deg, sac_thresh, peri_thresh, plot_flag);
%%
    if plot_intermediate==1
        pause
    end
%% do same processing of ET data as done for saccade detection
    if ET_Eyelink ==3
%        ET_trace_raw_1khz([1:4000]+4000*(tt-1),1:2)=cur_trial_ET_trace_full;
%        ET_trace_raw_1khz([1:4000]+4000*(tt-1),3) = ExperimentRecording{tt, 7}.ET_trace(3,:)';

        cur_trial_ETsamples_kofiko = ExperimentRecording{tt, 1}.m_afEyePositiontimes - ExperimentRecording{tt, 1}.m_afEyePositiontimes(1);
        cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, ExperimentRecording{tt, 1}.m_afEyeXPositionScreenCoordinates' - ExperimentRecording{tt, 1}.m_pt2iFixationSpot(1), linspace(0,4,240))';
        cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, ExperimentRecording{tt, 1}.m_afEyeYPositionScreenCoordinates' - ExperimentRecording{tt, 1}.m_pt2iFixationSpot(2), linspace(0,4,240))';    

    else
        ET_ad_up(1,:) = smooth(cur_trial_ET_trace_full(:,1)', eye_smooth(1), 'sgolay', sgolay_deg(1));
        ET_ad_up(2,:) = smooth(cur_trial_ET_trace_full(:,2)', eye_smooth(1), 'sgolay', sgolay_deg(1));
    
        et_params.eye_fs = 60; %ET_adfreq;
%        ET_times_up = cur_trial_ETsamples(1):(1/et_params.eye_fs):cur_trial_ETsamples(end);
        ET_times_up = 0:(1/et_params.eye_fs):0.2;
        cur_trial_ET_trace(1,:) = interp1(cur_trial_ETsamples, ET_ad_up(1,:), ET_times_up, 'spline');
        cur_trial_ET_trace(2,:) = interp1(cur_trial_ETsamples, ET_ad_up(2,:), ET_times_up, 'spline'); 
     
        eye_smooth2=3;
        cur_trial_ET_trace(1,:) = smooth(cur_trial_ET_trace(1,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
        cur_trial_ET_trace(2,:) = smooth(cur_trial_ET_trace(2,:), eye_smooth(2), 'sgolay', sgolay_deg(2));

    end
%%
%     cur_sac_start_inds=[]; cur_sac_stop_inds=[];
%     for sacs=1:length(saccade_times);
%         cur_sac_start_inds(sacs)=find(edges>sac_start_times(sacs),1,'first');
%         if isempty(find(edges>sac_stop_times(sacs),1,'first'))
%             cur_sac_stop_inds(sacs)=length(edges);
%         else
%         cur_sac_stop_inds(sacs)=find(edges>sac_stop_times(sacs),1,'first');
%         end
%     end

%     ET_trace_raw(cur_trlinds,1) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(ET_trialsec,1), ET_times_up, 'spline');
%     ET_trace_raw(cur_trlinds,2) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(ET_trialsec,2), ET_times_up, 'spline'); 
%     
%     ET_trace(cur_trlinds,1)=cur_trial_ET_trace(1,:)';
%     ET_trace(cur_trlinds,2)=cur_trial_ET_trace(2,:)';
    ET_trace(cur_trlinds,1)=cur_trial_ET_trace(1,:)';
    ET_trace(cur_trlinds,2)=cur_trial_ET_trace(2,:)';

%     sacc_inds=[sacc_inds; cur_trlinds([cur_sac_start_inds' cur_sac_stop_inds'])];
%     bad_inds_sac=[bad_inds_sac, cur_trlinds(cur_sac_start_inds(sacs)):cur_trlinds(cur_sac_stop_inds(sacs))];
    
    if tt>=2
        Block_offsetinds=[Block_offsetinds,cur_trlinds(1)-1]; %end of last block
        Block_onsetinds=[Block_onsetinds,cur_trlinds(1)]; %start of this block block
    end

stim(cur_trlinds,:)=repmat(127*ExperimentRecording{1, 1}.DiscprobeColorID(4:6),trlbins,1);
stim_polar(cur_trlinds,:)=repmat(ExperimentRecording{1, 1}.DiscprobeColorID(1:3),trlbins,1);

%}     
    if mod(tt,100)==0
        disp(['finished trial ' num2str(tt) ' of ' num2str(ntrls)])
    end

end

%%
tvec=[1:length(binned_SU1)];
Block_offsetinds=[Block_offsetinds,length(binned_SU1)];
Block_inds=[Block_onsetinds;Block_offsetinds];

bad_inds_fix=sort([Block_inds(1,:), Block_inds(1,:)+1, Block_inds(1,:)+2,Block_inds(1,:)+3,Block_inds(1,:)+4,Block_inds(1,:)+5,Block_inds(1,:)+6],1);
use_inds_fix=setdiff(tvec,unique(bad_inds_fix));

%%
disp('converting...') 
cd(strExperimentPath)

%%
nofix=1;

exptname=filenameP;
exptdate=filenameP(1:6);

cur_filename=['Jocamo_' filenameP(1:6) '_' arraylabel '_' 'DiscProbe_v08.mat'];
cur_filename_LFP=['Jocamo_' filenameP(1:6) '_' arraylabel '_' 'DiscProbe_v08_LFP.mat'];

stim=stim';
stim_polar=stim_polar';

%%
valid_data=use_inds_fix;
block_inds = Block_inds;

trial_start_ts = trialstart_plx';
trial_start_inds = round(trialstart_plx'.*1000);

sacc_inds=sacc_inds';
ETtrace_raw=ET_trace_raw_1khz';
ETtrace=ET_trace';
useLeye=useLeye';
useReye=useReye';

targchans_SU=targchans(targchans<=nSU);

SU_clusters=[];
Robs=binned_SU1(:,targchans_SU)';
Robs_probe_ID=spk_channels_SU(targchans_SU');
Robs_rating=spk_rating_SU(targchans_SU');

datafilts=ones(size(Robs));

targchansMU=targchans(targchans>nSU);
RobsMU=binned_SU1(:,targchansMU)';
RobsMU_probe_ID=spk_channels_MU(targchansMU'-nSU);
RobsMU_rating=spk_rating_MU(targchansMU'-nSU);

datafiltsMU=ones(size(RobsMU));
%%
spk_times=[]; spk_IDs=[];
for cc=1:length(spk_times_all)
    spk_times = [spk_times,spk_times_all{cc}'];
    spk_IDs = [spk_IDs, ones(1,length(spk_times_all{cc}))*cc];
end
%%
disp('saving...') 
%%
save(cur_filename, 'exptname', 'exptdate', 'stim', 'stim_polar', 'valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw', ...
    'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )

save(cur_filename_LFP, 'trial_start_ts', 'trial_start_inds', '-v7.3' )

disp('Done with Disc Probe model data') 
