%% process and package hartleys
% To run this script, you need the ExperimentRecording file from Step 1
% and the kilosort metadata loaded ins step 1

% this script outputs a processed HDF5 file containing the hartley data,
% including stim (packaged in the same way as clouds) and metadata

stimFilePath = '/media/felix/Internal_1/Data/BevilColor/Cloudstims_calib_01_2022/';
%% params
if useofflinesorting==1
    nChans=nSU+nMU;
else
nChans=256; %targchans=1:24;
end

targ_ETstimtype=7; % 7 for CC, 1 for 1D
targ_stimtype=6;

skipLP=0;
skipET=0;

plot_intermediate=0;

if any([strfind(filenameP,'220203'), strfind(filenameP,'220205'), strfind(filenameP,'220207')])
trlbins=160; dt=.024; repframes=3;
else
trlbins=240; dt=.016; repframes=2;
end
fixvarcutoff=6;
blink_thresh = 60;
cur_BlockID=1;         
% %%
% for cc=1:nChans;
% numSpks(cc)=sum(cell2mat(ExperimentRecording(:,11+3*(cc-1))));
% end
% targchans=find(numSpks>2000);
targchans
% %%
% if useofflinesorting==1
% targchans=1:nChans;
% end
%% for testing purposes
%/{
for tt=1:length(ExperimentRecording)
 %   test(tt)=ExperimentRecording{tt, 1}.DualstimPrimaryuseRGBCloud; 
    test(tt)=ExperimentRecording{tt, 1}.DualstimSecondaryUseCloud; 
    %test(tt)=ExperimentRecording{tt, 1}.m_aiStimulusArea; 
    %test(tt)=ExperimentRecording{tt, 1}.m_aiStimulusRect(1);
%    test(tt)=ExperimentRecording{tt, 1}.m_aiStimulusRect(2);
end
[unique(test), 0, mode(test) length(find(test==mode(test)))]
%}
%% get target indices

for switch_stimtype=0:8;

    targ_trials=[];
    for tt=1:length(ExperimentRecording)
    %    if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Dense Noise');
        if strcmp(ExperimentRecording{tt, 1}.m_strTrialType, 'Dual Stim') && ...
                ExperimentRecording{tt, 1}.m_bMonkeyFixated==1 && ...
                ExperimentRecording{tt, 1}.DualstimPrimaryuseRGBCloud==switch_stimtype &&...
                ExperimentRecording{tt, 1}.DualstimSecondaryUseCloud==targ_ETstimtype;%  &&...
                %ExperimentRecording{tt, 1}.m_aiStimulusRect(1)==975;
        targ_trials=[targ_trials,tt];
        end
    end
    
switch switch_stimtype
    case 0; curstimstype='GT';  ExperimentRecording_GT=ExperimentRecording(targ_trials,:);
%save([strExperimentPath,filesep,'ExperimentRecording4modeling' curstimstype '.mat'],'ExperimentRecording_GT','-v7.3')
    case 3; curstimstype='HL';  ExperimentRecording_HartleyLum=ExperimentRecording(targ_trials,:);
%save([strExperimentPath,filesep,'ExperimentRecording4modeling' curstimstype '.mat'],'ExperimentRecording_HartleyLum','-v7.3')
    case 6; curstimstype='HC';  ExperimentRecording_HartleyCol=ExperimentRecording(targ_trials,:);
%save([strExperimentPath,filesep,'ExperimentRecording4modeling' curstimstype '.mat'],'ExperimentRecording_HartleyCol','-v7.3')
    case 8; curstimstype='CC';    ExperimentRecording_ColCloud=ExperimentRecording(targ_trials,:);
%save([strExperimentPath,filesep,'ExperimentRecording4modeling' curstimstype '.mat'],'ExperimentRecording_ColCloud','-v7.3')
end

end

%%
%{
for tt=1:length(ExperimentRecording)
figure(1);
plot(ExperimentRecording{tt, 1}.m_afEyeXPositionScreenCoordinates-(ExperimentRecording{tt, 1}.m_pt2iFixationSpot(1))); hold on
plot(ExperimentRecording{tt, 1}.m_afEyeYPositionScreenCoordinates-(ExperimentRecording{tt, 1}.m_pt2iFixationSpot(2))); hold off

title(num2str(ExperimentRecording{tt, 1}.m_bMonkeyFixated));
%ExperimentRecording{4353, 1}.m_afEyeXPositionScreenCoordinates  

figure(2);
plot(ExperimentRecording{tt, 5}');

pause
end
%}


%% Load Hartleys
load('/media/felix/Internal_1/Data/BevilColor/hartleys_60.mat')
hartleys=hartleys60_DKL;
hartleys_metas=hartleys60_meta;
hartleys_metas(769:end,4)=3;

load([ExperimentFolder filenameP '.mat'], 'g_astrctAllParadigms')
DualstimETbars = int8(squeeze((g_astrctAllParadigms{1, 1}.DualstimETbars-128)/127)');
%%
switch targ_stimtype
    case 3; curstimstype='HL';  ExperimentRecording_mod=ExperimentRecording_HartleyLum;
    case 6; curstimstype='HC';  ExperimentRecording_mod=ExperimentRecording_HartleyCol;
    case 7; curstimstype='LC';    ExperimentRecording_mod=ExperimentRecording_LumCloud;
    case 8; curstimstype='CC';    ExperimentRecording_mod=ExperimentRecording_ColCloud;
end
switch targ_ETstimtype
    case 1; curETstimtype='1D';
    case 7; curETstimtype='CC';
end
%%
ntrls = size(ExperimentRecording_mod,1);
NT=trlbins*ntrls;

trlsecs = 4;

if ET_Eyelink==3
    numETtraces=2;
    ET_trace_raw_1khz=zeros(ntrls*trlsecs*1000,3);
else
numETtraces=size(ExperimentRecording_mod{end, 7}.ET_trace,1);
ET_trace_raw_1khz=zeros(ntrls*trlsecs*1000,numETtraces);
end
ET_ad_up=[];

if any([strfind(filenameP,'220203'), strfind(filenameP,'220205'), strfind(filenameP,'220207')])
edges=linspace(0,2.67,trlbins); edges_hist=linspace(0,2.67,trlbins+1);
else
edges=linspace(0,trlsecs,trlbins); edges_hist=linspace(0,trlsecs,trlbins+1);
end

stim=int8(zeros(NT,60,60,3)); stimtype=zeros(NT,1);  

if ~skipET
    if targ_ETstimtype==1
    stimET=int8(zeros(NT,60)); stimtypeET=zeros(NT,1); stimETori=zeros(NT,1); 
    else
    stimET=int8(zeros(NT,60,60,3)); stimtypeET=zeros(NT,1); 
    end
end

trialstart_plx=zeros(ntrls,1);

hartleystim_metas=zeros(NT,4); 
binned_SU1=int8(zeros(NT,nChans)); binned_MUA=[]; ET_trace=zeros(NT,2); ET_trace_raw=zeros(NT,2);
%binned_SU2=int8(zeros(NT,nChans)); binned_SU3=int8(zeros(NT,nChans)); %binned_SU4=zeros(NT,nChans); 
bad_inds_fix=[1:10]; bad_inds_sac = []; sacc_inds=[];
Block_onsetinds=1; Block_offsetinds=[]; Block_onsetinds_blink=1; Block_offsetinds_blink=[];
BlockID=[]; TrialID=[]; useLeye=zeros(NT,1); useReye=zeros(NT,1);
spk_times_all=cell([nSU+nMU,1]);

cloud_scale=zeros(NT,1); 
cloud_binary=zeros(NT,1);     

cd(stimFilePath)
cur_scale=g_astrctAllParadigms{1, 1}.DualstimScale.Buffer(find(g_astrctAllParadigms{1, 1}.DualstimScale.TimeStamp<ExperimentRecording_mod{1,2},1,'last'));
load(sprintf(['Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID));
%load(sprintf([stimFilePath 'Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID)) 
DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
%%

for tt=1:ntrls
    cur_trlinds=[1:trlbins]+trlbins*(tt-1);
    pixelscaf=round(ExperimentRecording_mod{tt, 1}.m_aiStimulusArea/60);

    trialstart_plx(tt) = ExperimentRecording_mod{tt, 1}.PlexonOnsetTime;

    if ~isfield(ExperimentRecording_mod{tt,1}, 'usebinary')
        ExperimentRecording_mod{tt, 1}.usebinary=0;
    end
    
    if targ_stimtype==8
        if ExperimentRecording_mod{tt,1}.BlockID~=cur_BlockID;
            cur_BlockID = ExperimentRecording_mod{tt,1}.BlockID;
            cur_scale=g_astrctAllParadigms{1, 1}.DualstimScale.Buffer(find(g_astrctAllParadigms{1, 1}.DualstimScale.TimeStamp<ExperimentRecording_mod{tt,2},1,'last'));
            if ExperimentRecording_mod{tt, 1}.usebinary  
            load(sprintf(['Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat'], cur_scale, cur_BlockID));            
            else
            load(sprintf(['Cloudstims_Chrom_size60_scale%d_%02d.mat'], cur_scale, cur_BlockID));
            end
            DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        end
        cur_TrialID = ExperimentRecording_mod{tt,1}.TrialID;
        BlockID(cur_trlinds)=cur_BlockID;
        TrialID(cur_trlinds)=cur_TrialID;
    
        cloud_scale(cur_trlinds) = cur_scale;
        cloud_binary(cur_trlinds) = ExperimentRecording_mod{tt, 1}.usebinary; 
    end

    if ~isfield(ExperimentRecording_mod{tt,1}, 'UseLeye')
        ExperimentRecording_mod{tt, 1}.UseLeye=1;
        ExperimentRecording_mod{tt, 1}.UseReye=1;
    end
    useLeye(cur_trlinds)=ExperimentRecording_mod{tt,1}.UseLeye;
    useReye(cur_trlinds)=ExperimentRecording_mod{tt,1}.UseReye;
    
%/{ 
% comment out to only get ET info    
    for channel=1:nChans
        if ~isempty(ExperimentRecording_mod{tt,10+3*(channel-1)})
            try
                cur_spks=ExperimentRecording_mod{tt,10+3*(channel-1)}.unit1;
                cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
                binned_SU1(cur_trlinds,channel)=histcounts(ExperimentRecording_mod{tt,10+3*(channel-1)}.unit1,edges_hist);
                spk_times_all{channel}=[spk_times_all{channel}; cur_spks+(trlsecs*(tt-1))];
            catch
                cur_spks=ExperimentRecording_mod{tt,10+3*(channel-1)};
                cur_spks(cur_spks<0)=[]; cur_spks(cur_spks>trlsecs)=[];
                binned_SU1(cur_trlinds,channel)=histcounts(ExperimentRecording_mod{tt,10+3*(channel-1)},edges_hist);
                spk_times_all{channel,1}=[spk_times_all{channel,1}; cur_spks+(trlsecs*(tt-1))];
            end
        end   
    end
%}

%%
    cur_trial_ET_trace_full=[]; cur_trial_ET_trace=[];
    cur_trial_ETsamples=[0:.001:3.999];
    cur_trial_ET_trace_full(:,1)=(ExperimentRecording_mod{tt, 7}.ET_trace(1,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
    cur_trial_ET_trace_full(:,2)=(ExperimentRecording_mod{tt, 7}.ET_trace(2,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);

    if numETtraces>2
        cur_trial_ET_trace_full(:,3)=(ExperimentRecording_mod{tt, 7}.ET_trace(3,:)' )*(g_strctEyeCalib.GainX.Buffer(end)./1000);
        cur_trial_ET_trace_full(:,4)=(ExperimentRecording_mod{tt, 7}.ET_trace(4,:)' )*(g_strctEyeCalib.GainY.Buffer(end)./1000);
    end
% % if needing to median-adjust:
%     cur_trial_ET_trace_full(:,1)=(ExperimentRecording_mod{tt, 7}.ET_trace(1,:)' - median(ExperimentRecording_mod{tt, 7}.ET_trace(1,:)'))*(g_strctEyeCalib.GainX.Buffer(end)./1000);
%     cur_trial_ET_trace_full(:,2)=(ExperimentRecording_mod{tt, 7}.ET_trace(2,:)' - median(ExperimentRecording_mod{tt, 7}.ET_trace(2,:)'))*(g_strctEyeCalib.GainY.Buffer(end)./1000);    

% % if needing to interpolate time positions
%         cur_trial_ET_trace(:,1)=interp1(cur_trial_ETsamples,cur_trial_ET_trace_full(:,1),edges,'linear');
%         cur_trial_ET_trace(:,2)=interp1(cur_trial_ETsamples,cur_trial_ET_trace_full(:,2),edges,'linear');

%%
if ET_Eyelink==1
    eye_smooth = [15, 3, 3]; sgolay_deg = [2,2,2];
    sac_thresh = 9; %3.5; %threshold eye speed % default 6 for 0314
    peri_thresh = 3; %threshold eye speed for defining saccade boundary inds % default 2.5 for 0314
else
     eye_smooth = [181, 3, 5]; sgolay_deg = [3,2,4];
     sac_thresh = 2.5; %3.5; %threshold eye speed
     peri_thresh = 1.2; %threshold eye speed for defining saccade boundary inds
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
    ET_trace_raw_1khz([1:4000]+4000*(tt-1),1:2)=cur_trial_ET_trace_full;
    ET_trace_raw_1khz([1:4000]+4000*(tt-1),3) = ExperimentRecording_mod{tt, 7}.ET_trace(3,:)';

    cur_trial_ETsamples_kofiko = ExperimentRecording_mod{tt, 1}.m_afEyePositiontimes - ExperimentRecording_mod{tt, 1}.m_afEyePositiontimes(1);
    cur_trial_ET_trace(1,:)=interp1(cur_trial_ETsamples_kofiko, ExperimentRecording_mod{tt, 1}.m_afEyeXPositionScreenCoordinates' - ExperimentRecording_mod{tt, 1}.m_pt2iFixationSpot(1), linspace(0,4,240))';
    cur_trial_ET_trace(2,:)=interp1(cur_trial_ETsamples_kofiko, ExperimentRecording_mod{tt, 1}.m_afEyeYPositionScreenCoordinates' - ExperimentRecording_mod{tt, 1}.m_pt2iFixationSpot(2), linspace(0,4,240))';    

else
    ET_trace_raw_1khz([1:4000]+4000*(tt-1),:)=cur_trial_ET_trace_full;
    ET_ad_up(1,:) = smooth(cur_trial_ET_trace_full(:,1)', eye_smooth(1), 'sgolay', sgolay_deg(1));
    ET_ad_up(2,:) = smooth(cur_trial_ET_trace_full(:,2)', eye_smooth(1), 'sgolay', sgolay_deg(1));

    et_params.eye_fs = 60; %ET_adfreq;
    ET_times_up = cur_trial_ETsamples(1):(1/et_params.eye_fs):cur_trial_ETsamples(end);
    cur_trial_ET_trace(1,:) = interp1(cur_trial_ETsamples, ET_ad_up(1,:), ET_times_up, 'spline');
    cur_trial_ET_trace(2,:) = interp1(cur_trial_ETsamples, ET_ad_up(2,:), ET_times_up, 'spline'); 
 
    eye_smooth2=3;
    cur_trial_ET_trace(1,:) = smooth(cur_trial_ET_trace(1,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
    cur_trial_ET_trace(2,:) = smooth(cur_trial_ET_trace(2,:), eye_smooth(2), 'sgolay', sgolay_deg(2));

end
%%
    cur_sac_start_inds=[]; cur_sac_stop_inds=[];
    for sacs=1:length(saccade_times);
        cur_sac_start_inds(sacs)=find(edges>sac_start_times(sacs),1,'first');
        if isempty(find(edges>sac_stop_times(sacs),1,'first'))
            cur_sac_stop_inds(sacs)=length(edges);
        else
        cur_sac_stop_inds(sacs)=find(edges>sac_stop_times(sacs),1,'first');
        end
    end

%     ET_trace_raw(cur_trlinds,1)=cur_trial_ET_trace_full(:,1);
%     ET_trace_raw(cur_trlinds,2)=cur_trial_ET_trace_full(:,2);
%    ET_trace(cur_trlinds,1)=smooth(ET_trace_raw(cur_trlinds,1),5);
%    ET_trace(cur_trlinds,2)=smooth(ET_trace_raw(cur_trlinds,2),5);

%    ET_trace_raw(cur_trlinds,1) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(:,1), ET_times_up, 'spline');
%    ET_trace_raw(cur_trlinds,2) = interp1(cur_trial_ETsamples, cur_trial_ET_trace_full(:,2), ET_times_up, 'spline'); 
    
    ET_trace(cur_trlinds,1)=cur_trial_ET_trace(1,:)';
    ET_trace(cur_trlinds,2)=cur_trial_ET_trace(2,:)';

    sacc_inds=[sacc_inds; cur_trlinds([cur_sac_start_inds' cur_sac_stop_inds'])];
    bad_inds_sac=[bad_inds_sac, cur_trlinds(cur_sac_start_inds(sacs)):cur_trlinds(cur_sac_stop_inds(sacs))];

    %{
tt=3;
figure;
subplot(2,1,1); plot(ExperimentRecording_ColCloud{tt, 7}.ET_trace'); title('Plexon ET traces')
subplot(2,1,2); plot(ExperimentRecording_ColCloud{tt, 1}.ETthisTrialRawEyeData); title('Kofiko ET traces'); axis tight

    %}
% previous ET integration code    
%{
    try
    cur_trial_ETsamples=ExperimentRecording_mod{tt, 1}.m_afEyePositiontimes-ExperimentRecording_mod{tt, 1}.m_afEyePositiontimes(1)+g_strctStatistics.m_strctEyeData.m_fEyeIntegrationPeriod(1);
    catch
%        cur_trial_ETsamples=2*[1:length(ExperimentRecording_mod{tt,1}.m_afEyeXPositionScreenCoordinates)]./length(ExperimentRecording_mod{tt,1}.m_afEyeXPositionScreenCoordinates);
        sprintf('ET data missing? \n')
    end
    cur_trial_ET_trace(:,1)=([interp1(cur_trial_ETsamples,ExperimentRecording_mod{tt,1}.m_afEyeXPositionScreenCoordinates,edges,'linear')-ExperimentRecording_mod{tt,1}.m_pt2iFixationSpot(1)])./pixelscaf;
    cur_trial_ET_trace(:,2)=([interp1(cur_trial_ETsamples,ExperimentRecording_mod{tt,1}.m_afEyeYPositionScreenCoordinates,edges,'linear')-ExperimentRecording_mod{tt,1}.m_pt2iFixationSpot(2)])./pixelscaf;


    ET_trace_raw(cur_trlinds,1)=cur_trial_ET_trace(:,1);
    ET_trace_raw(cur_trlinds,2)=cur_trial_ET_trace(:,2);

    ET_trace(cur_trlinds,1)=ET_trace_raw(cur_trlinds,1);
    ET_trace(cur_trlinds,2)=ET_trace_raw(cur_trlinds,2);
    
    if ~isempty(ExperimentRecording_mod{tt,6}.saccades);
%         cur_sac_inds=[1];
%         sacc_inds=[sacc_inds; [cur_trlinds(1)-1, cur_trlinds(1)]];
%         sac_edges=[0,1];
        cur_sac_inds=[]; cur_sac_start_inds=[]; cur_sac_stop_inds=[];
        
        if length(ExperimentRecording_mod{tt,6}.saccade_stop)<length(ExperimentRecording_mod{tt,6}.saccade_start)
            ExperimentRecording_mod{tt,6}.saccade_start=ExperimentRecording_mod{tt,6}.saccade_start(1:length(ExperimentRecording_mod{tt,6}.saccade_stop));
            ExperimentRecording_mod{tt,6}.saccades=ExperimentRecording_mod{tt,6}.saccades(1:length(ExperimentRecording_mod{tt,6}.saccade_stop));
        end
        if any(ExperimentRecording_mod{tt,6}.saccades-ExperimentRecording_mod{tt,2}>max(edges))
            ExperimentRecording_mod{tt,6}.saccades(ExperimentRecording_mod{tt,6}.saccades>max(edges))=[];
        end
        
        for sacs=1:length(ExperimentRecording_mod{tt,6}.saccades);
            curtrl_sactimes = ExperimentRecording_mod{tt,6}.saccades-ExperimentRecording_mod{tt,2};
            curtrl_sactimes_start = ExperimentRecording_mod{tt,6}.saccade_start-ExperimentRecording_mod{tt,2};
            curtrl_sactimes_stop = ExperimentRecording_mod{tt,6}.saccade_stop-ExperimentRecording_mod{tt,2};
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
%     %             ET_trace(cur_trlinds(cur_fixs_adj_inds),1)=mean(ExperimentRecording_mod{tt,3}(1,cur_fixs_adj_inds));
%     %             ET_trace(cur_trlinds(cur_fixs_adj_inds),2)=mean(ExperimentRecording_mod{tt,3}(2,cur_fixs_adj_inds));
%         ET_trace(cur_trlinds(cur_fixs_adj_inds),1)=mean(ET_trace(cur_trlinds(cur_fixs_adj_inds),1));
%         ET_trace(cur_trlinds(cur_fixs_adj_inds),2)=mean(ET_trace(cur_trlinds(cur_fixs_adj_inds),2));
%     end
    
    if tt>=2
        Block_offsetinds=[Block_offsetinds,cur_trlinds(1)-1]; %end of last block
        Block_onsetinds=[Block_onsetinds,cur_trlinds(1)]; %start of this block block
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
    switch ExperimentRecording_mod{tt,1}.DualstimPrimaryuseRGBCloud
        case 8
            stim(cur_trlinds,:,:,:)=permute( DensenoiseChromcloud_DKlspace(:,:,ExperimentRecording_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:),[3 1 2 4]);
            %stimulus(cur_trlinds_stim,:,:,:)=permute(ExperimentRecording_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
            stimtype(cur_trlinds)=8;
    
        case 7
            stim(cur_trlinds,:,:,:)=permute( DensenoiseAchromcloud_binned(:,:,ExperimentRecording_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:),[3 1 2 4]);
            %stimulus(cur_trlinds_stim,:,:,:)=permute(ExperimentRecording_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
            stimtype(cur_trlinds)=8;
            
        case 6
            stim(cur_trlinds,:,:,:)=hartleys60_DKL(ExperimentRecording_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:,:,:);
            stimtype(cur_trlinds)=6;
            hartleystim_metas(cur_trlinds,:)=hartleys_metas(ExperimentRecording_mod{tt,1}.stimseq(1:repframes:repframes*trlbins),:);

        case 3
            stim(cur_trlinds_stim,:,:,:)=hartleys60_DKL(ExperimentRecording_mod{tt,1}.stimseq(1:repframes:end),:,:,:);
            stimtype(cur_trlinds)=3;
            hartleystim_metas(cur_trlinds_stim,:)=hartleys_metas(ExperimentRecording_mod{tt,1}.stimseq(1:2:end),:);
    end
    
    if ~skipET
        switch ExperimentRecording_mod{tt, 1}.DualstimSecondaryUseCloud
            case 1
                stimET(cur_trlinds,:)=DualstimETbars(ExperimentRecording_mod{tt, 1}.stimseq_ET_bars(1:repframes:end),:);
                %stimulus(cur_trlinds_stim,:,:,:)=permute(ExperimentRecording_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
                stimETori(cur_trlinds,:)=ExperimentRecording_mod{tt, 1}.stimseq_ET_baroris(1:repframes:end);
                stimtypeET(cur_trlinds)=1;
            case 7
                stimET(cur_trlinds,:,:,:)=permute( DensenoiseChromcloud_DKlspace(:,:,ExperimentRecording_mod{tt,1}.stimseq_ET_Cclouds(1:repframes:repframes*trlbins),:),[3 1 2 4]);
                %stimulus(cur_trlinds_stim,:,:,:)=permute(ExperimentRecording_mod{tt,1}.stimuli(:,:,1:2:end,:),[2 3 1 4]);
                stimtypeET(cur_trlinds)=7;            
        end
    end
%}     
    if mod(tt,50)==0
        disp(['finished trial ' num2str(tt) ' of ' num2str(ntrls)])
    end

end


%%
tvec=[1:length(binned_SU1)];
Block_offsetinds=[Block_offsetinds,length(binned_SU1)];
%Block_offsetinds_blink=[Block_offsetinds_blink,length(binned_SU1)];
%Block_inds=[unique([Block_onsetinds,Block_onsetinds_blink]);unique([Block_offsetinds,Block_offsetinds_blink])];
Block_inds=[Block_onsetinds;Block_offsetinds];

bad_inds_fix=sort([Block_inds(1,:), Block_inds(1,:)+1, Block_inds(1,:)+2,Block_inds(1,:)+3,Block_inds(1,:)+4,Block_inds(1,:)+5,Block_inds(1,:)+6],1);
use_inds_fix=setdiff(tvec,unique(bad_inds_fix));

%% Hartley analysis
%{
% get target indices
targchans=find(sum(binned_SU1)>200);

targchans_SU=targchans(targchans<=nSU);
targchans_MU=targchans(targchans>nSU);
Robs_probe_ID=spk_channels_SU(targchans_SU');
RobsMU_probe_ID=spk_channels_MU(targchans_MU'-nSU);
targchans_probe_ID = [Robs_probe_ID; RobsMU_probe_ID];
spike_ID=double([spk_ID_SU; spk_ID_MU]);

% targchans_probeindices = find(targchans_probe_ID<24);
% targchans=targchans(targchans_probeindices);
% targchans_probe_ID = targchans_probe_ID(targchans_probeindices);
% spike_ID = spike_ID(targchans_probeindices);

%% STA along hartley axes
hart_use_inds=find(stimtype<=6);
hart_use_inds(end-10:end)=[];

freqvec=unique(hartleystim_metas(:,1));
rotvec=unique(hartleystim_metas(:,2));
shiftvec=unique(hartleystim_metas(:,3));
colvec=unique(hartleystim_metas(:,4));

binned_SU=double(binned_SU1);
%%
base_use_inds1=intersect(hart_use_inds,use_inds_fix);

targlag = 3;

ii=1; all_cells_tuning=[];
tuningfig = figure('position', [300 500 1200 500]);

for cc=targchans;
    cur_tuning_all=[];
    base_use_inds=base_use_inds1;
    for ff=1:length(freqvec)
        for oo=1:length(rotvec)
            for ss=1:length(shiftvec)
                for ccol=1:length(colvec)
                    cur_freqinds=find(hartleystim_metas(:,1)==freqvec(ff));
                    cur_rotinds=find(hartleystim_metas(:,2)==rotvec(oo));
                    cur_shiftinds=find(hartleystim_metas(:,3)==shiftvec(ss));
                    cur_colinds=find(hartleystim_metas(:,4)==colvec(ccol));
                    
                    cur_use_inds1=intersect(cur_freqinds,cur_rotinds);
                    cur_use_inds2=intersect(cur_use_inds1,cur_shiftinds);
                    cur_use_inds3=intersect(cur_use_inds2,cur_colinds);
                    cur_use_inds=intersect(base_use_inds,cur_use_inds3);

                    for curlag=1:6
    %                    cur_STA_all(ff,oo,ss,curlag)=sum(binned_MUA(cur_use_inds+curlag,cc))./length(cur_use_inds);
                        cur_tuning_all(ff,oo,ss,ccol,curlag)=sum(binned_SU(cur_use_inds+curlag,cc))./length(cur_use_inds);
                    end
                end
            end
        end
    end

    [~,targlag]=max(var(var(var(var(cur_tuning_all,1,4),1,3),1,2),1,1));

    cur_tuning_lum = mean(squeeze(cur_tuning_all(:,:,:,1,targlag)),3);
    cur_tuning_LM = mean(squeeze(cur_tuning_all(:,:,:,2,targlag)),3);
    cur_tuning_S = mean(squeeze(cur_tuning_all(:,:,:,3,targlag)),3);
    maxlim=max([cur_tuning_lum(:); cur_tuning_LM(:); cur_tuning_LM(:)]);
    
    subplot(1,3,1); imagesc(cur_tuning_lum); title('Lum');
        ylabel('SF'); yticklabels(freqvec); xlabel('Orientation'); xticklabels(rotvec(2:2:end)); caxis([0 maxlim]);
    subplot(1,3,2); imagesc(cur_tuning_LM); title('L-M');  xticklabels(rotvec(2:2:end)); yticklabels(freqvec); caxis([0 maxlim]); 
    subplot(1,3,3); imagesc(cur_tuning_S); title('S'); xticklabels(rotvec(2:2:end)); yticklabels(freqvec); caxis([0 maxlim]); colorbar
    
    figtitle(['Spike ID: ' num2str(spike_ID(ii)) ' - Probe ' num2str(targchans_probe_ID(ii)) ' - Spikes: ' num2str(sum(binned_SU(base_use_inds1,cc))) ' - Firing rate: ' num2str(mean(binned_SU(base_use_inds1,cc))/dt)])
    saveas(tuningfig,[strExperimentPath 'Spike ID ' num2str(spike_ID(ii)) ' - Probe ' num2str(targchans_probe_ID(ii)) '_HartleyTuning.png'])
    %{
    figure(1);
    subplot(2,1,1)
%    plot(rotvec,mean(mean(mean(cur_STA_all,4),3),1)*10); title('Orientation')
    cur_STAs_rot=squeeze(mean(mean(cur_STA_all,3),1));
%    cur_STAvars_rot=squeeze(mean(var(cur_STA_all,1,3),1));
    plot(rotvec,cur_STAs_rot*60); title(['SU Unit ' num2str(cc)])
    xlabel('Orientation (deg)'); ylabel('Firing rate (Hz)')
     
    subplot(2,1,2)
    [~,bestlag_rot]=max(var(cur_STAs_rot));
    [~,bestori]=max(cur_STAs_rot(:,bestlag_rot));
%    [~,bestori]=max(var(cur_STAs_rot'));
%    plot(freqvec/(pi/3),mean(mean(mean(cur_STA_all,4),3),2)*10); title('Spatial Frequency');set(gca, 'XScale', 'log')
    cur_STAs_freq=squeeze(max(cur_STA_all(:,bestori,:,:),[],3));
    %    cur_STAs_freqFull=squeeze(mean(mean(cur_STA_all,3),2));
    cur_STAvars_freq=squeeze(mean(var(cur_STA_all,1,2),3));

    plot(freqvec/(pi/3),cur_STAs_freq*60); title(['Spatial Frequency at ori' num2str(rotvec(bestori))]);%set(gca, 'XScale', 'log')
    xlabel('Freq (cpd)'); ylabel('Firing rate (Hz)')
    
%    figtitle(['SU Unit ' num2str(cc)])
    
%     figure(2)
%     imagesc(squeeze(mean(cur_STA_all(:,:,:,3),3))')
%     xlabel('Frequency steps')
%     ylabel('Orientations')
%     title(['SU Unit ' num2str(cc)])
%     
    figure(3)
    subplot(2,1,1)

    plot(rotvec,cur_STAs_rot(:,bestlag_rot)*40); title(['SU Unit ' num2str(cc)])
    xlabel('Orientation (deg)'); ylabel('Firing rate (Hz)')
%    ylim([0 max(cur_STAs_rot(:,bestlag_rot)*40)])
    
    subplot(2,1,2)
    [~,bestlag_freq]=max(var(cur_STAs_freq));
%    plot(freqvec/(pi/3),cur_STAs_freq*15); title(['Spatial Frequency at ori' num2str(rotvec(bestori))]);%set(gca, 'XScale', 'log')
    plot(freqvec/(pi/3),cur_STAs_freq(:,bestlag_freq)*40); title(['Spatial Frequency at ori' num2str(rotvec(bestori))]);%set(gca, 'XScale', 'log')
    xlabel('Freq (cpd)'); ylabel('Firing rate (Hz)')
%    ylim([0 max(cur_STAs_freq(:,bestlag_freq)*40)])
    
%    figtitle(['Best lags for SU Unit ' num2str(cc)])
    
%     figure(4)
%     plot(freqvec/(pi/3),cur_STAvars_freq*15); title('Spatial Frequency');%set(gca, 'XScale', 'log')
%}
    all_cells_tuning(cc).probe = targchans_probe_ID(ii);
    all_cells_tuning(cc).unitID = spike_ID(ii);
    all_cells_tuning(cc).all_tuning=cur_tuning_all;
    all_cells_tuning(cc).Lum=cur_tuning_lum;
    all_cells_tuning(cc).LM=cur_tuning_LM;
    all_cells_tuning(cc).S=cur_tuning_S;
    
    ii=ii+1;
end  

cur_filename=[strExperimentPath 'Jocamo_' filenameP(1:6) '_' arraylabel '_hartleytuning.mat'];
save(cur_filename, 'all_cells_tuning');

%% STA across several lags
%cur_use_inds=intersect(find(stimtype==6),use_inds_fix);
cur_use_inds=use_inds_fix;
cur_use_inds(end-10:end)=[];
%stim2=reshape(stim,size(stim,1),3*60*60);
stim2=(double(reshape(stim,size(stim,1),3*60*60))./127);

%%
%{
teststim=reshape(squeeze(stimulus(1,:,:,:)),1,3*60*60);
teststim2=reshape(teststim,60,180);
figure; 
subplot(2,3,1); imagesc(squeeze(stimulus(1,:,:,1)))
subplot(2,3,2); imagesc(squeeze(stimulus(1,:,:,2)))
subplot(2,3,3); imagesc(squeeze(stimulus(1,:,:,3)))
teststim2(:,[60 120])=max(teststim2(:));
subplot(2,1,2); imagesc(teststim2)
%}
ii=1;
STAfig = figure('position', [300 500 1200 500]);
for cc=targchans    
    tic
%     for curlag=1:5
% %    cur_STA{curlag}=binned_MUA(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
% %    cur_STA{curlag}=binned_SU1(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
% %    cur_STA{curlag}=binned_SU1(cur_use_inds+curlag,cc)'*stim3(cur_use_inds,:);
% %    cur_STA{curlag}=binned_SU2(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
% %    cur_STA{curlag}=binned_SUSort(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
%     end

    for curlag=1:6
        cur_STA1(curlag,:)=binned_SU(cur_use_inds+curlag,cc)'*stim2(cur_use_inds,:);
    end
    for curlag=1:6
%    cur_StA2=reshape(cur_STA{curlag},3,25,25);
%    cur_StA3=reshape(permute(cur_StA2,[3 2 1]),25,75);
%   cur_StA2=reshape(cur_STA{curlag},60,180);
    cur_StA2=reshape(cur_STA1(curlag,:),60,180);

    %curlim=max(abs(cell2mat(cur_STA(:)')));
    curlim=max(abs(cur_STA1(:)'));
    cur_StA2(:,[60 120])=curlim;
    subplot(1,6,curlag)
    imagesc(cur_StA2'); caxis([-curlim curlim]); pbaspect([1 3 1])
    colormap(gray); xlabel(['Lag ' num2str(curlag)]);
    if curlag==1; ylabel('S          L-M          Lum'); end
    end
%     figtitle(['Unit ' num2str(cc) ' - Probe ' num2str(targchans_probe_ID(ii)) ' - Spikes: ' num2str(sum(binned_SU(base_use_inds1,cc))) ' - Firing rate: ' num2str(mean(binned_SU(base_use_inds1,cc))/dt)])
%     saveas(STAfig,[strExperimentPath 'Unit ' num2str(cc) ' - Probe ' num2str(targchans_probe_ID(ii)) '_HartleySTA.png'])
    figtitle(['Spike ID: ' num2str(spike_ID(ii)) ' - Probe ' num2str(targchans_probe_ID(ii)) ' - Spikes: ' num2str(sum(binned_SU(base_use_inds1,cc))) ' - Firing rate: ' num2str(mean(binned_SU(base_use_inds1,cc))/dt)])
    saveas(STAfig,[strExperimentPath 'Spike ID ' num2str(spike_ID(ii)) ' - Probe ' num2str(targchans_probe_ID(ii)) '_HartleySTA.png'])
    ii=ii+1;

    toc
    %pause
end  
%}
%%
%%
disp('converting...') 
cd(strExperimentPath)

%% now package the data into HDF5
nofix=1;

exptname=filenameP;
exptdate=filenameP(1:6);

if nofix==1
cur_filename=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_nofix_v08.mat'];
%cur_filename_targchans=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_nofix_v08_cloudSUinds.mat'];
cur_filename_targchans=['Jocamo_' filenameP(1:6) '_' arraylabel '_CC_ETCC_nofix_v08_cloudSUinds.mat'];
else
cur_filename=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_v08.mat'];
cur_filename_targchans=['Jocamo_' filenameP(1:6) '_' arraylabel '_CC_ETCC_nofix_v08_cloudSUinds.mat'];
end

%dt=.016;
pixel_size=1;
electrode_info=[];
stim_location=ExperimentRecording_mod{end, 1}.m_aiStimulusRect;
if isfield(ExperimentRecording_mod{end,1}, 'm_aiTiledStimulusRect')
    stim_location=ExperimentRecording_mod{end, 1}.m_aiTiledStimulusRect;
end
ETstim_location=[ExperimentRecording_mod{end, 1}.secondarystim_bar_rect; 
    ExperimentRecording_mod{end, 1}.tertiarystim_bar_rect];
fix_location=ExperimentRecording_mod{end, 1}.m_pt2iFixationSpot;
fix_size=ExperimentRecording_mod{end, 1}.m_fFixationSizePix-1;

stim=permute(stim,[2 3 4 1]);
stimtype=stimtype';
stimscale=(stim_location(3)-stim_location(1))/60;

%%
fixlocs1=[fix_location(1)-fix_size:fix_location(1)+fix_size];
fixlocs2=[fix_location(2)-fix_size:fix_location(2)+fix_size];

if nofix
    fixdot=[];
else
    overlap_d1=any(fixlocs1>stim_location(1) & fixlocs1<stim_location(3));
    dist_d1=fix_location(1)-stim_location(1);
    
    overlap_d2=any(fixlocs2>stim_location(2)&fixlocs2<stim_location(4));
    dist_d2=fix_location(2)-stim_location(2);
    
    if overlap_d1 && overlap_d2
        if stimscale>1
            %todo: get this for new fixation calc
            if mod(fix_size,2)
                if mod(dist_d1,2)
                    fix_d1=ceil((dist_d1-fix_size-1)/stimscale):floor((dist_d1+fix_size-1)/stimscale);
                else
                    fix_d1=ceil((dist_d1-fix_size)/stimscale):floor((dist_d1+fix_size)/stimscale);
                end
    
                if mod(dist_d2,2)
                    fix_d2=ceil((dist_d2-fix_size-1)/stimscale):floor((dist_d2+fix_size-1)/stimscale);
                else
                    fix_d2=ceil((dist_d2-fix_size)/stimscale):floor((dist_d2+fix_size)/stimscale);
                end
            else
                if mod(dist_d1,2)
                    fix_d1=ceil((dist_d1-fix_size)/stimscale):floor((dist_d1+fix_size)/stimscale);
                else
                    fix_d1=ceil((dist_d1-fix_size-1)/stimscale):floor((dist_d1+fix_size-1)/stimscale);
                end
                if mod(dist_d2,2)
                    fix_d2=ceil((dist_d2-fix_size)/stimscale):floor((dist_d2+fix_size)/stimscale);
                else
                    fix_d2=ceil((dist_d2-fix_size-1)/stimscale):floor((dist_d2+fix_size-1)/stimscale);
                end    
            end
            fixdot=[fix_d1;fix_d2];
        else
    %         fix_d1=(dist_d1-fix_size):(dist_d1+fix_size);
    %         fix_d2=(dist_d2-fix_size):(dist_d2+fix_size);    
            fix_d1=(fixlocs1-stim_location(1)); fix_d1=fix_d1(fix_d1>0); fix_d1=fix_d1(fix_d1<61);
            fix_d2=(fixlocs2-stim_location(2)); fix_d2=fix_d2(fix_d2>0); fix_d2=fix_d2(fix_d2<61);
            fixdot=[{fix_d2;fix_d1}];     
        end
        stim(fix_d2,fix_d1,1,:)=0;
        stim(fix_d2,fix_d1,2,:)=0;
        stim(fix_d2,fix_d1,3,:)=0;

    else
        fixdot=[0,0;0,0];
    end
end
%%

%%
if ~skipET
    if targ_ETstimtype==1;
        stimET=stimET';
        stimETori=stimETori';
    else
        stimET=permute(stimET,[2 3 4 1]);
        
        if nofix
            fixdotET=[];
        else
            overlapET1_d1=any(fixlocs1>ETstim_location(1,1) & fixlocs1<ETstim_location(1,3));
            overlapET1_d2=any(fixlocs2>ETstim_location(1,2) & fixlocs2<ETstim_location(1,4));
            
            overlapET2_d1=any(fixlocs1>ETstim_location(2,1) & fixlocs1<ETstim_location(2,3));
            overlapET2_d2=any(fixlocs2>ETstim_location(2,2) & fixlocs2<ETstim_location(2,4));
            
            if any([(overlapET1_d1 && overlapET1_d2), (overlapET2_d1 && overlapET2_d2)])
                fixET_d1=unique([fixlocs1-ETstim_location(1,1),fixlocs1-ETstim_location(2,1)]); fixET_d1=fixET_d1(fixET_d1>0); fixET_d1=fixET_d1(fixET_d1<61);
                fixET_d2=unique([fixlocs2-ETstim_location(1,2),fixlocs2-ETstim_location(2,2)]); fixET_d2=fixET_d2(fixET_d2>0); fixET_d2=fixET_d2(fixET_d2<61);
                fixdotET={fixET_d2;fixET_d1};  
    
                stimET(fixET_d2,fixET_d1,1,:)=0;
                stimET(fixET_d2,fixET_d1,2,:)=0;
                stimET(fixET_d2,fixET_d1,3,:)=0;
            else
                fixdotET=[];
            end
        end
    end
    stimtypeET=stimtypeET';
else
    curETstimtype='none';
end
%%
valid_data=use_inds_fix;
block_inds = Block_inds;

trial_start_ts = trialstart_plx';
trial_start_inds = round(trialstart_plx'.*1000);

sacc_inds=sacc_inds';
ETtrace_raw=ET_trace_raw_1khz';
ETtrace=ET_trace';
ETtrace_plex_calib = PlexET_ad_calib';
ETgains = [g_strctEyeCalib.GainX.Buffer,g_strctEyeCalib.GainY.Buffer];
useLeye=useLeye';
useReye=useReye';

switch targ_stimtype
    case 8
        blockID = BlockID;
        trialID = TrialID;
        cloud_scale = cloud_scale';
        cloud_binary = cloud_binary';

        targchans=find(sum(binned_SU1)>2000);
    case 6
        hartley_metas = hartleystim_metas';
        load(cur_filename_targchans)
end

targchans_SU=targchans(targchans<=nSU);
%targchans=1:nSU;

SU_clusters=[];
Robs=binned_SU1(:,targchans_SU)';
Robs_probe_ID=spk_channels_SU(targchans_SU');
Robs_rating=spk_rating_SU(targchans_SU');

%Robs_sort_inds=[1, 3, 4, 6, 7, 52, 67, 77, 81, 109, 112, 139, 152, 161, 163,171,176,188,190,191,194,194,205,214,215,217,218,222,224,225,226,235,236,240,241,244,249,252];
datafilts=ones(size(Robs));
%
%targchansMU=[4 5 10 28 47 72 80 83 84 86 95 97 98 107 114 119]; %for 220209
%targchansMU=find(sum(binned_SU2)>500);
%RobsMU=binned_SU2(:,targchansMU)';

%targchansMU=[1:nMU]+nSU;
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

switch targ_stimtype
    case 8
        if ~skipET
            if targ_ETstimtype==1;
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
                'stim', 'stimET', 'stimETori', 'stimtype', 'stimtypeET', 'stimscale','valid_data', 'cloud_scale', 'cloud_binary', 'block_inds', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
            else
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot','fixdotET',...
                'stim', 'stimET', 'stimtype', 'stimtypeET', 'stimscale','valid_data', 'cloud_scale', 'cloud_binary', 'block_inds', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
            end
        else
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
                'stim', 'stimtype', 'stimscale','valid_data', 'block_inds', 'cloud_scale', 'cloud_binary', 'blockID', 'trialID', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETtrace_plex_calib', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
        end
        save(cur_filename_targchans, 'targchans')
    case 6
        if ~skipET
            if targ_ETstimtype==1;
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
                'stim', 'stimET', 'stimETori', 'stimtype', 'hartley_metas', 'stimtypeET', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
            else
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
                'stim', 'stimET', 'stimtype', 'hartley_metas', 'stimtypeET', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw',  'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
            end
        else
            save(cur_filename, 'exptname', 'exptdate', 'dt','electrode_info','stim_location','ETstim_location', 'pixel_size', 'fix_location', 'fix_size','fixdot',...
                'stim', 'stimtype', 'hartley_metas', 'stimscale','valid_data', 'block_inds', 'sacc_inds', 'ETtrace', 'ETtrace_raw', 'ETgains', 'trial_start_ts', 'useLeye', 'useReye',...
                'Robs', 'Robs_probe_ID', 'Robs_rating', 'spk_times', 'spk_IDs', 'datafilts', 'RobsMU',  'RobsMU_probe_ID', 'RobsMU_rating', 'datafiltsMU', '-v7.3' )
        end
end
disp('Done with model data') 

%%
disp('saving LFPs...') 

if ~skipLFP
    switch targ_stimtype
        case 8
            cur_filename_LFP=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_v08_LFP.mat'];
                save(cur_filename_LFP,'LFP_ad', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
        case 6
            cur_filename_LFP=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_v08_LFP.mat'];
                save(cur_filename_LFP, 'trial_start_ts', 'trial_start_inds', '-v7.3' )
    end

end
% cur_filename_spikes=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.mat'];
% save(cur_filename_spikes, 'spk_times_all');
% cur_filename_spikes_dat=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.dat'];
% writecell(spk_times_all, cur_filename_spikes_dat);
% cur_filename_spikes_py=['Jocamo_' filenameP(1:6) '_' arraylabel '_' curstimstype '_ET' curETstimtype '_spikes_v07.npy'];
% writeNPY(spk_times_all, cur_filename_spikes_py);
disp('Done with LFPs') 
%%

disp('Done with hartleys!')
close all