%%
disp('converting...') 
cd(strExperimentPath)

%%
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
