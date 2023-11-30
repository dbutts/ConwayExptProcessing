function [saccade_inds, sac_start_inds, sac_stop_inds, saccade_times, sac_start_times, sac_stop_times] = detect_saccades_bevfel(ET_trace, ET_times, upsample, eye_smooth, sgolay_deg, sac_thresh, peri_thresh, plot_flag)
%detect_saccades_bevfel Identify saccade indices in Bevil's Rig C data
%   Detailed explanation goes here

%% Parameters
%     ET_times=[0:.001:3.999];
%     upsample=3;
%     eye_smooth = [81, 3, 3]; sgolay_deg = [2,2,2];
%     sac_thresh = 6; %3.5; %threshold eye speed
%     peri_thresh = 2.5; %threshold eye speed for defining saccade boundary inds
  
eye_dt=1; 
min_isi = 0.15; max_isi = inf; %min/max inter-saccade intervals

%for amplitude correction
amp_cutoff=6; %arcmins combined for X and Y

% win=[3 8];
%% upsample traces
if upsample==1 % upsample 60Hz kofiko traces to to 1000Hz
    et_params.eye_fs = 1000;
    ET_times_up = ET_times(1):(1/et_params.eye_fs):ET_times(end);
    ET_ad_up(1,:) = interp1(ET_times, ET_trace(1,:), ET_times_up, 'spline');
    ET_ad_up(2,:) = interp1(ET_times, ET_trace(2,:), ET_times_up, 'spline');
    
    ET_ad_up_s(1,:) = smooth(ET_ad_up(1,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
    ET_ad_up_s(2,:) = smooth(ET_ad_up(2,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
    sm_avg_eyepos = ET_ad_up_s';

elseif upsample==2 % keep high srate of plexon traces
    et_params.eye_fs = 1000; %ET_adfreq;
    ET_times_up = ET_times(1):(1/et_params.eye_fs):ET_times(end);
    ET_ad_up_s(1,:) = smooth(ET_trace(1,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
    ET_ad_up_s(2,:) = smooth(ET_trace(2,:), eye_smooth(1), 'sgolay', sgolay_deg(1));

    sm_avg_eyepos = ET_ad_up_s';

elseif upsample==3 % keep high srate of plexon traces for smoothing, then downsample to stim resolution
    ET_ad_up(1,:) = smooth(ET_trace(1,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
    ET_ad_up(2,:) = smooth(ET_trace(2,:), eye_smooth(1), 'sgolay', sgolay_deg(1));

    et_params.eye_fs = 60; %ET_adfreq;
    ET_times_up = ET_times(1):(1/et_params.eye_fs):ET_times(end);
    ET_ad_up_s(1,:) = interp1(ET_times, ET_ad_up(1,:), ET_times_up, 'spline');
    ET_ad_up_s(2,:) = interp1(ET_times, ET_ad_up(2,:), ET_times_up, 'spline'); 
 
    ET_ad_up_s(1,:) = smooth(ET_ad_up_s(1,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
    ET_ad_up_s(2,:) = smooth(ET_ad_up_s(2,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
    sm_avg_eyepos = ET_ad_up_s';
%     dbstop if error
%     error stop

elseif upsample==4 % use low-pass filter, then downsample
%     ET_ad_up(1,:) = smooth(ET_trace(1,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
%     ET_ad_up(2,:) = smooth(ET_trace(2,:), eye_smooth(1), 'sgolay', sgolay_deg(1));
% ET_ad_up(1,:) = ET_trace(1,:);
% ET_ad_up(2,:) = ET_trace(2,:);

    passband_fs=120;    
    ET_ad_up(1,:) = lowpass(ET_trace(1,:), passband_fs, 1000);
    ET_ad_up(2,:) = lowpass(ET_trace(2,:), passband_fs, 1000);


    et_params.eye_fs = 60; %ET_adfreq;
    ET_times_up = ET_times(1):(1/et_params.eye_fs):ET_times(end);
    ET_ad_up_s(1,:) = interp1(ET_times, ET_ad_up(1,:), ET_times_up, 'spline');
    ET_ad_up_s(2,:) = interp1(ET_times, ET_ad_up(2,:), ET_times_up, 'spline'); 
 
    ET_ad_up_s(1,:) = smooth(ET_ad_up_s(1,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
    ET_ad_up_s(2,:) = smooth(ET_ad_up_s(2,:), eye_smooth(2), 'sgolay', sgolay_deg(2));
    sm_avg_eyepos = ET_ad_up_s';
%     dbstop if error
%     error stop

else %use kofiko eyetracking samples at 60hz
    et_params.eye_fs = 60; %ET_adfreq;
    sm_avg_eyepos = ET_trace';
    sm_avg_eyepos(:,1) = smooth(sm_avg_eyepos(:,1),eye_smooth(1));
    sm_avg_eyepos(:,2) = smooth(sm_avg_eyepos(:,2),eye_smooth(1));
end

%% SACCADE DETECTION (from detect_saccades_v2
%disp('Detecting saccades');

eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]./eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]./eye_dt;
all_eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

if length(eye_smooth)>2; %smooth eye velocity for better peri-saccade detection?
    all_eye_speed = smooth(all_eye_speed, eye_smooth(3), 'sgolay', sgolay_deg(3));
end
ET_n=length(all_eye_speed);
% %
%parameters
isi_cutoff=min_isi*et_params.eye_fs;

%find local maxima of eye speed signal exceeding sac_thresh
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);
%saccade_inds(saccade_inds<isi_cutoff)=[];
%saccade_inds(saccade_inds>(ET_n-isi_cutoff))=[];

%saccade_inds = [1; saccade_inds; ET_n];

%%
%find times when speed signal crossed above and below the peri-saccade
%threshold
thresh_cross_up = [1; 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh)];
thresh_cross_down = [1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh); ET_n];
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
%sac_start_inds(1) = 1;
%sac_stop_inds(end)=ET_n;

%% enforce inter--saccade interval to get rid of double-peaks
%/{
isis = [Inf; diff(sac_start_inds)]/et_params.eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
%}
%{
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
%}

%% kludge: identify artifacts (saccade stayed above peri_thresh for too long)
% bad_sac_inds=find(diff(sac_start_inds)==0)+1;
% sac_start_inds(bad_sac_inds)=[];
% sac_stop_inds(bad_sac_inds)=[];
% saccade_inds(bad_sac_inds)=[];
% 
% bad_sac_inds=find(saccade_inds>[ET_n-win(2)-1]);
% sac_start_inds(bad_sac_inds)=[];
% sac_stop_inds(bad_sac_inds)=[];
% saccade_inds(bad_sac_inds)=[];

%% correct for amplitude
%{
nsacs=length(saccade_inds);
ssi=1; cur_ssi=2;
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
%/{
nsacs=length(saccade_inds);
ssi=1; cur_ssi=2;
while ssi<nsacs-2
    cur_preampX=median(sm_avg_eyepos(saccade_inds(cur_ssi-1):saccade_inds(cur_ssi),1));
    cur_postampX=median(sm_avg_eyepos(saccade_inds(cur_ssi):saccade_inds(cur_ssi+1),1));
    cur_preampY=median(sm_avg_eyepos(saccade_inds(cur_ssi-1):saccade_inds(cur_ssi),2));
    cur_postampY=median(sm_avg_eyepos(saccade_inds(cur_ssi):saccade_inds(cur_ssi+1),2));
    if [abs(cur_postampX-cur_preampX)+abs(cur_postampY-cur_preampY)] < amp_cutoff;
        saccade_inds(cur_ssi)=[];
        sac_start_inds(cur_ssi) = []; 
        sac_stop_inds(cur_ssi) = [];
        ssi=ssi+1;
    else
        cur_ssi=cur_ssi+1;
        ssi=ssi+1;
    end
%     [sac_start_inds,sac_stop_inds]
%     pause

end
%}


%%
if plot_flag>0
    dbstop if error
    figure(plot_flag);
    plot(sm_avg_eyepos); hold on
    plot(all_eye_speed, 'k'); 
    if ~isempty(saccade_inds)
    vline(sac_start_inds, 'k')
    vline(sac_stop_inds, 'r')
    end
    %xlabel('Seconds'); 
    ylabel('arcmin'); hold off
end

if upsample>=1
saccade_times = ET_times_up(saccade_inds); %saccade peak times
sac_start_times = ET_times_up(sac_start_inds); %saccade start times
sac_stop_times = ET_times_up(sac_stop_inds); %saccade end times
else
%    [sac_start_inds, sac_stop_inds];
saccade_times = ET_times(saccade_inds); %saccade peak times
sac_start_times = ET_times(sac_start_inds); %saccade start times
sac_stop_times = ET_times(sac_stop_inds); %saccade end times    
end

end