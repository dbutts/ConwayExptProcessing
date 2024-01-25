function [ETtrace] = fnAlignEyeTrace(ETfilename, plx_trial_ts, trial_bins)
%
% Usage: ETtrace = fnAlignEyeTrace(ETfilename, plx_trial_ts, <trial_bins>)
% 
% Takes full ET trace and align relevant sections to trial onsets
% 
% ex ETfilepath: ('/home/conwaylab/Data/230510_141725_Jacomo/Analysis/230510_141725_Jacomo_FullExpt_ET.mat')
% plx_trial_ts is currently saved in PackageCloudData as trial_start_ts

if nargin < 3
    trial_bins = 4000; % default to 4 second trials
end

load(ETfilename);
num_trials = length(plx_trial_ts);
num_ETchans = size(PlexET_ad_calib, 2); %

if exist('DDPI_trace','var')
    ETtrace=zeros([num_trials*trial_bins,2]); % DDPI only works on one eye, so two traces (V and H) 
    % but we could extract additional channels for other variables like pupil size and include them here too
   
    % DDPI signal processing goes here
else
    % if there is no DDPI data, we just use the Plexon signals
    ETtrace=zeros([num_trials*trial_bins,num_ETchans]);

    for tt = 1:num_trials
        cur_trial_start = find(PlexET_times>=plx_trial_ts(tt),1,'first');
        cur_inds = cur_trial_start:cur_trial_start+trial_bins-1;
        ETtrace([1:trial_bins]+trial_bins*(tt-1),:) = PlexET_ad_calib(cur_inds,:);
    end
end
end