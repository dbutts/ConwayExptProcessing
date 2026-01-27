%% Master file
disp('Setup starting')

% make sure paths for dependencies exist
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilosort2/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')

%directory for pregenerated stimulus files
stimpath = '/home/conwaylab/Processing/Cloudstims_calib_01_2022/';

% Switch into data directory
dirpath = '/home/conwaylab/Data/'; %location of the Kofiko trial data folder
pl2path = '/mnt/bc9/Data/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
cd(dirpath)

% this is the name of the experiment you want to run
filenameP = '230510_141725_Jacomo'; 

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')

%% First run the drop test check and see if plexon dropped any measurements
disp('Drop test check starting')
try
    [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] );
end
if n<2
    try
        [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC01'] );
    end
end
[~, n_aux, ts_aux, fn_aux, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['AI01'] );
droptestcheck = [n/40000- n_aux/1000];
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end
disp('Drop test check complete')

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces
disp('Kofiko alignment starting')
which_computer = 3; % (default=2) 2 = LSR, room 2A58

ks.use_online = 1; % set to 1 to use on-line sorting, should be 0 if you want to use kilosort
ks.onlinechans = [1:24]; % which channels of on-line sorted spikes should we go through? 
ks.stitched=0; % if you combined kilosort outputs for multiple arrays
ks.arraylabel ='arrays';
%ks.filepath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.pl2path = pl2path;

opts.eye_tracker = 3; % (default=3) 0=eyescan, 1=monoc eyelink, 2=binoc eyelink, 3=monocular dDPI 
opts.is_cloud = 1; % (default=1) indicates processing for cloud data. set to 0 to skip cloud-specific variables and align task data or other paradigms 
opts.trialwindow = [0 4]; % 4 second trials
if abs(droptestcheck)>0.1;
    opts.spk_offset = droptestcheck; % cTODO: check if this is also true for negative values
else
    opts.spk_offset = 0;
end

[ExptTrials, ExptInfo] = ExtractAndAlignData( filenameP, dirpath, which_computer, ks, opts ); 
% to see options (including a different datadirectory, type 'help ExtractAndAlignData'
% this will automatically generate and save CSD/LFP, but you can call the
% function again and analyze them

% ExtractAndAlign saves FullExpt_ks1_lam_v08.mat in Analysis directory, and
% can be loaded or its outputs (in memory can be used directly
disp('Kofiko alignment complete')
