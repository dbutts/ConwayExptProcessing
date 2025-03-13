%% Master file
clear all

disp('Setup starting')

addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/iCSD/')
addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilosort2/')
addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK/')
% addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')
addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/kilo2Tools-master/')
addpath('/home/bizon/Git/ConwayExptProcessing')
addpath('/home/bizon/Git/ConwayExptProcessing/Tools')

%directory for pregenerated stimulus files
stimpath = '/home/bizon/Processing/Cloudstims_calib_04_2024/';

% Switch into data directory
dirpath = '/home/bizon/Data/';
pl2path = '/home/bizon/Data/'; %'/mnt/bc9/Data/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
cd(dirpath)

% this is the name of the experiment you want to run
filenameP = '250303_175423_Jacomo';
monkey_name = 'Jacamo';

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')
%% first kilosort the data for the laminar probe
disp('Kilosorting Starting')

opts.preconverted = 0; % 0 if sorting the file for the first time, 1 if the recording has already been converted into dat format
opts.nChans = 64; % laminar probes will usually have 24 channels, but this can be set to 32 or 64 depending on the probe
opts.ChnOffset=0;
opts.batch_size = 5e9; % batch size prevents memory overflow errors. default = 300000000 for IT
opts.monkey_name = monkey_name;

%
[datPath, droptestcheck] = Step0_KilosortLaminar(dirpath,filenameP,pl2path,opts); %this will take some time to run! Packages and then kilosorts laminar probe data
% 
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted go curate in phy!'); end
disp('Kilosorting complete')

disp('In terminal type:')
disp('conda activate phy2')
disp(['cd ' dirpath filenameP filesep 'kilosorting_laminar' filesep])
disp('phy template-gui params.py')

%% Once this prints "DONE", go curate the file in Phy!

%% If you have multiple arrays, now kilosort the data for the other arrays
    % %% Utah 1 - Serial 0071
    % opts.ArrayLabel = 'UT1'; %load channel map
    % opts.chInfo = load('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/Vinny/V_UT1_chanMap_nogroup.mat');
    % opts.ChnOffset=64;
    % 
    % opts.curchannels = [1:32]; outputFolder = [dirpath filenameP '/kilosorting_UT1_1to32/'];
    % [datPaths.UT1a, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
    % 
    % opts.curchannels = [33:64]; outputFolder = [dirpath filenameP '/kilosorting_UT1_33to64/'];
    % [datPaths.UT1b, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
    % 
    % opts.curchannels = [65:96]; outputFolder = [dirpath filenameP '/kilosorting_UT1_65to96/'];
    % [datPaths.UT1c, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
    %% Utah 2 - Serial 0072
    % opts.ArrayLabel = 'UT2'; %load channel map
    % opts.chInfo = load('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/Vinny/V_UT2_chanMap_nogroup.mat');
    % opts.ChnOffset=160;
    % 
    % opts.curchannels = [1:32]; outputFolder = [dirpath filenameP '/kilosorting_UT2_1to32/'];
    % [datPaths.UT2a, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
    % 
    % opts.curchannels = [33:64]; outputFolder = [dirpath filenameP '/kilosorting_UT2_33to64/'];
    % [datPaths.UT2b, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
    % 
    % opts.curchannels = [65:96]; outputFolder = [dirpath filenameP '/kilosorting_UT2_65to96/'];
    % [datPaths.UT2c, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);

     %% NForm array
%     opts.ArrayLabel = 'NF'; %load channel map
%     opts.chInfo = load('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/Vinny/V_NF_chanMap.mat');
%     opts.ChnOffset=256;
% 
%     opts.curchannels = [1:32]; outputFolder = [dirpath filenameP '/kilosorting_NF_1to32/'];
%     [datPaths.UT1, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
% 
%     opts.curchannels = [33:64]; outputFolder = [dirpath filenameP '/kilosorting_NF_33to64/'];
%     [datPaths.UT1, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);
% 
%     opts.curchannels = [65:96]; outputFolder = [dirpath filenameP '/kilosorting_NF_65to96/'];
%     [datPaths.UT1, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);    
% 
%     opts.curchannels = [97:128]; outputFolder = [dirpath filenameP '/kilosorting_NF_97to128/'];
%     [datPaths.UT1, ~] = Step0_KilosortArray(dirpath,filenameP,pl2path,outputFolder,opts);   
    %%
    % disp('Done with arrays!! Now ready to combine kilosorting files.')

    %% If you have multiple arrays, combine kilosort outputs for multiple arrays

%     strStitchPath = [dirpath filenameP '/kilosorting_stitched/'];
%     [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_laminar/'],0);
% 
% %    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_1to32/'],160);
%     [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_1to32/'],160, spk_info, spk_times, spk_clusters);
%     [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_33to64/'],192, spk_info, spk_times, spk_clusters);
%     [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_65to96/'],224, spk_info, spk_times, spk_clusters);
% 
%     if ~exist(strStitchPath,'dir');
%         mkdir(strStitchPath);
%     end
%     save([strStitchPath 'KS_stitched.mat'], "spk_clusters", "spk_times", "spk_info")
% 
%     disp('saving ks_stich - done')
%% if you are rerunning the stuff below, uncomment this cell to quickly rerun the drop test check and see if plexon dropped any measurements
disp('Drop test check starting')
try
    [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] ); % raw voltages probes with >100 channels 
end
if n<2
    try
        [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC01'] ); % raw voltages probes with <100 channels
    end
end

[~, n_aux, ts_aux, fn_aux, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['AI01'] ); % AIO1 is analog input for strobes / LFPs
droptestcheck = [n/40000- n_aux/1000];
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end
disp('Drop test check complete')

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces
disp('Kofiko alignment starting')
which_computer = 4; % 4 = Cameron's bizon, other options within ExtractAndAlignData.m


ks.use_online = 1; % set to 1 to use on-line sorting, should be 0 if you want to use kilosort
ks.onlinechans = [1:64]; % which channels of on-line sorted spikes should we go through? 

ks.stitched = 0; % if you combined kilosort outputs for multiple arrays 0=one array/probe 1=more than one array/probe
ks.arraylabel ='laminar'; % ='laminar or ='Utah'; goes into metadata
ks.filepath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
% ks.filepath = [dirpath filenameP filesep 'kilosorting_stitched' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.pl2path = pl2path;

opts.eye_tracker = 4; % (default=3) 0=eyescan, 1=monoc eyelink, 2=binoc eyelink, 3=monocular dDPI, 4=binocular dDPI 
opts.is_cloud = 0; % (default=1) indicates processing for cloud data. set to 0 to skip cloud-specific variables and align task data or other paradigms. cloud data processing generates iCSDs and LFPs
opts.trialwindow = [0 4]; % change trial window for both clouds and mturk - cloud trial window [0 4] - mturk1 trial window [-0.5 6]
opts.trl_fix_thresh = 2/(opts.trialwindow(2)-opts.trialwindow(1)); % include trials with at least 2 seconds of fixation - flag is used in package data step
opts.monkey_name = monkey_name;

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
%% next, package the data
% disp('Data packaging starting')
% %load([outputdir filenameP '_FullExpt_ks1_lam_v09.mat'])
% 
% % first package the hartley stimuli
% data_hartleys = PackageCloudData_v9( ExptTrials, ExptInfo, 6, [], stimpath, outputdir );
% disp('Data packaging complete')

%% Generate hartleys 
% disp('Hartley generation starting')
% % get_hartleys(input data, SUs we want to plot, lag at which we want to
% % plot)
% % to plot all SUs target_SUs = size(data_hartleys.Robs, 1)
% % to plot specific SUs target_SUs = [ 1 33 45 ]
% target_SUs = [1:5];
% lag = 4;
% save_vars.to_save = 1; % saves the STAs as pdf when set to 1
% save_vars.outputdir = outputdir;
% save_vars.titlestr = filenameP; 
% 
% get_hartleys(data_hartleys, target_SUs, lag, save_vars, data);
% 
% disp('Hartley generation complete')
%% Package cloud data
disp('Cloud packaging starting')
% now package cloud data
ExptInfo.trialdur = 4;
ExptInfo.monkey_name = monkey_name;
data = PackageCloudData_v9( ExptTrials, ExptInfo, [], [], stimpath, outputdir );
% data = PackageCloudData_v9( ExptTrials, ExptInfo );

disp('Cloud packaging complete')

%% finally, generate STAs for all recorded cells
disp('STA generation starting')
stas = get_sta(data);
disp('STA generation complete')

%% if you want to generate STAs for only a subset of the data, use this code instead:
disp('STA subset generation starting')
target_SUs = [3 5 8 10 13 17 18 size(data.Robs, 1)+4 ]; % Change selection to a specific number of SUs Total number of SUs = 1:(size(data.Robs, 1)+size(data.RobsMU, 1)); To show MUs need to do (size(data.Robs, 1) + MU_number
apply_ETshifts = 1; % 1= account for eye position 0 = don't account for eye position
stim_deltas = 0; 
num_lags = 6;
%stim_deltas = data.stim_location_deltas;
use_inds = intersect(find(data.cloud_area==60), data.valid_data);
% %use_inds = intersect(1:48000, data.valid_data); % cut off part of the experiment for sta generation
%     use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts
use_inds = intersect(find(data.cloud_binary<2), data.valid_data); % 0 = full 1 = matched 2 = hybrid
    use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts

save_vars.to_save = 1; % saves the STAs as pdf when set to 1
save_vars.outputdir = outputdir;
%save_vars.titlestr = [filenameP '_Plexon_Cloud60_']; 
save_vars.titlestr = [filenameP '_Cloud60_']; 

stas = get_sta(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', num_lags, save_vars);
disp('STA subset generation complete')

%% plot the L,M,and S of rectangular regions
target_SUs = [9];
rectangle1 = [10 10 12 28]; % [x y w h]
rectangle2 = [2 2 8 8];
rectangle3 = [25 20 20 20];
rectangle4 = [26 19 20 6];
get_rec_opponency(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', num_lags, save_vars, rectangle1, rectangle2, rectangle3, rectangle4);

%% plot the L,M, and S of circular regions
target_SUs = [9];
type1 = ["circle", "r"]; %[type, color]
parameters1 = [35, 28, 5]; %[xcenter ycenter radius]
type2 = ["rectangle","b"];
parameters2 = [10, 2, 30, 10]; % [x y width height]
type3 = ["polygon","g"];
parameters3 = [12 21 17 20 13 8; 15 15 27 37 37 30]; % [x1 x2... ; y1 y2...] clockwise
type4 = ["polygon","y"];
parameters4 = [23 34 48 37 25; 23 17 21 22 27];
save_vars.titlestr = [filenameP '_blue_normalized_response_']; 

get_regional_response(data, stas, target_SUs, apply_ETshifts, num_lags, save_vars, ...
    type1, parameters1, type2, parameters2, type3, parameters3, type4, parameters4);

%% now we plot the best lag of our STAs by depth
%TODO: add stas_depth code that shows select STAs in order of laminar
%arrangement

%% Delete large temporary files
disp('Deleting temporary files started')
dirpath2 = ['/home/bizon/Data/',filenameP,'/kilosorting_laminar/'];
cd(dirpath2)
dat = [filenameP, '.dat'];
% Move .dat, SampToSecsMap.mat, and res2.mat file to trash
delete (dat) 
delete sampsToSecsMap.mat rez2.mat
disp('Deleting temporary files complete')
