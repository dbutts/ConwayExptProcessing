%% Master file
clear all

disp('Setup starting')

bizon = true;

if bizon
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/iCSD/')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilosort2/')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK/')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/kilo2Tools-master/')
    addpath('/home/bizon/Git/ConwayExptProcessing')
    addpath('/home/bizon/Git/ConwayExptProcessing/Tools')
    addpath('/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')
    
    % Switch into data directory
    dirpath = '/home/bizon/Data/V1_Fovea/';
    pl2path = '/home/bizon/Data/V1_Fovea/'; %'/mnt/bc9/Data/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
    stimpath = '/home/bizon/Processing/Cloudstims_calib_04_2024/';
else
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/')
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilosort2/')
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK/')
    addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')

    % Switch into data directory
    dirpath = '/home/conwaylab/Data/';
    pl2path = '/mnt/isilon/DATA/monkey_ephys/Jocamo/2024_Singleprobe/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
    stimpath = '/home/conwaylab/Processing/Cloudstims_calib_04_2024/';
end

cd(dirpath)

% this is the name of the experiment you want to run
filenameP = '250529_152043_Jacomo';
monkey_name = 'Jocamo';

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')
%% first kilosort the data for the laminar probe
disp('Kilosorting Starting')

opts.preconverted = 0; % 0 if sorting the file for the first time, 1 if the recording has already been converted into dat format
opts.nChans = 64; % laminar probes will usually have 24 channels, but this can be set to 32 or 64 depending on the probe
opts.ChnOffset=0;
opts.batch_size = 5e9; % batch size prevents memory overflow errors. default = 300000000 for IT
opts.monkey_name = monkey_name;

[datPath, droptestcheck] = Step0_KilosortLaminar(dirpath,filenameP,pl2path,opts); %this will take some time to run!

disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end
disp('Kilosorting complete, go curate in phy')


disp('In terminal type:')
disp('conda activate phy2')
disp(['cd ' dirpath filenameP filesep 'kilosorting_laminar' filesep])
disp('phy template-gui params.py')

%% Once this prints "DONE", go curate the file in Phy!

%% Once this prints "DONE", go curate the file in Phy!
%{
%% now kilosort the data for the other arrays
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
    disp('Done with arrays!! Now ready to combine kilosorting files.')

    %% to combine kilosort outputs for multiple arrays

    strStitchPath = [dirpath filenameP '/kilosorting_stitched/'];
    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_laminar/'],0);

%    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_1to32/'],160);
    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_1to32/'],160, spk_info, spk_times, spk_clusters);
    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_33to64/'],192, spk_info, spk_times, spk_clusters);
    [spk_info, spk_times, spk_clusters] = fn_kiloappend([dirpath filenameP '/kilosorting_UT2_65to96/'],224, spk_info, spk_times, spk_clusters);

    if ~exist(strStitchPath,'dir');
        mkdir(strStitchPath);
    end
    save([strStitchPath 'KS_stitched.mat'], "spk_clusters", "spk_times", "spk_info")

    disp('saving ks_stich - done')
%}    
%% if you are rerunning the stuff below, uncomment this cell to quickly rerun the drop test check and see if plexon dropped any measurements
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
droptestcheck = [n_aux/1000-n/40000];
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end
disp('Drop test check complete')

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces
disp('Kofiko alignment starting')
which_computer = 2; % (default=2) 2 = Caneron's Bizon

ks.use_online = 0; % set to 1 to use on-line sorting or no sorting at all. Should be 0 if you want to use kilosort
ks.onlinechans = [1:64]; % which channels of on-line sorted spikes should we go through? 

ks.stitched=0; % if you combined kilosort outputs for multiple arrays
ks.arraylabel ='lam';
ks.filepath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
%ks.filepath = [dirpath filenameP filesep 'kilosorting_stitched' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.pl2path = pl2path;

opts.eye_tracker = 4;   % (default=3) 0=eyescan, 1=monoc eyelink, 2=binoc eyelink, 3=monocular dDPI, 4=binocular dDPI
opts.is_cloud = 1;      % (default=1) indicates processing for cloud data. set to 0 to skip cloud-specific variables and align task data or other paradigms 
opts.extractfixinfo = 1;% (default=1) extracts fixinfo file containing fivedot and dotgrid info
opts.trialwindow = [0 4]; % 4 second trials
opts.trl_fix_thresh = 0.6; % include trials with at least 3 seconds of fixation
opts.monkey_name = monkey_name;

if abs(droptestcheck)>0.1;
    opts.spk_offset = droptestcheck;
else
    opts.spk_offset = 0;
end

[ExptTrials, ExptInfo, ETdata] = ExtractAndAlignData( filenameP, dirpath, which_computer, ks, opts ); 
% to see options (including a different datadirectory, type 'help ExtractAndAlignData'
% this will automatically generate and save CSD/LFP, but you can call the
% function again and analyze them

% ExtractAndAlign saves FullExpt_ks1_lam_v08.mat in Analysis directory, and
% can be loaded or its outputs (in memory can be used directly
% This step also makes the fixinfo file
disp('Kofiko alignment complete')

%% Package cloud data
disp('Cloud packaging starting')
% now package cloud data
ExptInfo.trialdur = 4;
ExptInfo.monkey_name = monkey_name;
data = PackageCloudData_v9( ExptTrials, ExptInfo, [], [], stimpath, outputdir );
% data = PackageCloudData_v9( ExptTrials, ExptInfo );

disp('Cloud packaging complete')
%% Check to see if all of the relevant files were saved
checkOutput(outputdir,filenameP)

%% finally, generate STAs for all recorded cells
disp('STA generation starting')
stas = get_sta(data);
disp('STA generation complete')

brake
 %% if you want to generate STAs for only a subset of the data, use this code instead:
% disp('STA subset generation starting')
% n_all = size(data.Robs, 1)+size(data.RobsMU, 1);
% target_SUs = [35:n_all]; % Change selection to a specific number of SUs Total number of SUs = size(data.Robs, 1)+size(data.RobsMU, 1)
% apply_ETshifts = 0;
% stim_deltas = ones(size(data.ETtrace))'; 
% % stim_deltas(1,:) = stim_deltas(1,:).*10;
% % stim_deltas(2,:) = stim_deltas(2,:).*20;
% num_lags = 6;
% %stim_deltas = data.stim_location_deltas;
% %use_inds = data.valid_data;
% use_inds = intersect(find(data.cloud_binary<2), data.valid_data); % 0 = full 1 = matched 2 = hybrid
%     use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts
% 
% save_vars.to_save = 0; % saves the STAs as pdf when set to 1
% save_vars.outputdir = outputdir;
% %save_vars.titlestr = [filenameP '_Plexon_Cloud60_']; 
% save_vars.titlestr = [filenameP '_FullClouds_']; 
% 
% stas = get_sta(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', num_lags, save_vars);
% disp('STA subset generation complete')

%% and if you want to run a forward correlation on the same info, set up the info in the same way
% disp('Forward Correlation subset generation starting')
% n_all = size(data.Robs, 1)+size(data.RobsMU, 1);
% target_SUs = [9]; %[1:n_all]; % Change selection to a specific number of SUs Total number of SUs = size(data.Robs, 1)+size(data.RobsMU, 1)
% apply_ETshifts = 0;
% stim_deltas = ones(size(data.ETtrace))'; 
% % stim_deltas(1,:) = stim_deltas(1,:).*10;
% % stim_deltas(2,:) = stim_deltas(2,:).*20;
% num_lags = 6;
% thresh=0.15;
% %stim_deltas = data.stim_location_deltas;
% %use_inds = data.valid_data;
% use_inds = intersect(find(data.cloud_binary==2), data.valid_data);
%     use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts
% 
% save_vars.to_save = 1; % saves the STAs as pdf when set to 1
% save_vars.outputdir = outputdir;
% %save_vars.titlestr = [filenameP '_BaseClouds_']; 
% save_vars.titlestr = [filenameP '_MatchedClouds_v2_']; 
% 
% [fwdc_L, fwdc_M, fwdc_diff] = get_fwdc(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', thresh, num_lags, save_vars);
% disp('Forward Correlation subset generation complete')

%% next, package the data, starting with Hartleys
% disp('Data packaging starting')
% %load([outputdir filenameP '_FullExpt_ks1_lam_v09.mat'])
% 
% % first package the hartley stimuli
% data_hartleys = PackageCloudData_v9( ExptTrials, ExptInfo, 6, [], stimpath, outputdir );
% disp('Data packaging complete')

%% Generate hartleys 
% disp('Hartley generation starting')
% % get_hartleys(input data, SUs we want to plot, lag at which we want to plot)
% % to plot all SUs target_SUs = size(data_hartleys.Robs, 1)
% % to plot specific SUs target_SUs = [ 1 33 45 ]
% target_SUs = [1:67];
% lag = 4;
% save_vars.to_save = 1; % saves the STAs as pdf when set to 1
% save_vars.outputdir = outputdir;
% save_vars.titlestr = filenameP; 
% 
% % data=[]; % if you want to skip plotting STA info alongside hartley tuning plots
% get_hartleys(data_hartleys, target_SUs, lag, save_vars, data);
% 
% disp('Hartley generation complete')
%% now we plot the best lag of our STAs by depth
%TODO: add stas_depth code that shows select STAs in order of laminar
%arrangement

%% Move files to server (very slow, probably better to do manually)
% remoteBase = '/run/user/1000/gvfs/smb-share:server=lsr-isilon.nei.nih.gov,share=lsr-conway';
% remotePath = fullfile(remoteBase, 'projects', 'V1_Fovea', 'processing', filenameP);
% copyfile(outputdir, [remotePath '/Analysis']);
% copyfile([dirpath filenameP '/kilosorting_laminar'],[remotePath '/kilosorting_laminar'])


%% Delete large temporary files
disp('Deleting temporary files started')
dirpath2 = ['/home/bizon/Data/V1_Fovea/',filenameP,'/kilosorting_laminar/'];
cd(dirpath2)
dat = [filenameP, '.dat'];
% Move .dat, SampToSecsMap.mat, and res2.mat file to trash
delete (dat) 
delete sampsToSecsMap.mat rez2.mat
disp('Deleting temporary files complete')
