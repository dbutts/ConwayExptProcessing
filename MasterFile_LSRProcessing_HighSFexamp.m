%% Master file

addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilosort2/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK/')
addpath('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')

%directory for pregenerated stimulus files
stimpath = '/home/conwaylab/Processing/Cloudstims_calib_01_2022/';

% Switch into data directory
dirpath = '/home/conwaylab/Data/';
pl2path = '/mnt/bc9/Data/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
cd(dirpath)

% this is the name of the experiment you want to run
%filenameP = '220620_131047_Jacomo'; 
filenameP = '230501_160618_Jacomo'; 

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')
%% first kilosort the data for the laminar probe

preconverted = 0; % 0 if sorting the file for the first time, 1 if the recording has already been converted into dat format

[datPath, droptestcheck] = Step0_KilosortLaminar(dirpath,filenameP,pl2path,preconverted); %this will take some time to run!

disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end

%% Once this prints "DONE", go curate the file in Phy!

%% now kilosort the data for the other arrays
% TODO

%% if you are rerunning the stuff below, uncomment this cell to quickly rerun the drop test check and see if plexon dropped any measurements
%/{
[fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] );
[~, n_aux, ts_aux, fn_aux, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['AI01'] );
droptestcheck = [n/40000- n_aux/1000];
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end
%}

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces

which_computer = 2; % (default=2) 2 = LSR, room 2A58

ks.use_online = 0; % set to 1 to use on-line sorting
ks.onlinechans = [1:24]; % which channels of on-line sorted spikes should we go through? 
ks.stitched=0; % if you combined kilosort outputs for multiple arrays
ks.arraylabel ='lam';
ks.filepath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; 
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

%% next, package the data
%load([outputdir filenameP '_FullExpt_ks1_lam_v09.mat'])

% first package the hartley stimuli
data_hartleys = PackageCloudData_v9( ExptTrials, ExptInfo, 6, [], stimpath, outputdir );

% now package cloud data
ExptInfo.trialdur = 4;
data = PackageCloudData_v9( ExptTrials, ExptInfo, [], [], stimpath, outputdir );
% data = PackageCloudData_v9( ExptTrials, ExptInfo );

disp('All Done')

%% finally, generate STAs for all recorded cells
stas = get_sta(data);
disp('STAs finished!')

%% if you want to generate STAs for only a subset of the data, use this code instead:
target_SUs = [1:63];
apply_ETshifts = 1;
stim_deltas = 0; 
num_lags = 6;
%stim_deltas = data.stim_location_deltas;
use_inds = intersect(find(data.cloud_area==180), data.valid_data);
    use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts

save_vars.to_save = 1; % saves the STAs as pdf when set to 1
save_vars.outputdir = outputdir;
save_vars.titlestr = [filenameP '_Plexon_Cloud180_']; 

stas = get_sta(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', num_lags, save_vars)

%% compare firing rates
use_inds180 = intersect(find(data.cloud_area==180), data.valid_data);
use_inds60 = intersect(find(data.cloud_area==60), data.valid_data);
FRs(1,:) = mean(data.Robs(:,use_inds60)');
FRs(2,:) = mean(data.Robs(:,use_inds180)');
figure; plot(FRs'); legend({'Cloud size 60', 'Cloud size 180'})

%% now we plot the best lag of our STAs by depth
%TODO: add stas_depth code that shows select STAs in order of laminar
%arrangement
%%

