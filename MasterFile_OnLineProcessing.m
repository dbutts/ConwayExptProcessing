%% Master file
disp('Setup starting')

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
filenameP = '230510_141725_Jacomo'; 

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')

%% First run the drop test check and see if plexon dropped any measurements
disp('Drop test check starting')
[fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] );
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
%% next, package the data
disp('Data packaging starting')
%load([outputdir filenameP '_FullExpt_ks1_lam_v09.mat'])

% first package the hartley stimuli
data_hartleys = PackageCloudData_v9( ExptTrials, ExptInfo, 6, [], stimpath, outputdir );
disp('Data packaging complete')

%% Generate hartleys 
disp('Hartley generation starting')
% get_hartleys(input data, SUs we want to plot, lag at which we want to
% plot)
% to plot all SUs target_SUs = size(data_hartleys.Robs, 1)
% to plot specific SUs target_SUs = [ 1 33 45 ]
target_SUs = [5];
lag = 3;
save_vars.to_save = 1; % saves the STAs as pdf when set to 1
save_vars.outputdir = outputdir;
save_vars.titlestr = filenameP; 


get_hartleys(data_hartleys, target_SUs, lag, save_vars, data);

disp('Hartley generation complete')
%% Package cloud data
disp('Cloud packaging starting')
% now package cloud data
ExptInfo.trialdur = 4;
data = PackageCloudData_v9( ExptTrials, ExptInfo, [], [], stimpath, outputdir );
% data = PackageCloudData_v9( ExptTrials, ExptInfo );

disp('Cloud packaging complete')

%% finally, generate STAs for all recorded cells
disp('STA generation starting')
stas = get_sta(data);
disp('STA generation complete')

%% if you want to generate STAs for only a subset of the data, use this code instead:
disp('STA subset generation starting')
target_SUs = [5]; % Change selection to a specific number of SUs Total number of SUs = size(data.Robs, 1)
apply_ETshifts = 1;
stim_deltas = 0; 
num_lags = 6;
%stim_deltas = data.stim_location_deltas;
use_inds = intersect(find(data.cloud_area==60), data.valid_data);
    use_inds(end-num_lags:end) = []; %cut last few indices to avoid artifacts

save_vars.to_save = 1; % saves the STAs as pdf when set to 1
save_vars.outputdir = outputdir;
%save_vars.titlestr = [filenameP '_Plexon_Cloud60_']; 
save_vars.titlestr = [filenameP '_Cloud60_']; 

stas = get_sta(data, target_SUs, use_inds, apply_ETshifts, stim_deltas', num_lags, save_vars);
disp('STA subset generation complete')
%% now we plot the best lag of our STAs by depth
%TODO: add stas_depth code that shows select STAs in order of laminar
%arrangement

%% Delete large temporary files
disp('Deleting temporary files started')
dirpath2 = ['/home/conwaylab/Data/',filenameP,'/kilosorting_laminar/'];
cd(dirpath2)
dat = [filenameP, '.dat'];
% Move .dat, SampToSecsMap.mat, and res2.mat file to trash
delete (dat) 
delete sampsToSecsMap.mat rez2.mat
disp('Deleting temporary files complete')
