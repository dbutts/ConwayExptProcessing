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
%dirpath = '/mnt/bc9/Data';
cd(dirpath)

% this is the name of the experiment you want to run
filenameP = '230501_160618_Jacomo'; 

outputdir = [dirpath filenameP '/Analysis/'];
disp('setup complete')
%% first kilosort the data for the laminar probe

preconverted = 0; % 0 if sorting the file for the first time, 1 if the recording has already been converted into dat format

[datPath, droptestcheck] = Step0_KilosortLaminar(dirpath,filenameP,preconverted); %this will take some time to run!

disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present
if abs(droptestcheck)>0.1; warning("Danger - Plexon might have dropped frames! Check pl2 file."); else; disp('Experiment kilosorted and ready for curation!'); end

%% Once this prints "DONE", go curate the file in Phy!

%% now kilosort the data for the other arrays
% TODO

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces

which_computer = 2; % (default=2) 2 = LSR, room 2A58

ks.use_plexon = 1; % set to 1 to use on-line sorting
ks.stitched=0; % if you combined kilosort outputs for multiple arrays
ks.arraylabel ='lam';
ks.FilePath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; 

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
target_SUs = [1 2 4 5 6 9 11];
apply_ETshifts = 1;
stim_deltas = 0;

stas = get_sta(data, target_SUs, apply_ETshifts, stim_deltas)

%% now we plot the best lag of our STAs by depth
%TODO: add stas_depth code that shows select STAs in order of laminar
%arrangement
%%

