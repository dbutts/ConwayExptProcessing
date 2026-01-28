%% Master file
disp('Setup starting')

% make sure paths for dependencies exist
addpath/lsr-conway/code/ConwayExptProcessing/Dependencies/Plexon-Matlab Offline Files SDK')
% these are only necessary for CSD analyses
addpath('/mnt/isilon/code/ConwayExptProcessing/Dependencies/iCSD/')
addpath('/mnt/isilon/code/ConwayExptProcessing/Dependencies/iCSD/CSD_functions/')
%these are only necessary for kilosort
addpath('/mnt/isilon/code/ConwayExptProcessing/Dependencies/Kilosort2/')
addpath('/mnt/isilon/code/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/')
addpath('/mnt/isilon/code/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')


% Switch into data directory
dirpath = '/mnt/isilon/DATA/monkey_ephys/Pollux/CausalGlobs/240913_164554_Pollux/240913_164554_Pollux'; %location of the Kofiko trial data folder
pl2path = '/mnt/isilon/DATA/monkey_ephys/Pollux/CausalGlobs/240913_164554_Pollux'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
cd(dirpath)

% this is the name of the experiment you want to run
filenameP = '240913_164554_Pollux'; 
monkey_name = 'Pollux';

outputdir = ['/home/conwaylab/Data/' filenameP '/'];
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
ks.onlinechans = [1]; % which channels of on-line sorted spikes should we go through? ex: [1:24] for a 24-channel laminar probe. leave as 1 for single units
ks.stitched=0; % if you combined kilosort outputs for multiple arrays
ks.arraylabel ='single';
%ks.filepath = [dirpath filenameP filesep 'kilosorting_laminar' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.pl2path = pl2path;

opts.eye_tracker = 3; % (default=3) 0=eyescan, 1=monoc eyelink, 2=binoc eyelink, 3=monocular dDPI 
opts.is_cloud = 1; % (default=1) indicates processing for cloud data. set to 0 to skip cloud-specific variables and align task data or other paradigms 
opts.trialwindow = [0 4]; % 4 second trials
opts.trl_fix_thresh = 0.9; % include trials where the monkey fixated at least 90% of the trial (saves trials that were excluded for minor artifacts)
opts.monkey_name = monkey_name;
opts.outputdir = outputdir;

if abs(droptestcheck)>0.1;
    opts.spk_offset = droptestcheck;
else
    opts.spk_offset = 0;
end

[ExptTrials, ExptInfo] = ExtractAndAlignData( filenameP, dirpath, which_computer, ks, opts ); 
% to see options (including a different datadirectory, type 'help ExtractAndAlignData'

disp('Kofiko alignment complete')
