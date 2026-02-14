%% Master file
clear; clc

% this is the name of the experiment you want to run
% filenameP = '251024_154348_Jacomo';
exptDate = '260204';
filenameP = '260204_161143_Sprout';
monkey_name = 'Sprout';

outputdir = '/home/bizon/Data/V1_Fovea/';
disp('Setup starting')

if strcmpi(getenv('USER'), 'bizon') % if this is true, we are on bizon
    % set path
    processingPath = '/home/bizon/Git/ConwayExptProcessing';
    % Switch into data directory
    dirpath = fullfile('/home/bizon/Data/V1_Fovea/',monkey_name, exptDate);
    pl2path = dirpath; %'/mnt/bc9/Data/'; % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
    stimpath = '/home/bizon/Processing/Cloudstims_calib_04_2024/';
else
    isilonPath = fullfile([filesep 'mnt'], 'isilon');
    if ~isdir(isilonPath)
        isilonPath = fullfile([filesep 'Volumes'], 'isilon');
    end

    processingPath = fullfile(isilonPath, 'code', 'ConwayExptProcessing');

    % Switch into data directory
    dirpath = uigetdr(pwd, 'Select directory with data');
    pl2path = fullfile(isilonPath, 'monkey_ephys', monkey_name); % you can load the plexon file directly from the server, which may make loading data slower but save you data transfer complications
    stimpath = fullfile(isilonPath, 'PROJECTS', 'V1_Fovea', 'stimuli', 'Cloudstims_calib_04_2024');
    
end

addpath(processingPath); % add necessary dependencies
addpath(genpath(fullfile(processingPath, 'Dependencies')));
addpath(fullfile(processingPath, 'Tools'));

plexon_fname = fullfile(pl2path, [filenameP '.pl2']);

%cd(dirpath)

disp('setup complete')

opts.monkey_name = monkey_name;
opts.batch_size = 5e9;
%% multiple arrays

arrayLabels = {'UT1', 'UT2', 'lam'};
nChans = [96, 96, 64];

if numel(nChans) > 1
    chnOffsets = cumsum([0 nChans(1:end-1)]);
else
    chnOffsets = 0;
end
curChannels = {{{1:32}, {33:64}, {65:96}}, {{1:32}, {33:64}, {65:96}}, {{1:64}}};
preconverted = {{zeros(1,3)}, {zeros(1,3)}, {0}};
arraySpacing = [1, 1, 0    % x
                1, 1, 50]; % y

disp('Kilosorting Starting')
for a = 1:numel(arrayLabels) % for each array
    opts.ArrayLabel = arrayLabels{a};

    try
        opts.chInfo = ...
            load(['/home/bizon/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/Sprout/V1_'...
            arrayLabels{a} '_chanMap_nogroup.mat']);
    catch
        if contains(opts.ArrayLabel, 'UT', 'IgnoreCase', true)
            error('No channel map info')
        end

        if isfield(opts, 'chInfo')
            opts = rmfield(opts, 'chInfo');
        end
    end

    for c = 1:numel(curChannels{a})
        opts.preconverted = preconverted{a}{:}(c);
        opts.nChans = length(curChannels{a}{c}{:});
        opts.ChnOffset = chnOffsets(a);
        opts.curchannels = curChannels{a}{c}{:};
        opts.arraySpacing = arraySpacing(:,a);

        [datPath{a,c}, droptestcheck{a,c}] = Step0_Kilosort(dirpath,filenameP,pl2path,opts); %this will take some time to run!

        disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck{a,c})]) % this will tell us if the plexon time alignment issue is present
        if abs(droptestcheck{a,c})>0.1
            warning("Danger - Plexon might have dropped frames! Check pl2 file.")
        else
            disp('Experiment kilosorted and ready for curation!')
        end
    end
end
disp('Kilosorting complete, go curate in phy')
disp('In terminal type:')
disp('conda activate phy2')
disp(['cd ' dirpath filenameP filesep 'kilosorting_laminar' filesep])
disp('phy template-gui params.py')


%% if you are rerunning the stuff below, uncomment this cell to quickly rerun the drop test check and see if plexon dropped any measurements
pl2 = PL2ReadFileIndex(plexon_fname);
numDigitsInLastSpkChan = ceil(log10(length(pl2.SpikeChannels)));

numAIchans = sum(cellfun(@(x) contains(x.Name, 'AI'), pl2.AnalogChannels));
numDigitsInLastAIchan = ceil(log10(numAIchans));
disp('Drop test check starting')

[fs, n, ts, fn, ~] = plx_ad_v(plexon_fname,...
    ['SPKC' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);

[fs_aux, n_aux, ts_aux, fn_aux, ~] = plx_ad_v(plexon_fname, ['AI' num2str(1, ['%0' num2str(numDigitsInLastAIchan), '.f'])] );

spikeChannel1SignalDurSec = n/fs; % samples / (samples/sec)
dpiSyncSignalDurSec = n_aux/fs_aux;

droptestcheck = n_aux/fs_aux - n/fs;
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present

if abs(droptestcheck)>0.1 
    warning("Danger - Plexon might have dropped frames! Check pl2 file.") 
else 
    disp('Experiment kilosorted and ready for curation!')
end

disp('Drop test check complete')

%% now, align kofiko information about stimuli with plexon data about spikes, LFPs, and eye traces
disp('Kofiko alignment starting')
which_computer = 2; % (default=2) 2 = Cameron's Bizon

ks.use_online = 0; % set to 1 to use on-line sorting or no sorting at all. Should be 0 if you want to use kilosort
%ks.onlinechans = [1:64]; % which channels of on-line sorted spikes should we go through? 

ks.stitched=0; % if you combined kilosort outputs for multiple arrays
%ks.arraylabel = opts.ArrayLabel;
ks.arraylabel = arrayLabels;

ii = 1;
for a = 1:numel(arrayLabels)
    for c = 1:numel(curChannels{a})
        if contains(arrayLabels{a}, 'lam')
            ksFolders{ii} = fullfile(dirpath, filenameP, 'kilosorting_laminar');
            ks_nChans(ii) = length(curChannels{a}{c}{:});
        else
            ksFolders{ii} = ['kilosorting_' arrayLabels{a} '_' num2str(min(curChannels{a}{c}{:})) 'to' num2str(max(curChannels{a}{c}{:}))];
            ksFolders{ii} = fullfile(dirpath, filenameP, ksFolders{ii});
            ks_nChans(ii) = length(curChannels{a}{c}{:});
            ii = ii + 1;
        end
    end
end

ks.filepath = ksFolders; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
%ks.filepath = [dirpath filenameP filesep 'kilosorting_stitched' filesep]; % point this at array folders or the "stiched" folder if you want to sort data from multiple arrays
ks.pl2path = pl2path;
ks.ks_nChans = ks_nChans;
opts.eye_tracker = 4;   % (default=3) 0=eyescan, 1=monoc eyelink, 2=binoc eyelink, 3=monocular dDPI, 4=binocular dDPI
opts.is_cloud = 1;      % (default=1) indicates processing for cloud data. set to 0 to skip cloud-specific variables and align task data or other paradigms 
opts.extractfixinfo = 1;% (default=1) extracts fixinfo file containing fivedot and dotgrid info
opts.trialwindow = [0, 4]; % 4 second trials
opts.trl_fix_thresh = 0.6; % include trials with at least 3 seconds of fixation
opts.monkey_name = monkey_name;

if abs(droptestcheck)>0.1
    opts.spk_offset = droptestcheck;
else
    opts.spk_offset = 0;
end

[ExptTrials, ExptInfo, ETdata, LFP_ad] = ExtractAndAlignData(filenameP, dirpath, which_computer, ks, opts); 
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
addpath(genpath('/home/bizon/Data/V1_Fovea/output_greenemj')); 

data = PackageCloudData_v9mod( ExptTrials, ExptInfo, [], [], stimpath, fullfile(pl2path, filenameP, 'Analysis') , 0, which_computer); % temporarily got rid of LFP_ad as final arg

disp('Cloud packaging complete')
%% Check to see if all of the relevant files were saved
checkOutput(outputdir,filenameP)

%% finally, generate STAs for all recorded cells
disp('STA generation starting')
%%force save STA plots [RamonBartolo 20250709]
data.to_save = true;
data.outputdir = fullfile(pl2path, filenameP, 'Analysis', 'STAs'); %[outputdir filenameP filesep 'Analysis' filesep, 'STAs'];

%%%%%% Comment to skip saving
if ~exist(data.outputdir,'dir'); mkdir(data.outputdir); end
stas = get_sta(data, ExptInfo);
disp('STA generation complete')

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
