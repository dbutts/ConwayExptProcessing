%% Kilo pipeline
clear all
%%
%cd('C:\Users\ConwayColor\Desktop\Kilotools')
%cd('/home/felix/Dropbox/Project_BevilColor')

%filenameP = '220314_143059_Jacomo'; fs=40000; datPath = '/media/felix/Internal_1/Data/BevilColor/220314_143059_Jacomo/kilosorting_Utah/220314_143059_Jacomo.dat';
%filenameP = '221003_143523_Jacomo';
filenameP = '221007_141759_Jacomo'; %rerun utah

%datPath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah/' filenameP '.dat'];

%%
for current_run = [1:4]

%%
cd('/media/felix/Internal_1/Data/BevilColor')

switch current_run
    case 0
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah'];
        opts.nChans = 96;
        opts.ChnOffset=160;        
    case 1
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_1to24/'];
        opts.nChans = 24;
        opts.ChnOffset=160;
    case 2
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_25to48/'];
        opts.nChans = 24;
        opts.ChnOffset=184;
    case 3
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_49to72/'];
        opts.nChans = 24;
        opts.ChnOffset=208;
    case 4
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_73to96/'];
        opts.nChans = 24;
        opts.ChnOffset=232;
    case 5
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_1to48/'];
        opts.nChans = 48;
        opts.ChnOffset=160;
    case 6
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_49to96/'];
        opts.nChans = 48;
        opts.ChnOffset=208;
    case 11
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_1to16/'];
        opts.nChans = 16;
        opts.ChnOffset=160;
    case 12
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_17to32/'];
        opts.nChans = 16;
        opts.ChnOffset=176;
    case 13
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_33to48/'];
        opts.nChans = 16;
        opts.ChnOffset=192;
    case 14
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_49to64/'];
        opts.nChans = 16;
        opts.ChnOffset=208;
    case 15
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_65to80/'];
        opts.nChans = 16;
        opts.ChnOffset=224;
    case 16
        outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_81to96/'];
        opts.nChans = 16;
        opts.ChnOffset=240;
end       
disp('Starting conversion')

opts.commonAverageReferencing = false; % TODO: set to true as experiment later
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;

opts.batch_size = 300000000;

%opts.specificChannels = [];  % user can select which plexon channels to use for
%                               conversion. remember, this must be in
%                               plexon-numbering, eg SPKC1 is usually ch num 65.
%
%% 
datPath = [outputFolder filenameP '.dat'];
fs = 40000; 

%%
[fs, n, ts, fn, ~] = plx_ad_v([filenameP '.pl2'], ['SPKC' num2str(161,'%03.f')] );
if n<2; error('error reading Utah data'); end
%%
[samples, datPath] = convertRawToDatv3_Utah([filenameP '.pl2'], opts);
disp('Done with conversion')

%% Probe geometry
%[xcoords, ycoords, kcoords] = probeGeometry2coords('linear50', 24); % now
%lives in mastermegafile

%% edit the Megafile, save it into the target directory
%edit masterMegaFile_base_FB

%datPath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_laminar/' filenameP '.dat'];

% if ~exist(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/'],'dir')
%     mkdir(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/']);
% end

%fs = 40000; %from above
nCh = opts.nChans; % = size(samples,1);

% extracted from MasterMegaFile code to set here
switch current_run
    case 0; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_nogroup.mat');
    case 1; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_1to24nogroup.mat');
    case 2; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_25to48nogroup.mat');
    case 3; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_49to72nogroup.mat');
    case 4; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_73to96nogroup.mat');
    case 5; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_1to48nogroup.mat');
    case 6; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_49to96nogroup.mat');
    case 11; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_1to16nogroup.mat');
    case 12; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_17to32nogroup.mat');
    case 13; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_33to48nogroup.mat');
    case 14; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_49to64nogroup.mat');
    case 15; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_65to80nogroup.mat');
    case 16; ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap_81to96nogroup.mat');
end

% set Parameters: 
% frequency for high pass filtering (150)
ops.fshigh = 150; % KS2.5  

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.01; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 5];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 3;  

% splitting a cluster at the end requires at least this much isolation for
% each sub-cluster (max = 1) [0.8]
ops.AUCsplit = 0.95; 

% minimum spike rate (Hz), if a cluster falls below this for too long it
% gets removed [1/50]
ops.minFR = 0; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; %[8]
%% KS3 only:
%ops.fshigh = 300; % KS3: high-pass more aggresively

% spatial smoothness constant for registration
ops.sig        = 20;  % default = 20

% blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
ops.nblocks    = 0;  % default = 5

%%
dbstop if error
Nfilt=256;
step0b_masterMegaFile_Utah_FBv2(datPath, fs, nCh, ['linear50'], ops, Nfilt)
%step0b_masterMegaFile_ks3_Utah_FB(datPath, fs, nCh, ['linear50'], ops, Nfilt)
%%
disp('DONE with the first sorting pass WOOOOO')
%%
end