%% Kilo pipeline
clear all
%%
%cd('C:\Users\ConwayColor\Desktop\Kilotools')
%cd('/home/felix/Dropbox/Project_BevilColor')

filenameP = '221007_141759_Jacomo'; fs=40000; 
%filenameP = '220320_131300_Jacomo';

%datPath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah/' filenameP '.dat'];

%%
cd('/media/felix/Internal_1/Data/BevilColor')
outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_NForm/'];
opts.nChans = 28;
opts.specificChannels = [33,40,46,47,52,53,54,59,65,67,71,81,83,89,90,95,98,102,103,109,112,131,138,139,145,146,152,158];
opts.ChnOffset = 33; 

disp('Starting conversion')

opts.commonAverageReferencing = false; % TODO: set to true as experiment later
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;
opts.batch_size = 300000000;

%% 
datPath = [outputFolder filenameP '.dat'];
fs = 40000; 

%%
[fs, n, ts, fn, ~] = plx_ad_v([filenameP '.pl2'], ['SPKC' num2str(33,'%03.f')] );
if n<2; error('error reading N-Form'); end
%%
[samples, datPath] = convertRawToDatv3_NForm([filenameP '.pl2'], opts);
disp('Done with conversion')

%% edit the Megafile, save it into the target directory

nCh = opts.nChans; % = size(samples,1);

% extracted from MasterMegaFile code to set here
ops.chanMap = fullfile('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Nform_chanMap.mat');

% set Parameters: 
% frequency for high pass filtering (150)
%ops.fshigh = 150; % KS2.5  
ops.fshigh = 300; % KS3: high-pass more aggresively

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.001; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 5];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 2;  

% splitting a cluster at the end requires at least this much isolation for
% each sub-cluster (max = 1) [0.8]
ops.AUCsplit = 0.85; 

% minimum spike rate (Hz), if a cluster falls below this for too long it
% gets removed [1/50]
ops.minFR = 0; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

%% KS3 only:
% spatial smoothness constant for registration
ops.sig        = 20;  

% blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
ops.nblocks    = 5;  

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 7; 

%%
dbstop if error
Nfilt=256;
step0b_masterMegaFile_Plexon_FBv3(datPath, fs, nCh, ['linear150'], ops, Nfilt)
%step0b_masterMegaFile_ks3_Plexon_FB(datPath, fs, nCh, ['linear150'], ops, Nfilt)
%%
disp('DONE with the first sorting pass for NForm vWOOOOO')
%%
