%% Kilo pipeline
clear all
%%
filenameP = '230331_144125_Jacomo'; iso_SUs=[];
%filenameP = '220921_140119_Jacomo';
%filenameP = '220323_123513_Jacomo';

datPath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_laminar/' filenameP '.dat'];
%%
cd('/media/felix/Internal_1/Data/BevilColor')
outputFolder = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_laminar/'];

% cd('/home/felix/NTlab_dataserver3/Conway/raw')
% outputFolder = ['/home/fellixbartsch/NTlab_dataserver2/Conway/raw/' filenameP '/kilosorting_laminar/'];

disp('Starting conversion')

opts.commonAverageReferencing = false; % TODO: set to true as experiment later
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;

opts.nChans = 24;
opts.ChnOffset=0;
opts.batch_size = 300000000;
%opts.specificChannels = [];  % user can select which plexon channels to use for
%                               conversion. remember, this must be in
%                               plexon-numbering, eg SPKC1 is usually ch num 65.
%
%%
[fs, n, ts, fn, ~] = plx_ad_v([filenameP '.pl2'], ['SPKC' num2str(1,'%03.f')] );
if n<2; error('error reading laminar probe'); end
%%
[samples, datPath] = convertRawToDat_Laminar([filenameP '.pl2'],opts);
%[samples, datPath] = convertRawToDatv3_Utah([filenameP '.pl2'],opts);

disp('Done with conversion')
%% Probe geometry
%[xcoords, ycoords, kcoords] = probeGeometry2coords('linear50', 24); % now
%lives in mastermegafile

%% edit the Megafile, save it into the target directory
%edit masterMegaFile_base_FB

%
% if ~exist(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/'],'dir')
%     mkdir(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/']);
% end
fs = 40000; %from above
nCh = opts.nChans; % = size(samples,1);

% set Parameters: 
% frequency for high pass filtering (150)
%ops.fshigh = 150; % KS2.5  
ops.fshigh = 300; % KS3: high-pass more aggresively

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4]; % KS2.5
%ops.Th = [9 9];  % KS3
%ops.Th = [9 8];

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 10;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.8; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 0; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask  = 30; 

%% KS3 only:
% spatial smoothness constant for registration
ops.sig        = 20;  % default = 20

% blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
ops.nblocks    = 10;  % default = 5

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; %defauilt=[8]

%%
step0b_masterMegaFile_base_FBv2(datPath, fs, nCh, ['linear50'], ops, 128)
%step0b_masterMegaFile_ks3_FB(datPath, fs, nCh, ['linear50'], ops, 128)
%%
disp('DONE with the first sorting pass WOOOOO')
%%
