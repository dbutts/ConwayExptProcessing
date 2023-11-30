%% Master file

addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/')
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/iCSD/CSD_functions/')
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilosort2/')
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/')
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Plexon-Matlab Offline Files SDK/')
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/Dependencies/Kilotools_FB_2023/kilo2Tools-master/npy-matlab/npy-matlab')
% my directory
addpath('/Users/dbutts/Projects/ColorV1/ColorProcessing_Packaged/DAB')

% Switch into data directory
cd /Users/dbutts/Data/Conway/
%dirpath = '/Volumes/DanData/ConwayData/';
dirpath = '/Users/dbutts/Data/Conway/';

% After running step0 and kilosorting, now replace step 1 with the following command (in matlab):
% Running step0 will have generated 230510_141725_Jacomo
[ExptTrials, ExptInfo] = ExtractAndAlignData( '230510_141725_Jacomo', dirpath ); 
% to see options (including a different datadirectory, type 'help ExtractAndAlignData'
% this will automatically generate and save CSD/LFP, but you can call the
% function again and analyze them

% ExtractAndAlign saves FullExpt_ks1_lam_v08.mat in Analysis directory, and
% can be loaded or its outputs (in memory can be used directly


% if step 1 done
load Analysis/230510_141725_Jacomo_FullExpt_ks1_lam_v09.mat

data = PackageCloudData_v9( ExptTrials, ExptInfo, [], [], stimpath, dirpath );
% data = PackageCloudData_v9( ExptTrials, ExptInfo );
