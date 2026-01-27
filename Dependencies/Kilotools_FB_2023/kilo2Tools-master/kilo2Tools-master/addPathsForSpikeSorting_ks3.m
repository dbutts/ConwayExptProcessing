function [paths] = addPathsForSpikeSorting
%
% Different machines have different paths... This script adds the paths for
% your spike sorting needs given a particular machine.
% Fee free to add your own machine to the list...
%
% Necessary paths:
%   kiloTools
%   KiloSort-master
%   npy-matlab
%
% INPUT
%   none
% OUTPUT
%   paths - struct of folder locations

%% Define paths

% get hostname:
[~, hostName] = system('hostname');

% addpaths
if contains(hostName, 'DESKTOP-KEJGC64', 'IgnoreCase', 1)
    % Dell spike sorter in NIH krauzlis rig space 
    paths.lnkToolbox    = 'D:\Code\Toolboxes\0-lnkToolbox';
    paths.kilo2Tools    = 'D:\Code\Toolboxes\kilo2Tools';
    paths.kiloSort2     = 'D:\Code\Toolboxes\Kilosort2-master';
    paths.npymatlab     = 'D:\Code\Toolboxes\npy-matlab';
    
elseif contains(hostName, 'IT', 'IgnoreCase', 1) || ... % UMD
        contains(hostName, 'NEIK2A79LK07A', 'IgnoreCase', 1)
    
%    paths.lnkToolbox     = '~/Dropbox/Code/spike_sorting/0-lnkToolbox';
    paths.kilo2Tools     = '/home/felix/Dropbox/Project_BevilColor/kilo2Tools-master/kilo2Tools-master';
    paths.spikes         = '/home/felix/Dropbox/Project_BevilColor/kilo2Tools-master/spikes-master';
    paths.sortingQuality = '/home/felix/Dropbox/Project_BevilColor/kilo2Tools-master/sortingQuality-master';
    paths.kiloSort      = '/home/felix/Dropbox/LibFel/Git/Kilosort';
    paths.npymatlab      = '/home/felix/Dropbox/LibFel/Git/npy-matlab';
    paths.plexonSdk      = '/home/felix/Dropbox/LibFel/Bevil_RainnieLabInVivo/Plexon-Matlab Offline Files SDK';
elseif contains(hostName, 'mt', 'IgnoreCase', 1) || ... % UMD
        contains(hostName, 'NEIK2A79LK07A', 'IgnoreCase', 1)
    
%    paths.lnkToolbox     = '~/Dropbox/Code/spike_sorting/0-lnkToolbox';
    paths.kilo2Tools     = '/home/fellixbartsch/Dropbox/Project_BevilColor/kilo2Tools-master/kilo2Tools-master';
    paths.spikes         = '/home/fellixbartsch/Dropbox/Project_BevilColor/kilo2Tools-master/spikes-master';
    paths.sortingQuality = '/home/fellixbartsch/Dropbox/Project_BevilColor/kilo2Tools-master/sortingQuality-master';
    paths.kiloSort      = '/home/fellixbartsch/Dropbox/LibFel/Git/Kilosort';
    paths.npymatlab      = '/home/fellixbartsch/Dropbox/LibFel/Git/npy-matlab';
    paths.plexonSdk      = '/home/fellixbartsch/Dropbox/LibFel/Bevil_RainnieLabInVivo/Plexon-Matlab Offline Files SDK';
     
else
    error('Unrecognized hostname. Could not add necessary paths for sorting')
end

%% add paths
flds = fieldnames(paths);
for iF = 1:numel(flds)
    addpath(genpath(paths.(flds{iF})));
    disp(['added - ' flds{iF}]);
end

