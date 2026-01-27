function extractSpikesTDT(sessionName, varargin) 
%
% usage: extractSpikesTDT_20141205( sessionName, varargin )
%
% This function will                                                  % ADD IN THE FUNCTION DESCRIPTION HERE
%
% INPUTS:
%   sessionName - full path to the session to be analyzed
%
%   'tetrodelist' - list of tetrodes on which to perform the analysis
%       (default is to analyze all non-EEG/EMG channels). Enter as 't01',
%       'r03', etc. in a cell array
%   'threshold'
%   'numsigmasegments'                                                  % NEED TO ADD IN THE REST OF THE VARARGS

tetrodeList      = {};
rel_threshold    = 8;   % in units of standard deviation
numSigmaSegments = 60;  % number of segments to use to calculate the standard deviation of the signal on each wire
sigmaChunkLength = 1;   % duration in seconds of data chunks to use to extract the standard deviations of the wavelet-filtered signals
snle_window      = 48;    % Alex's default
r_upsample       = 2;     % the upsampling ratio
waveLength       = 24;    % width of waveform sampling window in A-D clock ticks
peakLoc          = 8;     % location of waveform peak in A-D clock ticks
deadTime         = 16;    % dead time in A-D clock ticks
overlapTolerance = 16;     % amount waveforms can overlap on different wires
% of the same tetrode (in A-D clock ticks) and
% still be counted as the same waveform
maxSNLE          = 10^5;  % maximum allowable non-linear energy. Any values
% greater than this are assumed to be noise and
% are not included in the standard deviation
% calculations (but are thresholded and
% extracted as potential spikes)
maxLevel = r_upsample + 4;    % max wavelet filtering level = upsampling ratio + 5
cutoff_Fs        = 9000; % cutoff of the anti-aliasing filter, needed for the sincInterp function
startTime        = 0;
endTime          = 0;
snle_T = 10;
override_nas_path = '';
override_conf_path = '';

for iarg = 1 : 2 : nargin - 1
    switch lower(varargin{iarg})
        case 'tetrodelist',
         tetrodeList = varargin{iarg + 1};
        case 'threshold',
         rel_threshold = varargin{iarg + 1};
        case 'numsigmasegments',
         numSigmaSegments = varargin{iarg + 1};
        case 'machineformat',
         machineFormat = varargin{iarg + 1};
        case 'upsampleratio',
         r_upsample = varargin{iarg + 1};
        case 'wavelength',
         waveLength = varargin{iarg + 1};
        case 'peakloc',
         peakLoc = varargin{iarg + 1};
        case 'deadtime',
         deadTime = varargin{iarg + 1};
        case 'overlaptolerance',
         overlapTolerance = varargin{iarg + 1};
        case 'maxsnle',
         maxSNLE = varargin{iarg + 1};
        case 'sincinterp_cutoff_fs',
        cutoff_Fs = varargin{iarg + 1};
        % added in to make it possible to set start and end times for 
        %processing (for example, in case there was a long period before 
        %recording was initiated after turning on the software).
        % this functionality has not been added to the code yet,
        % though.
        case 'starttime',
         startTime = varargin{iarg + 1};   
         % parallel processing doesn't like endTime being created in
         % two places
%         case 'endtime',
%          endTime = varargin{iarg + 1};
        case 'snlewindow',
         snle_window = varargin{iarg + 1};
        case 'snle_period',
         snle_T = varargin{iarg + 1};
        case 'override_nas_path'
         override_nas_path = varargin{iarg + 1};
        case 'override_conf_path'
         override_conf_path = varargin{iarg + 1};
    end
end

if isempty(override_conf_path) % use MySQL
    [~, ratID] = sql_getSubjectFromSession(sessionName);
    chMap = sql_getChannelMap(ratID);
    validMasks = sql_getAllTetChannels(sessionName);
    nasPath = sql_findNASpath(ratID);
else
    load(override_conf_path); % use session config
end
% will cause error if conf file is used and override_nas_path is not
% passed in, could handle more gracefully in future
if ~isempty(override_nas_path)
    nasPath = override_nas_path;
end

if isempty(tetrodeList); tetrodeList = chMap.tetNames; end
if ~iscell(tetrodeList); tetrodeList = {tetrodeList}; end
             
sessionTDTpath = fullfile(nasPath, ratID, [ratID '-rawdata'], sessionName, sessionName);
processedSessionPath = fullfile(nasPath, ratID, [ratID '-processed']);


if ~exist(processedSessionPath, 'dir')
    % defined above after nasPath extracted from the sql db
    mkdir(processedSessionPath);
end

tetChannels = zeros(length(tetrodeList), size(chMap.chMap, 2));
numValidTets = 0;
for iTet = 1 : length(tetrodeList)
    tetIdx = find(strcmpi(tetrodeList{iTet}, chMap.tetNames));
    if ~isempty(tetIdx)
        numValidTets = numValidTets + 1;
        tetChannels(numValidTets,:) = chMap.chMap(tetIdx,:);
    end
end

% turn on parallel worker pool
% clusterInfo = parcluster('local');
% poolobj = parpool(clusterInfo, clusterInfo.NumWorkers);
% disp(['Pool size: ',num2str(clusterInfo.NumWorkers)]);

tetChannels = tetChannels(1:numValidTets, :);
tetWireStd = zeros(numValidTets, 4);
for iTet = 1 : numValidTets %debug line
% parfor iTet = 1 : numValidTets
    % check for valid channels
    if ~any(validMasks(iTet,:))
        disp(['Skipping ' sessionName ', tetrode ' tetrodeList{iTet} ' - no valid channels']);
        continue;
    end
    disp(['calculating single wire standard deviations for tetrode ' tetrodeList{iTet}]);
    % it's not the standard deviation, but the median of the snle / 0.6745 
    tetWireStd(iTet, :) = extractSigma_snle_TDT(sessionTDTpath, ...            
                          tetChannels(iTet,2:end), ...
                          validMasks(iTet,:), ...
                          numSigmaSegments, sigmaChunkLength, r_upsample, ...
                          'snlewindow', snle_window, ...
                          'maxlevel', maxLevel, ...
                          'maxsnle', maxSNLE, ...
                          'sincinterp_cutoff_fs', cutoff_Fs, ...
                          'starttime', startTime, ...
                          'endtime', endTime, ...
                          'snle_period', snle_T);
end

tet_thresholds  = rel_threshold * tetWireStd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at this point, should have standard deviations for the wavelet filtered
% signal on each relevant wire - now time to do the thresholding!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTet = 1 : length(tetrodeList) %debug line
% parfor iTet = 1 : length(tetrodeList)
     % check for valid channels
    if ~any(validMasks(iTet,:))
        disp(['Skipping ' sessionName ', tetrode ' tetrodeList{iTet} ' - no valid channels']);
        continue;
    end
    
    extract_PLXtimestamps_sincInterp_TDT(sessionTDTpath, tetrodeList{iTet},...
        tetChannels(iTet,2:end), tet_thresholds(iTet, :), validMasks(iTet,:), ...
        'wavelength', waveLength, ...
        'peakloc', peakLoc, ...
        'snlewindow', snle_window, ...
        'deadtime', deadTime, ...
        'upsampleratio', r_upsample, ...
        'overlaptolerance', overlapTolerance, ...
        'starttime', startTime, ...
        'endtime', endTime, ...
        'snle_period', snle_T);
end

% delete(poolobj);