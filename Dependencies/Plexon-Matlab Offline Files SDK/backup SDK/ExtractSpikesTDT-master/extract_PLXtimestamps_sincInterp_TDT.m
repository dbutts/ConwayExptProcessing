function extract_PLXtimestamps_sincInterp_TDT( rawdataPath, tetName, chList, thresholds, validMask, varargin )
%
% usage: extract_PLXtimestamps_sincInterp_TDT_20141208( hsdFile, chList, thresholds, varargin )
%
% INPUTS:
%   rawdataPath - name of the path to the raw data for the current session
%   chList - list of channels for this tetrode
%   thresholds - vector containing the thresholds for each wire in chList
%   validMask - boolean vector the same length as chList. "0"s indicate bad
%       channels that should be ignored; "1"s indicate good channels that
%       should be included in the calculations
%
% VARARGs:
%   maxlevel - maximum level for the wavelet filter (default 5 + the
%       upsampling ratio)
%   wavelength - duration of waveforms in samples
%   peakloc - location of peaks within waveforms in samples
%   deadtime - dead time required before detecting another spike, in
%       samples
%   sinclength - length of the sinc function to use for upsampling
%   upsampleratio - ratio by which to upsample the signal (ie, a value of 2
%       would take Fs from, for example, 30 kHz to 60 kHz)
count=0;
deadTime   = 16;     % dead time after a spike within which another spike
                     % cannot be detected on the same wire (in samples)
peakLoc    = 8;      % number of samples to look backwards from a peak (ie,
                     % peaks should be aligned peakLoc samples into each
                     % waveform
waveLength = 24;     % duration of each waveform in number of samples
sincLength = 13;     % length of sinc function for sinc interpolation
r_upsample = 2;      % upsampling ratio - ie, if original Fs = 31250 and
                     % r_upsample = 2, the sampling rate for each waveform
                     % will be 62500 Hz.
maxLevel   = 0;      % max level for wavelet filtering
snle_T     = 5;
overlapTolerance = 1;
cutoff_Fs = 9000;
windowSize = 12 * r_upsample;

% note that deadTime, peakLoc, waveLength, and sincLength are in units of
% number of samples for the ORIGINAL signal. That is, if the signal is
% upsampled by a factor of 2, the deadTime, etc. written to the .plx file
% will be 2 * the deadTime supplied above (or as a varargin).

for iarg = 1 : 2 : nargin - 5
    switch lower(varargin{iarg})
        case 'maxlevel',
            maxLevel = varargin{iarg + 1};
        case 'wavelength',
            waveLength = varargin{iarg + 1};
        case 'peakloc',
            peakLoc = varargin{iarg + 1};
        case 'deadtime',
            deadTime = varargin{iarg + 1};
        case 'sinclength',
            sincLength = varargin{iarg + 1};
        case 'upsampleratio',
            r_upsample = varargin{iarg + 1};
        case 'overlaptolerance',
            overlapTolerance = varargin{iarg + 1};
        case 'sincinterp_cutoff_fs',
            cutoff_Fs = varargin{iarg + 1};
        case 'snlewindow',
            windowSize = varargin{iarg + 1};
        case 'starttime',
            startTime = varargin{iarg + 1};            % added in to make it possible to set start and end times for processing (for example, in case there was a long period before recording was initiated after turning on the software).
        case 'endtime',                                % this functionality has not been added to the code yet, though. -DL 12/15/2014
            endTime = varargin{iarg + 1};
        case 'snle_period',
            snle_T = varargin{iarg + 1};
    end
end

r_upsample = round(r_upsample);

if maxLevel == 0
    maxLevel = r_upsample + 4;      % cutoff frequency = samplingrate/(2^(maxlevel+1))
                                    % this should make the cutoff frequency
                                    % ~230 Hz for an initial sampling rate
                                    % of ~30 kHz. For an initial sampling
                                    % rate of 20 kHz, the cutoff will be
                                    % ~150 Hz. May want to use r_upsample +
                                    % 4 if Fs = 20 kHz (cutoff ~300 Hz)
end

% load the header from the first .sev file for this tetrode to get
% recording parameters
cd(rawdataPath);
[rawDataSessionPath, sessionName, ~] = fileparts(rawdataPath);
[rawDataTopPath, ~, ~] = fileparts(rawDataSessionPath);
[ratTopPath, ~, ~] = fileparts(rawDataTopPath);
ratID = sessionName(1:5);
processedDataPath = fullfile(ratTopPath, [ratID '-processed'], sessionName);
if ~exist(processedDataPath, 'dir')
    mkdir(processedDataPath);
end

sessionDateStr = sessionName(7:14);
sessionDateVec = datevec(sessionDateStr, 'yyyymmdd');

% find the names of the .sev files
sevNames = cell(length(chList), 1);
sevInfo = dir('*.sev');
if isempty(sevInfo)
    error('extractSigma_snle_TDT_20141207:noSevFiles',['Cannot find sev files for ' sessionName]);
end
firstValidSev = 0;
for iCh = 1 : length(chList)
    for iSev = 1 : length(sevInfo)
        header = getSEVHeader(sevInfo(iSev).name);
        if header.channelNum == chList(iCh)
            sevNames{iCh} = sevInfo(iSev).name;
            if firstValidSev == 0
                firstValidSev = iSev;
            end
        end
    end
end
if ~firstValidSev
    error('extract_PLXtimestamps_sincInterp_TDT_20141208:noSevFilesForRequestedChannels',['Cannot find sev files for current tetrode in ' sessionName]);
end

Fs         = header.Fs;
numCh      = length(chList);
datalength = (header.fileSizeBytes - header.dataStartByte) / header.sampleWidthBytes;

if Fs == 0
    Fs = 24414.0625;
end

blockSize   = round(Fs * 10);    % process 10 sec at a time
overlapSize = round(Fs * 0.1);   % 100 ms overlap between adjacent blocks 
                                 % to avoid edge effects
final_Fs         = r_upsample * Fs;
final_peakLoc    = r_upsample * peakLoc;
final_waveLength = r_upsample * waveLength;

% make sure chList and thresholds are column vectors
if size(chList, 1) < size(chList, 2); chList = chList'; end
if size(thresholds, 1) < size(thresholds, 2); thresholds = thresholds'; end

% is it more efficient to read single wires in sequence, or read in a big
% chunk of data including all wires, then pull out the one to four wires of
% interest? I think the latter... - DL 3/27/2012

numBlocks = ceil(datalength / blockSize);

% numBlocks = 3;                          % just for debugging
% datalength = round(blockSize * (numBlocks - 0.5));    % just for debugging

% write the .plx header
spikeParameterString = sprintf('WL%02d_PL%02d_DT%02d', waveLength, peakLoc, deadTime);
PLX_fn = fullfile(processedDataPath, [sessionName '_' tetName '_' spikeParameterString '.plx']);

plxInfo.comment    = '';
plxInfo.ADFs       = final_Fs;           % record the upsampled Fs as the AD freq for timestamps
plxInfo.numWires   = length(chList);
plxInfo.numEvents  = 0;
plxInfo.numSlows   = 0;
plxInfo.waveLength = final_waveLength;
plxInfo.peakLoc    = final_peakLoc;

plxInfo.year  = sessionDateVec(1);
plxInfo.month = sessionDateVec(2);
plxInfo.day   = sessionDateVec(3);

timeVector = datevec('12:00', 'HH:MM');
plxInfo.hour       = timeVector(4);
plxInfo.minute     = timeVector(5);
plxInfo.second     = 0;
plxInfo.waveFs     = final_Fs;          % record the upsampled Fs as the waveform sampling frequency
plxInfo.dataLength = datalength * r_upsample;

plxInfo.Trodalness     = length(chList); %Trodalness - 0,1 = single electrode, 2 = stereotrode, 4 = tetrode
plxInfo.dataTrodalness = 0; %this is set to 0 in the Plexon tetrode sample file

plxInfo.bitsPerSpikeSample = 16;
plxInfo.bitsPerSlowSample  = 16;

plxInfo.SpikeMaxMagnitudeMV = 10000;    % +/- 1 V dynamic range on DAQ cards
plxInfo.SlowMaxMagnitudeMV  = 10000;    % +/- 1 V dynamic range on DAQ cards (probably not relevant for Berke lab systems)

%THIS SHOWS AS 65536 in OFFLINE SORTER, whats's wrong?
plxInfo.SpikePreAmpGain     = 1;        % gain before final amplification stage

PLXid = fopen(PLX_fn, 'w');
disp('PLX file opened...')
writePLXheader(PLXid, plxInfo );

% subjectName = strrep(header.name(1:end-9), '-', '');   % get rid of any hyphens in the subject name
% dateString  = sprintf('%04d%02d%02d', plxInfo.year, plxInfo.month, plxInfo.day);
% baseName    = [subjectName dateString];

for iCh=1:length(chList)
   
    chInfo.tetName  = [sessionName '_' tetName];
    chInfo.wireName = sprintf('%s_W%02d', chInfo.tetName, iCh);
    
    chInfo.wireNum   = iCh; %tetrode number
    chInfo.WFRate    = final_Fs;
    chInfo.SIG       = chList(iCh);  %channel number
    chInfo.refWire   = 0;     % not sure exactly what this is; Alex had it set to zero
    chInfo.gain      = 300;
    chInfo.filter    = 0;    % not sure what this is; Alex had it set to zero
    chInfo.thresh    = int32(thresholds(iCh));
    chInfo.numUnits  = 0;    % no sorted units
    chInfo.sortWidth = final_waveLength;
    chInfo.comment   = 'created by extractSpikesTDT.m';

    writePLXChanHeader(PLXid, chInfo );
end

% startSample = 1;
% stopSample = round(Fs*10);
% qtname = [sessionName '_' tetName '.mov'];

%             MakeQTMovie('start',qtname)
%             MakeQTMovie('quality', 0.1);

sevs = [];
for iCh = 1 : numCh
    if validMask(iCh)
        disp(['Building sev for ',sevNames{iCh}]);
        [sevs(iCh,:), ~] = read_tdt_sev(sevNames{iCh});
    end
end
    

PLXdata = {};
for iBlock = 1:numBlocks
    disp(sprintf('Extracting waveforms for %s, block %d of %d', ...
                 tetName, ...
                 iBlock, ...
                 numBlocks))
    %disp(['Finding timestamps and extracting waveforms for block ' num2str(iBlock) ' of ' num2str(numBlocks)]);
    
    rawData_curSamp   = (iBlock - 1) * blockSize;
    upsampled_curSamp = rawData_curSamp * r_upsample;
    
    % get overlapSize samples on either side of each block to prevent edge
    % effects (may not be that important, but it's easy to do)
    startSample = max(1, rawData_curSamp - overlapSize);
    if iBlock == 1
        numSamples  = blockSize + overlapSize;
    elseif iBlock == numBlocks
        numSamples = datalength - startSample + 1;
    else
        numSamples = blockSize + 2 * overlapSize;
    end
        
    stopSample = startSample + numSamples-1;
    
    rawData = zeros(numCh, (stopSample - startSample + 1));
    for iCh = 1 : numCh
        if validMask(iCh)
            rawData(iCh, :) = sevs(iCh,startSample : stopSample);
        end
    end

    if r_upsample > 1    
        interp_rawData = zeros(size(rawData, 1), size(rawData, 2) * r_upsample);
        for iCh = 1 : size(rawData, 1)
            if validMask(iCh)
                interp_rawData(iCh, :) = sincInterp(rawData(iCh, :), Fs, ...
                    cutoff_Fs, final_Fs, 'sinclength', sincLength);
            end
        end
        fdata = wavefilter(interp_rawData, maxLevel, 'validmask', validMask);
    else
        % wavelet filter the raw data
        % Don't bother to do the calculations for noisy wires.
        fdata = wavefilter(rawData, maxLevel, 'validmask', validMask);
    end
    
    SNLEdata = snle(fdata, validMask, 'windowsize', windowSize, ...
                     'snle_period', snle_T);
    
    % extract the timestamps of peaks in the smoothed non-linear energy
    % signal that are above threshold. Exclude wires with noisy recordings
    % from timestamp extraction.
    ts = gettimestampsSNLE(SNLEdata, thresholds, validMask, ...
                           'deadtime', deadTime * r_upsample, ...
                           'overlaptolerance', overlapTolerance * r_upsample);
                       
    % make sure peaks above threshold are not contained in the overlap
    % regions for adjacent blocks of data (and also that the first peak
    % location has enough data before it to extract a full waveform, and
    % the last spike has enough data after it to extract the full
    % waveform).
    switch iBlock
        case 1,
            block_ts = ts(ts > final_peakLoc & ts <= blockSize * r_upsample);
            ts = block_ts;
        case numBlocks,
            block_ts = ts((ts >= overlapSize * r_upsample + 2) & ...
                      (ts < (size(SNLEdata,2) - (final_waveLength - final_peakLoc))));
            ts = block_ts - (overlapSize * r_upsample + 1);
        otherwise,
            block_ts = ts((ts >= overlapSize * r_upsample + 2) & ...
                 (ts <= overlapSize * r_upsample + 1 + blockSize * r_upsample));
            ts = block_ts - (overlapSize * r_upsample + 1);
    end
    % NOTE: ts is timestamps in samples, not in real time. Divide by the
    % sampling rate to get real time
    
    if isempty(ts); continue; end

    waveforms = extractWaveforms(fdata, block_ts, final_peakLoc, final_waveLength);
    ts = ts + upsampled_curSamp;
    writePLXdatablock(PLXid, waveforms, ts);
end

fclose(PLXid);
disp('PLX file closed...');

end