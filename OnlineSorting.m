%% Set paths

% Add ConwayExptProcessing Dependencies to path
% whats needed from this are functions to read npy and tsv files

clc;

username = char(java.lang.System.getProperty('user.name'));
homeDir = char(java.lang.System.getProperty('user.home'));

switch username
    case 'bizon'
        codedir = fullfile(homeDir, 'Git');
        stimpath = fullfile(homeDir, 'Processing');
        computerLocation = 1;
        disp('Location: Bizon')
    case 'clab'
        codedir = fullfile(homeDir, 'Git');
        stimpath = fullfile(homeDir, 'Processing');
        computerLocation = 1;
        disp('Location: clab')

    case 'dab'
        [status, hostname] = system('hostname');   % need to know which computer
        codedir = '/home/dab/Code/';
        if hostname(2) == '1'
            disp('Location: M1')
            computerLocation = 10;
            dirpath =  '/Data/ColorV1/raw/';
            stimpath = '/media/dab/Internal_1/Data/BevilColor/Cloudstims_calib_01_2022/';
            savepath = '/Data/ColorV1/Processed/';
        else
            disp('Location: MST')
            computerLocation = 11;
            dirpath =  '/mnt/NTdata3/Conway/raw/';
            stimpath = '/mnt/NTdata1/Conway/Cloudstims_calib_01_2022/';
            savepath = '/home/dab/ColorV1/Processed/';
        end
    case 'greenemj' % Max Greene's laptop
        codedir = fullfile(homeDir, 'Git');
        stimpath = '/Volumes/lsr-conway/PROJECTS/V1_Fovea/stimuli/';
        computerLocation = 1;
        disp('Location: Max Greene''s laptop');
end

%% flags
fprintf('Flags \n')
fprintf('--------------\n')


% This only needs to run if there is one -- so should be detected
% you don't have any except old experiments, so I'll just set this based on
% config

% print defaults
% I don't need all these prompts so will just move the ones I do out and
% all the rest can be for how you want to enter information

fprintf('Set paths \n')
fprintf('--------------\n')

manualFnames = input('Input full file paths manually? \n [0] No \n [1] Yes \n');

 stimET = input(['stimET? \n' ...
        '[1] Yes \n' ...
        '[0] No \n ']); %1;

if manualFnames
    plexon_fname = input('Plexon file: ', 's');
    mainKofiko_fname = input('Main kofiko file: ', 's');
    kofiko_subfolder = input('Kofiko subfolder: ', 's');
else

    try
        kofiko_subfolder = input('Kofiko subfolder (no slash at the end): ', 's');
        mainKofiko_fname = [kofiko_subfolder '.mat'];
        plexon_fname = [kofiko_subfolder '.pl2'];
    catch
        fprintf('Files not organized as expected, specify full paths \n')
        plexon_fname = input('Plexon file: ', 's');
        mainKofiko_fname = input('Main kofiko file: ', 's');
        kofiko_subfolder = input('Kofiko subfolder: ', 's');
    end

end

stimfolder = input(['stimfolder \n' ...
    '[1] Cloudstims_calib_01_2022 \n'...
    '[2] Cloudstims_calib_04_2024 \n']);

if stimfolder == 1
    stimfolder = 'Cloudstims_calib_01_2022';
elseif stimfolder == 2
    stimfolder = 'Cloudstims_calib_04_2024';
end
stimpath = fullfile(stimpath, stimfolder);



% save time with these assertions
%assert(isdir(dirpath) & isdir(ks_path) & isdir(stimpath), 'Check that dirpath, ks_path, and stimpath exist!');

[~, filenameP] = fileparts(mainKofiko_fname);
exptdate  = str2num(filenameP(1:6));


%% Initialization for everyone given setup
addpath(genpath(fullfile(codedir, 'ConwayExptProcessing/Dependencies/')));

pl2 = PL2ReadFileIndex(plexon_fname);

% Hardcoded values
plexonAnalogScale = 1e-3;
LumScale = 0.1085;
minFixationDuration = 0.6;
maxFixationErrorPix = 45;
minSpikes = 2000;
targ_ETstimtype = 0;
ETdist_thresh=40;
nLags = 10;
%% Load kofiko data
tic;
fprintf('Loading Kofiko trial data\n')

% I think this code works with Dan lab -- dont need if here
%if computerLocation < 10
kofiko_dir = dir(fullfile(kofiko_subfolder, [filenameP '*.mat']));
kofiko_fname_re = [filenameP '(_\d*)?.mat'];
valid_kofiko_fname_idx = cellfun(@(x) ~isempty(x), ...
    cellfun(@(x) regexp(x, kofiko_fname_re), {kofiko_dir(:).name}, 'UniformOutput',false));
kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
kofiko_fname = {kofiko_dir(valid_kofiko_fname_idx).name};

% Get the Kofiko mat file numbers from their names
kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'),...
    kofiko_fname, 'UniformOutput',false);

kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);


% Load the g_strctLocalExperimentRecording structs which contain
% information about each trials
g_strctLocalExperimentRecording = [];
for i = 1:length(kofiko_file_num)
    idx = kofiko_file_num(i);
    temp_full_fname = fullfile(kofiko_folder{i}, kofiko_fname{i});
    temp_strct = load(temp_full_fname, 'g_strctLocalExperimentRecording');
    temp_cell = temp_strct.g_strctLocalExperimentRecording;
    g_strctLocalExperimentRecording{idx} = temp_cell;
end

% add main kofiko file

temp_strct = load(mainKofiko_fname, 'g_strctLocalExperimentRecording');
temp_cell = temp_strct.g_strctLocalExperimentRecording;
g_strctLocalExperimentRecording{end+1} = temp_cell;

% Load the file which should be the "main" Kofiko mat file (no number at the end)

%mainKofikoFolderIdx = strcmpi(kofiko_fname, [filenameP '.mat']);
%load(fullfile(kofiko_folder{mainKofikoFolderIdx}, kofiko_fname{mainKofikoFolderIdx}), 'g_strctDAQParams', ...
load(mainKofiko_fname, 'g_strctDAQParams', ...
    'g_astrctAllParadigms', 'g_strctLog', 'g_strctEyeCalib', ...
    'g_strctAppConfig', 'g_strctSharedParadigmData', ...
    'g_strctStimulusServer', 'g_strctSystemCodes');

% Concatenate the g_strctLocalExperimentRecorindg structs and remove empty cells
trialData = vertcat(g_strctLocalExperimentRecording{:});
trialData(cellfun(@isempty, trialData)) = [];

% Unify field names across g_strctLocalExperimentRecording to facilitate using cellfun
fieldNames = cellfun(@(x) fieldnames(x), trialData, 'UniformOutput', false);
uniqueFieldNames = unique(vertcat(fieldNames{:}));
missingFieldNames = cellfun(@(x) setdiff(uniqueFieldNames, fieldnames(x)),...
    trialData, 'UniformOutput',false);

for i = 1:numel(trialData)
    for  j = 1:numel(missingFieldNames{i})
        trialData{i}.(missingFieldNames{i}{j}) = [];
    end
end
toc;

%% %%%%%%%% Time alignment %%%%%%%%%%
tic;
fprintf('Loading Kofiko eye signals and converting to plexon time\n')

% Read kofiko strobes from plexon and kofiko files
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);

kofikoSyncStrobesTS = transpose(g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == g_strctSystemCodes.m_iSync));
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iSync);

% juice times in plexon time
REW_ON_ts = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iJuiceON);
REW_OFF_ts = events.timeStamps(events.strobeNumber == g_strctSystemCodes.m_iJuiceOFF);

% Linear regression to get kofiko time stamps into plexon time
B = [ones(size(kofikoSyncStrobesTS)) kofikoSyncStrobesTS]\plexonSyncStrobesTS;

sync_ch = 'AI01';
[adfreq, n, ts, fn, sync_ad] = plx_ad_v(plexon_fname, sync_ch);
t_plexon = (0:n-1)/adfreq;

%% Kofiko eye data
% Get the screen dimensions
ScreenSizeX_pix = g_strctStimulusServer.m_aiScreenSize(3);
ScreenSizeY_pix = g_strctStimulusServer.m_aiScreenSize(4);

% Convert Kofiko eye signal timestamps to Plexon time
Kofiko_ET_TS = g_strctEyeCalib.EyeRaw.TimeStamp';
Kofiko_ET_TS_PlexonTime = [ones(size(Kofiko_ET_TS)) Kofiko_ET_TS]*B;

% Get Kofiko raw eye signals
Kofiko_EyeRawX = g_strctEyeCalib.EyeRaw.Buffer(:,1);
Kofiko_EyeRawY = g_strctEyeCalib.EyeRaw.Buffer(:,2);
% Get Kofiko eye signal gains
Kofiko_GainX = g_strctEyeCalib.GainX.Buffer;
Kofiko_GainY = g_strctEyeCalib.GainY.Buffer;
% Get Kofiko eye signal gain timestamps & convert to plexon time
Kofiko_GainX_TS = g_strctEyeCalib.GainX.TimeStamp;
Kofiko_GainX_TS_PlexonTime = [ones(size(Kofiko_GainX_TS')) Kofiko_GainX_TS']*B;

Kofiko_GainY_TS = g_strctEyeCalib.GainY.TimeStamp;
Kofiko_GainY_TS_PlexonTime = [ones(size(Kofiko_GainY_TS')) Kofiko_GainY_TS']*B;

% Get Kofiko eye signal offsets
Kofiko_CenterX = g_strctEyeCalib.CenterX.Buffer;
Kofiko_CenterY = g_strctEyeCalib.CenterY.Buffer;
% Get Kofiko eye signal offset timestamps
Kofiko_CenterX_TS = g_strctEyeCalib.CenterX.TimeStamp;
Kofiko_CenterX_TS_PlexonTime = [ones(size(Kofiko_CenterX_TS')) Kofiko_CenterX_TS']*B;

Kofiko_CenterY_TS = g_strctEyeCalib.CenterY.TimeStamp;
Kofiko_CenterY_TS_PlexonTime = [ones(size(Kofiko_CenterY_TS')) Kofiko_CenterY_TS']*B;

% Determine the Kofiko gains and offsets for each eye signal time stamp

% Kofiko X gain
Kofiko_GainX_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_GainX_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_GainX)
    kidx = Kofiko_ET_TS >= Kofiko_GainX_TS(i);
    pidx = t_plexon >= Kofiko_GainX_TS_PlexonTime(i);

    Kofiko_GainX_forEachTimeStamp(kidx) = Kofiko_GainX(i);
    Kofiko_GainX_forEachPlexonSample(pidx) = Kofiko_GainX(i);
end

% Kofiko Y gain
Kofiko_GainY_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_GainY_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_GainY)
    kidx = Kofiko_ET_TS >= Kofiko_GainY_TS(i);
    pidx = t_plexon >= Kofiko_GainY_TS_PlexonTime(i);

    Kofiko_GainY_forEachTimeStamp(kidx) = Kofiko_GainY(i);
    Kofiko_GainY_forEachPlexonSample(pidx) = Kofiko_GainY(i);
end

% Kofiko X center
Kofiko_CenterX_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_CenterX_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_CenterX)
    kidx = Kofiko_ET_TS >= Kofiko_CenterX_TS(i);
    pidx = t_plexon >= Kofiko_CenterX_TS_PlexonTime(i);

    Kofiko_CenterX_forEachTimeStamp(kidx) = Kofiko_CenterX(i);
    Kofiko_CenterX_forEachPlexonSample(pidx) = Kofiko_CenterX(i);
end

% Kofiko Y center
Kofiko_CenterY_forEachTimeStamp = nan(size(Kofiko_ET_TS));
Kofiko_CenterY_forEachPlexonSample = nan(size(t_plexon));

for i = 1:numel(Kofiko_CenterY)
    kidx = Kofiko_ET_TS >= Kofiko_CenterY_TS(i);
    pidx = t_plexon >= Kofiko_CenterY_TS_PlexonTime(i);

    Kofiko_CenterY_forEachTimeStamp(kidx) = Kofiko_CenterY(i);
    Kofiko_CenterY_forEachPlexonSample(pidx) = Kofiko_CenterY(i);
end

% Convert raw eye signals to screen coordinates (pixels) by:
% [RawSignal-SignalOffset]*Gain + ScreenCenter = FixationPoint
% SignalOffset is CenterX or CenterY

Kofiko_Xpix = Kofiko_GainX_forEachTimeStamp.*(Kofiko_EyeRawX - Kofiko_CenterX_forEachTimeStamp) + ScreenSizeX_pix/2;
Kofiko_Ypix = Kofiko_GainY_forEachTimeStamp.*(Kofiko_EyeRawY - Kofiko_CenterY_forEachTimeStamp) + ScreenSizeY_pix/2;

toc;
%% %%%%%%%%%%%%%% Trial analysis %%%%%%%%%%%%%%
tic;
fprintf('Extracting important trial variables\n')

% Find unique trials by finding unique Flip ON timestamps
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
tempImageFlipON_TS_Kofiko_unique = unique(tempImageFlipON_TS_Kofiko);
[~, uniqueTrialIdx, ~] = intersect(tempImageFlipON_TS_Kofiko, tempImageFlipON_TS_Kofiko_unique);
trialData= trialData(uniqueTrialIdx);
% Order trials by flip ON
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trialData);
[~, trial_idx_for_sorting] = sort(tempImageFlipON_TS_Kofiko , 'ascend');
trialData = trialData(trial_idx_for_sorting);

allFieldNames = fieldnames(trialData{1});

for i = 1:numel(allFieldNames)
    temp = cellfun(@(x) x.(allFieldNames{i}), trialData, 'UniformOutput', false);
    sz = cellfun(@size, temp, 'UniformOutput', false);
    sz(cellfun(@isempty, sz)) = [];
    sz(cellfun(@(x) any(x==0), sz)) = [];
    try
        sz = mode(vertcat(sz{:}), 1);
        temp(cellfun(@isempty, temp)) = {nan(sz)};
    catch
    end
    vars.(allFieldNames{i}) = temp;
end

% Determine stimulus intervals
stimStartTimes = [ones(size([vars.m_fImageFlipON_TS_Kofiko{:}]')),...
    [vars.m_fImageFlipON_TS_Kofiko{:}]']*B;

stimStopTimes = stimStartTimes + [vars.m_fStimulusON_MS{:}]'/1e3;

stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames =  min([vars.numFrames{:}]', [vars.numFrames{:}]' ./ [vars.repframes{:}]');

% Expand relevant variables across frame

% Get trial numbers
trialNumPerFrame = cellfun(@(x, y) repelem(x, y), vars.TrialNum, num2cell(numFrames), 'UniformOutput', false);

% Get trial IDs (insane that there are trial numbers and IDs....)
TrialIDPerFrame = cellfun(@(x, y) repelem(x, y), vars.TrialID, num2cell(numFrames), 'UniformOutput', false);

useBinaryPerFrame = cellfun(@(x, y) repelem(x, y), vars.usebinary, num2cell(numFrames), 'UniformOutput', false);

try
    spatialscale = [vars.spatialscale{:}];
catch
    spatialscale_buffer = g_astrctAllParadigms{1}.DualstimScale.Buffer;
    spatialscale_TS = g_astrctAllParadigms{1}.DualstimScale.TimeStamp;
    spatialscale = nan(size(trialData));
    for i = 1:numel(spatialscale_buffer)
        spatialscale(ImageFlipON_TS_Kofiko >= spatialscale_TS(i)) = spatialscale_buffer(i);
    end

end

BlockIDPerFrame =  cellfun(@(x, y) repelem(x, y), vars.BlockID, num2cell(numFrames), 'UniformOutput', false);

% Matrix of unique conditions defined by cloud parameters and block
uniqueCloudConditions =  unique([[vars.usebinary{:}]' spatialscale' [vars.BlockID{:}]'], 'rows');
uniqueCloudConditions(any(isnan(uniqueCloudConditions),2),:) = [];

DualstimPrimaryuseRGBCloudPerFrame = cellfun(@(x, y) repelem(x, y), vars.DualstimPrimaryuseRGBCloud, num2cell(numFrames), 'UniformOutput', false);

StimulusAreaPerFrame =  cellfun(@(x, y) repelem(x, y), vars.m_aiStimulusArea, num2cell(numFrames), 'UniformOutput', false);

stimseq= cellfun(@(x,y) x(1:y:end), vars.stimseq, num2cell(cellfun(@(x) max(0, x), vars.repframes)), 'UniformOutput', false);

if stimET
    stimseqET = cellfun(@(x,y) x(1:y:end), vars.stimseq_ET_Cclouds, num2cell(cellfun(@(x) max(0, x), vars.repframes)), 'UniformOutput', false);

end
%s

X_fixationSpot = cellfun(@(x) x(1), vars.m_pt2iFixationSpot, 'UniformOutput',false);
Y_fixationSpot = cellfun(@(x) x(2), vars.m_pt2iFixationSpot, 'UniformOutput',false);

UseLeyePerFrame = cellfun(@(x, y) repelem(x, y), vars.UseLeye, num2cell(numFrames), 'UniformOutput', false);
UseReyePerFrame = cellfun(@(x, y) repelem(x, y), vars.UseReye, num2cell(numFrames), 'UniformOutput', false);

toc;

%% Bin Kofiko eye signal timestamps by stimulus intervals (i.e., by trial)
[~,~,Kofiko_ET_TS_PlexonTime_Bin] = histcounts(Kofiko_ET_TS_PlexonTime, stimIntervals);
%stimON_Kofiko_ET_TS_PlexonTime_Bin = uniqueKofiko_ET_TS_PlexonTime_Bin(1:2:end);

% I use accumarray to create cell array where even numbered cells give eye
% signal timestamps for stimulus on periods (stimulus onset to offset),
% and odd numbered (starting at 3) cells give eye signal timestamps between stimulus on
% periods (stimulus offset to subsequent onset).The first cell corresponds
% to invalid timestamps (i.e. outside of stimIntervals)

nBins = numel(stimIntervals) - 1;

Kofiko_ET_TS_PlexonTime_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:) + 1, ...
    Kofiko_ET_TS_PlexonTime(:), ...
    [nBins + 1, 1], ...   % force size (extra 1 for bin 0 → invalid)
    @(x){x}, ...
    {[]} ...              % fill empty bins with empty cells
    );

% Do as above for eye signal X and Y coordinates:
Kofiko_Xpix_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:)+1, ...
    Kofiko_Xpix(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});

Kofiko_Ypix_cellArray = accumarray( ...
    Kofiko_ET_TS_PlexonTime_Bin(:)+1, ...
    Kofiko_Ypix(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});


%% Determine trials with good fixation
goodFixationX = cellfun(@(x, x_fix) sum(abs(x - x_fix) < maxFixationErrorPix)./numel(x)...
    > minFixationDuration, Kofiko_Xpix_cellArray(2:2:end), X_fixationSpot) ;

goodFixationY = cellfun(@(y, y_fix) sum(abs(y - y_fix) < maxFixationErrorPix)./numel(y) ...
    > minFixationDuration, Kofiko_Ypix_cellArray(2:2:end), Y_fixationSpot) ;

goodFixationIdx = vertcat(vars.m_bMonkeyFixated{:}) | (goodFixationX & goodFixationY);

%% Bin stimulus sequences by trial
tic;
fprintf('Loading stimuli\n')
stimulus_cellArray = cell(size(trialData,1),1);
stimulusET_cellArray = cell(size(trialData,1),1);

for i = 1:size(uniqueCloudConditions,1)
    thisUsebinary = uniqueCloudConditions(i,1);
    thisSpatialscale = uniqueCloudConditions(i,2);
    thisBlockID = uniqueCloudConditions(i,3);

    trialIdx = [vars.usebinary{:}]' == thisUsebinary & spatialscale' == thisSpatialscale & [vars.BlockID{:}]' == thisBlockID;
    if thisUsebinary == 0 % full contrast
        load(fullfile(stimpath, sprintf('Cloudstims_Chrom_size60_scale%d_%02d.mat', thisSpatialscale, thisBlockID)));
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
        if stimET
            stimulusET_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseqET(trialIdx), 'UniformOutput', false);
        end
    elseif thisUsebinary == 1
        load(fullfile(stimpath,sprintf('Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat', thisSpatialscale, thisBlockID)))
        DensenoiseChromcloud_DKlspace=int8(127*(DensenoiseChromcloud_DKlspace));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
        if stimET
            stimulusET_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseqET(trialIdx), 'UniformOutput', false);
        end
    elseif thisUsebinary == 2 % matched contrast
        load(fullfile(stimpath,sprintf('Cloudstims_ContrastMatched_size60_scale%d_%02d.mat', thisSpatialscale, thisBlockID)));
        DensenoiseChromcloud_DKlspace(:,:,:,2:3)=int8(127*(DensenoiseChromcloud_DKlspace(:,:,:,2:3)));
        DensenoiseChromcloud_DKlspace(:,:,:,1)=int8((1/LumScale)*127*(DensenoiseChromcloud_DKlspace(:,:,:,1)));
        stimulus_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseq(trialIdx), 'UniformOutput', false);
        if stimET
            stimulusET_cellArray(trialIdx) = cellfun(@(x) DensenoiseChromcloud_DKlspace(:,:,x,:), stimseqET(trialIdx), 'UniformOutput', false);
        end
    end

    % reshape stimuli into matrix where each row is stimulus, each col is frame

    % make first dim frames
    %stimulus_cellArray(trialIdx) = cellfun(@(x) permute(x, [3 1 2 4]), stimulus_cellArray(trialIdx), 'UniformOutput', false);
    %stimulus_cellArray(trialIdx) = cellfun(@(x) reshape(x, [size(x,1), prod(size(x,2:4))]),stimulus_cellArray(trialIdx), 'UniformOutput', false);
    %stimulus_cellArray(trialIdx) = cellfun(@transpose, stimulus_cellArray(trialIdx), 'UniformOutput', false);

    % % modified by dab
    tI2 = find(trialIdx == 1);
    stimulus_cellArray(tI2) = cellfun(@(x) permute(x, [3 1 2 4]), stimulus_cellArray(tI2), 'UniformOutput', false);
    stimulusET_cellArray(tI2) = cellfun(@(x) permute(x, [3 1 2 4]), stimulusET_cellArray(tI2), 'UniformOutput', false);
    for ii = 1:length(tI2)
        sz = size(stimulus_cellArray{tI2(ii)});
        stimulus_cellArray{tI2(ii)} = reshape( stimulus_cellArray{tI2(ii)}, [sz(1), prod(sz(2:4))] );
        if stimET
            stimulusET_cellArray{tI2(ii)} = reshape( stimulusET_cellArray{tI2(ii)}, [sz(1), prod(sz(2:4))] );
        end
    end
    stimulus_cellArray(tI2) = cellfun(@transpose, stimulus_cellArray(tI2), 'UniformOutput', false);
    if stimET
        stimulusET_cellArray(tI2) = cellfun(@transpose, stimulusET_cellArray(tI2), 'UniformOutput', false);
    end

    % stimulus_cellArray(trialIdx) = cellfun(@(x) permute(x, [3 1 2 4]), stimulus_cellArray(trialIdx), 'UniformOutput', false);
    % stimulus_cellArray(trialIdx) = cellfun(@(x) reshape(x, [size(x,1), prod(size(x,2:4))]),stimulus_cellArray(trialIdx), 'UniformOutput', false);
    % stimulus_cellArray(trialIdx) = cellfun(@transpose, stimulus_cellArray(trialIdx), 'UniformOutput', false);


    %
end

% find dualstim elemetns
trialTypeOfInterest = 'Dual Stim';
trialTypeOfInterestIdx = strcmpi( vars.m_strTrialType, trialTypeOfInterest);

trlonset_diffs = [4; diff(stimStartTimes)];
CCidx = cellfun(@(x) any(x==8), vars.DualstimPrimaryuseRGBCloud);
areaOverZeroIdx = cellfun(@(x) x>0, vars.m_aiStimulusArea);

isTrialOfInterest = trialTypeOfInterestIdx & ...
    goodFixationIdx &...
    trlonset_diffs > 4 &...
    CCidx & ...
    areaOverZeroIdx & ...
    vertcat(vars.DualstimPrimaryuseRGBCloud{:}) == 8 ; % clouds

stimulus_matrix = horzcat(stimulus_cellArray{isTrialOfInterest});
stimulusET_matrix = horzcat(stimulusET_cellArray{isTrialOfInterest});

stimFrameBins = cellfun(@(x, y, z) linspace(x,y,z+1), num2cell(stimStartTimes), num2cell(stimStopTimes), num2cell(numFrames), 'UniformOutput', false);
valid_stimFrameBinIdx = cellfun(@(x) numel(x)>=2, stimFrameBins);
toc;


%%
[~, wfcounts, ~, ~] = plx_info(plexon_fname, true);

[onlineSortedUnits, onlineSortedChans] = find(wfcounts);

% for each channel, and for each unit, make vector of timestamps and vector
% of unit id

onlineUnitID = 1;
online_spk_times = [];
online_spk_clusters = [];

uniqueOnlineSortedChans = unique(onlineSortedChans);

for ch = 1:numel(uniqueOnlineSortedChans)
    onlineUnits = onlineSortedUnits(onlineSortedChans == uniqueOnlineSortedChans(ch));
    for u = 1:numel(onlineUnits)
        [n, npw, ts, wave] = plx_waves_v(plexon_fname,uniqueOnlineSortedChans(ch)-1, onlineUnits(u));
        online_spk_times = [online_spk_times; ts];
        online_spk_clusters = [online_spk_clusters; onlineUnitID*ones(size(ts))];
        online_mean_wave{ch,u} = mean(wave,1);
        onlineUnitID = onlineUnitID + 1;
    end

end

online_allUnit_clusterIDs = unique(online_spk_clusters);

% bin spikes by trial
[~,~,online_spk_times_bin] = histcounts(online_spk_times, stimIntervals);

% each cell gives spike times for stim ON and stim OFF periods

online_spk_times_cellArray =accumarray(...
    online_spk_times_bin(:)+1, ...
    online_spk_times(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});

online_clusterIDForEachSpk_cellArray = accumarray(...
    online_spk_times_bin(:)+1, ...
    online_spk_clusters(:), ...
    [nBins + 1, 1], ...
    @(x){x}, ...
    {[]});

online_clusterIDForEachSpk_stimON_cellArray = online_clusterIDForEachSpk_cellArray(2:2:end);
online_spk_times_stimON_cellArray =  online_spk_times_cellArray(2:2:end);
%    ks_batchForEachSpk_stimON_cellArray = ks_batchForEachSpk_stimON_cellArray{ks_batch}(2:2:end);

online_stimFrameNumForEachSpk_cellArray = cell(size(trialData,1),1);

[~,~,online_stimFrameNumForEachSpk_cellArray(valid_stimFrameBinIdx)] =...
    cellfun(@(x,y) histcounts(x, 'BinEdges', y), ...
    online_spk_times_stimON_cellArray(valid_stimFrameBinIdx),...
    stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
%spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

online_Robs = single(zeros(numel(unique(online_spk_clusters)),  sum(numFrames(isTrialOfInterest))));

for unit = 1:numel(online_allUnit_clusterIDs)
    unitID = online_allUnit_clusterIDs(unit);

    % each cell is a trial, and contains frame indicies where spikes
    % occurred (repeated frames == multiple spikes on that frame)

    online_stimFrameNumForEachSpk_thisUnit_cellArray = cell(size(trialData,1),1);
    online_stimFrameNumForEachSpk_thisUnit_cellArray(valid_stimFrameBinIdx) = ...
        cellfun(@(x,y) x(y==unitID), ...
        online_stimFrameNumForEachSpk_cellArray(valid_stimFrameBinIdx), ...
        online_clusterIDForEachSpk_stimON_cellArray(valid_stimFrameBinIdx),...
        'UniformOutput', false);

    online_spksPerFrame_cellArray =  cell(size(trialData,1),1);

    online_spksPerFrame_cellArray(isTrialOfInterest) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
        online_stimFrameNumForEachSpk_thisUnit_cellArray(isTrialOfInterest), num2cell(numFrames(isTrialOfInterest)) ,'UniformOutput',false);

    online_Robs(unit,:) = [online_spksPerFrame_cellArray{:}];
end

online_tempSTA = [];

S = single(stimulus_matrix)';
tic;
fprintf('Computing STAs\n');

for lag = 0:nLags-1
    online_tempSTA(:,:,:,:,lag+1) = online_Robs(:,lag+1:end) * S(1:end-lag,:);
end

online_STA = reshape(online_tempSTA, size(online_tempSTA,1), 60, 60, 3, nLags);

toc;



lags = 3:6;
chrom_chan_labels = {'Lum.', 'L-M', 'S'};
%if compute_stas && plotting
fprintf('Now plotting\n')

for unit = 1:size(online_STA,1)
    figure; i = 1;
    for chromatic_channel = 1:3
        for lag = lags
            subplot(3, length(lags), i );
            imagesc(circshift(squeeze(online_STA(unit,:,:,chromatic_channel,lag+1)), [30 30])), colormap gray
            i = i+1;

            if lag == lags(1)
                ylabel(chrom_chan_labels{chromatic_channel}, 'fontsize', 16)
            end

            if chromatic_channel == 1
                title(lag)
            end
        end
    end

    sgtitle(sprintf('channel: %i,unit: %i, total spikes: %i', onlineSortedChans(unit),onlineSortedUnits(unit), sum(online_Robs(unit,:))))
end
% end
