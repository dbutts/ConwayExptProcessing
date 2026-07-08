%%

% 3/2026 mjg -- wrote it
% cell arrays:

% For cell arrays of size (2*# trials) x 1, even numbered cells
% correspond to stimulus ON periods, odd numbered cells (starting at 3)
% correspond to interstimulus periods.

% Kofiko_ET_TS_PlexonTime_cellArray [(2*# trials) x 1]: Kofiko eye tracking timestamps. Even numbered
% cells correspond to stimulus on periods, odd numbered cells starting at 3
% correspond to interstimulus periods.

% Kofiko_Xpix_cellArray [(2*# trials) x 1]: Calibrated Kofiko horizontal eye trace (screen coordinates)

% Kofiko_Ypix_cellArray [(2*# trials) x 1]: Calibrated Kofiko vertical eye trace (screen coordinates)

% t_plexon_cellArray [(2*# trials) x 1]: Plexon sample times

% leftEyePupil_plexon_cellArray [(2*# trials) x 1]: Plexon left eye pupil
% trace

% rightEyePupil_plexon_cellArray: same as above for right eye.

% rightEyeX_plexon_calib_cellArray [(2*# trials) x 1]: Plexon right eye
% horizontal trace with Kofiko gains and offsets applied

% rightEyeY_plexon_calib_cellArray [(2*# trials) x 1]: same as above for
% vertical trace.

% leftEyeX_plexon_calib_cellArray [(2*# trials) x 1]: Plexon left eye
% horizontal trace with Kofiko gains and offsets applied

% leftEyeY_plexon_calib_cellArray [(2*# trials) x 1]: same as above for
% vertical trace.

% stimulus_cellArray [# trials x 1]: Each cell is a (stimulus width*stimulus
% height*# chromatic channels) x # frames matrix that contains the cloud
% pixel values (-127 to 127).

% Kofiko_Xpix_frameRate_cellArray [# trials of interest x 1]: Calibrated Kofiko
% horizontal eye trace downsampled to frame rate.

% Kofiko_Xpix_frameRate_cellArray [# trials of interest x 1]: Calibrated Kofiko
% vertical eye trace downsampled to frame rate.

% spk_times_cellArray [# spike sorting batches x 1]: Each cell contains a
% (2*# trials) x 1 cell array. Even cells give spikes that occurred during stimulus ON.
% Odd trials starting at 3 give spikes during interstimulus periods.

% clusterIDForEachSpk_cellArray [# spike sorting batches x 1]: Each cell contains a
% (2*# trials) x 1 cell array. Even cells give IDs corresponding to spikes that occurred during stimulus ON.
% Odd trials starting at 3 give IDs during interstimulus periods.

% ks_batchForEachSpk_cellArray [# spike sorting batches x 1]: Each cell
% contains a (2*# trials) x 1 cell array. Even cells give spike sorting batch corresponding to spikes that occurred during stimulus ON.
% Odd trials starting at 3 give batch during interstimulus periods.

%  stimFrameNumForEachSpk_cellArray  [# spike sorting batches x 1]: Each cell
% contains a (2*# trials) x 1 cell array. Even cells give frames shown during spikes that occurred during stimulus ON.
% Odd trials starting at 3 give frames during interstimulus periods.

% spksPerFrame_cellArray [# trials x 1]: Each cell contains # frames x 1
% vector giving spike counts per frame.

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

end

%% flags
fprintf('Flags \n')
fprintf('--------------\n')

% the one flag Dan needs (for now)
skipLFP = input(['skip LFP? \n' ...
    '[1] Yes \n' ...
    '[0] No \n ']); %1;

% This only needs to run if there is one -- so should be detected
% you don't have any except old experiments, so I'll just set this based on
% config

% print defaults
% I don't need all these prompts so will just move the ones I do out and
% all the rest can be for how you want to enter information


if computerLocation < 10   % NIH-specific settings and prompts
    
    % is there a situation where we do not want to save the data? I'll
    saving = input(['save data? \n' ...
        '[1] Yes \n' ...
        '[0] No \n ']); %1;


    online_sorting = input('online sorting? \n [0] No \n [1] Yes \n' );

    stimET = input(['stimET? \n' ...
        '[1] Yes \n' ...
        '[0] No \n ']); %1;

    fprintf('Set paths \n')
    fprintf('--------------\n')

    manualFnames = input('Input full file paths manually? \n [0] No \n [1] Yes \n');

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
    ks_path = input('ks_path: ', 's');

    monkey_name = input('monkey_name: ', 's');
    rig = input('Rig: ', 's');

    bino_ddpi = input('Binocular ddpi? \n [0] No \n [1] Yes \n');
    if ~bino_ddpi
        ET_eyelink = input('Specifiy eye tracker \n [0] Eyescan \n [1] Mono eyelink \n [2] Bino eyelink \n [3] Mono ddpi \n');
    end

    % save time with these assertions
    %assert(isdir(dirpath) & isdir(ks_path) & isdir(stimpath), 'Check that dirpath, ks_path, and stimpath exist!');

    [~, filenameP] = fileparts(mainKofiko_fname);
    exptdate  = str2num(filenameP(1:6));

    if saving
        savepath = input('savepath: ', 's');
    else
        savepath = '';
    end

else % Dan's lab specific settings and prompts

    stimET = 1;   % need for old experiments (if included)
    online_sorting = 0;
    saving = 1;
    rig = 'C'; % I don't know what this does

    % Enter date of experiment and figure out the filenameP
    exptdate = input( 'Enter date of experiment (YYMMDD): ' ); % this makes integer date
    files = dir(sprintf('%s*%d*', dirpath, exptdate));
    names = string({files.name});
    if isempty(names)
        if computerLocation == 11
            disp('On MST: make sure legacy data drives are mounted.')
        end
        error('Experiment with that date not found.')
    end
    lens = strlength(names);
    [~, idxMin] = min(lens);
    filenameP = char(names(idxMin));
    fprintf('Found %s\n', filenameP)
    plexon_fname = [dirpath filenameP '.pl2'];
    kofiko_subfolder = [dirpath filenameP '/'];
    ks_path = kofiko_subfolder;
    mainKofiko_fname = [dirpath filenameP '.mat'];
end

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
%else
%    all_names = {kofiko_dir.name};
%    kofiko_fname_re = ['^' filenameP '(?:_\d+)?\.mat$'];
%    matches = regexp(all_names, kofiko_fname_re, 'once');
%    valid_kofiko_fname_idx = ~cellfun(@isempty, matches);

%    kofiko_folder = {kofiko_dir(valid_kofiko_fname_idx).folder};
%    kofiko_fname  = {kofiko_dir(valid_kofiko_fname_idx).name};

% Get the Kofiko mat file numbers from their names
%    kofiko_file_numStr = cellfun(@(x) regexp(x, '\d*(?=(.mat))', 'match'), kofiko_fname, 'UniformOutput',false);
%    if sum(cellfun(@isempty, kofiko_file_numStr)) > 0
%    	kofiko_file_numStr{cellfun(@isempty, kofiko_file_numStr)} = {'nan'};
%    end
%    kofiko_file_num = cellfun(@(x) str2num(x{:}), kofiko_file_numStr);
%    kofiko_file_num(isnan(kofiko_file_num)) = max(kofiko_file_num)+1;
%end


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

%% Plexon eye data

if strcmpi(rig, 'B')
    % Rig B
    sync_ch = 'AI01';
    arc_ch = 'AI02';
    leftEyePupil_ch = 'AI07';
    rightEyePupil_ch = 'AI08';
    rightEyeX_ch = 'AI05';
    rightEyeY_ch = 'AI06';
    leftEyeX_ch ='AI03';
    leftEyeY_ch ='AI04';
elseif strcmpi(rig, 'C')
    % % Rig C
    if bino_ddpi
        sync_ch = 'AI01';
        arc_ch = 'AI02';
        leftEyePupil_ch = 'AI03';
        rightEyePupil_ch = 'AI04';
        rightEyeX_ch = 'AI05';
        rightEyeY_ch = 'AI06';
        leftEyeX_ch ='AI07';
        leftEyeY_ch ='AI08';
    else
        if ET_Eyelink == 1
            [~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI05');
            [~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI06');
            [~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
            [ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
            PlexET_ad_calib=PlexET_ad;
            PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./opts.plx_analogscale);
            PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./opts.plx_analogscale);
            PlexET_ad_calib(3,:) = (PlexET_ad_calib(3,:)-median(PlexET_ad_calib(3,:)))*(g_strctEyeCalib.GainX.Buffer(end)./opts.plx_analogscale);
            PlexET_ad_calib(4,:) = (PlexET_ad_calib(4,:)-median(PlexET_ad_calib(4,:)))*(g_strctEyeCalib.GainY.Buffer(end)./opts.plx_analogscale);

        elseif ET_Eyelink== 2
            [~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
            [ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
            PlexET_ad_calib=PlexET_ad;
            PlexET_ad_calib(1,:) = (PlexET_ad_calib(1,:)-median(PlexET_ad_calib(1,:)))*(g_strctEyeCalib.GainX.Buffer(end)./opts.plx_analogscale);
            PlexET_ad_calib(2,:) = (PlexET_ad_calib(2,:)-median(PlexET_ad_calib(2,:)))*(g_strctEyeCalib.GainY.Buffer(end)./opts.plx_analogscale);

        elseif ET_Eyelink == 3
            [~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI03');
            [~, ~, ~, ~, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI04');
            [~, ~, ~, ~, PlexET_ad(3,:)] = plx_ad_v(thisSessionFile, 'AI07');
            [ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(4,:)] = plx_ad_v(thisSessionFile, 'AI08');
            PlexET_ad_calib=PlexET_ad;
            PlexET_ad_calib(3,:) = PlexET_ad_calib(3,:)*(g_strctEyeCalib.GainX.Buffer(end)./opts.plx_analogscale);
            PlexET_ad_calib(4,:) = PlexET_ad_calib(4,:)*(g_strctEyeCalib.GainY.Buffer(end)./opts.plx_analogscale);

        elseif ET_Eyelink == 0
            [~, ~, ~, ~, PlexET_ad(1,:)] = plx_ad_v(thisSessionFile, 'AI07');
            [ET_adfreq, ET_n, ET_ts, ET_fn, PlexET_ad(2,:)] = plx_ad_v(thisSessionFile, 'AI08');
            PlexET_ad(1,:)=PlexET_ad(1,:)-median(PlexET_ad(1,:));
            PlexET_ad(2,:)=PlexET_ad(2,:)-median(PlexET_ad(2,:));
        end
    end

end

[adfreq, n, ts, fn, sync_ad] = plx_ad_v(plexon_fname, sync_ch);
t_plexon = (0:n-1)/adfreq;

[~, ~, ~, arc_ad] = plx_ad_v(plexon_fname, arc_ch);
[~, ~, ~, ~, leftEyePupil_plexon] = plx_ad_v(plexon_fname, leftEyePupil_ch);
[~, ~, ~, ~, rightEyePupil_plexon] = plx_ad_v(plexon_fname, rightEyePupil_ch);
[~, ~, ~, ~, rightEyeX_plexon] = plx_ad_v(plexon_fname,rightEyeX_ch);
[~, ~, ~, ~, rightEyeY_plexon] = plx_ad_v(plexon_fname, rightEyeY_ch);
[~, ~, ~, ~, leftEyeX_plexon] = plx_ad_v(plexon_fname, leftEyeX_ch);
[~, ~, ~, ~, leftEyeY_plexon] = plx_ad_v(plexon_fname, leftEyeY_ch);



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

%% Apply calibration to plexon signals
rightEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*rightEyeX_plexon';
rightEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*rightEyeY_plexon';
leftEyeX_plexon_calib = plexonAnalogScale.*Kofiko_GainX_forEachPlexonSample.*leftEyeX_plexon';
leftEyeY_plexon_calib = plexonAnalogScale.*Kofiko_GainY_forEachPlexonSample.*leftEyeY_plexon';

% create PlexET_ad structure
if bino_ddpi
PlexET_times = t_plexon;
PlexET_ad_calib = [];
PlexET_ad_calib(:,1) = sync_ad;
PlexET_ad_calib(:,2) = arc_ad;
PlexET_ad_calib(:,3) = leftEyePupil_plexon;
PlexET_ad_calib(:,4) = rightEyePupil_plexon;
PlexET_ad_calib(:,5) = rightEyeX_plexon_calib;
PlexET_ad_calib(:,6) = rightEyeY_plexon_calib;
PlexET_ad_calib(:,7) = leftEyeX_plexon_calib;
PlexET_ad_calib(:,8) = leftEyeY_plexon_calib;
end


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

%% %%%%%%%%%%%%% Extract fixinfo %%%%%%%%%%%%%
calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), vars.m_strTrialType);

ETdata.fixloc = vertcat(vars.m_pt2iFixationSpot{isCalibrationTrial});
ETdata.plxonset = stimStartTimes(isCalibrationTrial);

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

%% Bin Plexon eye signal by stimulus intervals
[~,~, t_plexon_bin] = histcounts(t_plexon, stimIntervals);

t_plexon_cellArray = accumarray(t_plexon_bin(:)+1, t_plexon(:), [], @(x){x});

% pupils
leftEyePupil_plexon_cellArray = accumarray(t_plexon_bin(:)+1, leftEyePupil_plexon(:), [], @(x){x});
rightEyePupil_plexon_cellArray = accumarray(t_plexon_bin(:)+1, rightEyePupil_plexon(:), [], @(x){x});

% right eye
rightEyeX_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, rightEyeX_plexon_calib(:), [], @(x){x});
rightEyeY_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, rightEyeY_plexon_calib(:), [], @(x){x});

% left eye
leftEyeX_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, leftEyeX_plexon_calib(:), [], @(x){x});
leftEyeY_plexon_calib_cellArray = accumarray(t_plexon_bin(:)+1, leftEyeY_plexon_calib(:), [], @(x){x});

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

%% resample eye signal at frame rate
% NOTE: Without 'extrap' argument, there will be NANs:

Kofiko_Xpix_frameRate_cellArray = ...
    cellfun(@(t,x, start, stop, fixspot, n) interp1(t, x - fixspot, linspace(start, stop, n), 'linear'),...
    Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    Kofiko_Xpix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimStartTimes(isTrialOfInterest)),...
    num2cell(stimStopTimes(isTrialOfInterest)),...
    X_fixationSpot(isTrialOfInterest),...
    num2cell(numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);

Kofiko_Ypix_frameRate_cellArray = ...
    cellfun(@(t,y, start, stop, fixspot, n) interp1(t, y - fixspot, linspace(start, stop, n), 'linear'),...
    Kofiko_ET_TS_PlexonTime_cellArray(2*find(isTrialOfInterest)),...
    Kofiko_Ypix_cellArray(2*find(isTrialOfInterest)),...
    num2cell(stimStartTimes(isTrialOfInterest)),...
    num2cell(stimStopTimes(isTrialOfInterest)),...
    Y_fixationSpot(isTrialOfInterest),...
    num2cell(numFrames(isTrialOfInterest)), ...
    'UniformOutput', false);

%% load in spike data
tic;
fprintf('Loading spike data and computing Robs\n');
% drop test check

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
    opts.spk_offset = droptestcheck;
else
    opts.spk_offset = 0;
end

disp('Drop test check complete')
% end drop test check %%%%
%% cluster label strings

goodStr = 'good ';
muaStr = 'mua  ';
noiseStr = 'noise';
blankStr = '     ';

%% %%%%%%%%%%%%% Load and organize spike data %%%%%%%%%%%%%
% if using kilosort
% Find folders with kilosort output
spike_times_dir = dir(fullfile(ks_path, '**/spike_times.npy'));
spike_clusters_dir = dir(fullfile(ks_path, '**/spike_clusters.npy'));
cluster_info_dir = dir(fullfile(ks_path, '**/cluster_info.tsv'));

spike_times_folders = {spike_times_dir(:).folder};
spike_clusters_folders = {spike_clusters_dir(:).folder};
cluster_info_folders = {cluster_info_dir(:).folder};

% For each folder with kilosort outputs
num_ks_batch = length(spike_times_dir);
assert(~isempty(spike_times_dir));

[~, ks_folders, ~] =  cellfun(@fileparts, spike_times_folders, 'UniformOutput', false);
tokens = regexp(ks_folders, '^[^_]+_([^_]+)', 'tokens');
array_labels = cellfun(@(t) t{1}{1}, tokens, 'UniformOutput', false);
%unique_array_labels = unique(array_labels);
chan_offset = 0;
cluster_offset = 0;
%chan_offsets(1) = chan_offset;

nSU = [];
nMU = [];
RobsSU =[];
RobsMU =[];

for ks_batch = 1:num_ks_batch
    this_array_label = array_labels{ks_batch};
    disp(ks_batch)
    try
        % Get label of kilosort batch (often corresponding to array name and range of channels processed)

        % Read in kilosort outputs
        spk_times = readNPY(fullfile(spike_times_folders{ks_batch}, 'spike_times.npy'));
        spk_times = double(spk_times)./fs + opts.spk_offset; % convert to seconds
        spk_clusters = readNPY(fullfile(spike_times_folders{ks_batch}, 'spike_clusters.npy'))...
            + cluster_offset;
        cluster_KSLabel = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_KSLabel.tsv'));

        chan_map = readNPY(fullfile(spike_times_folders{ks_batch}, 'channel_map.npy'));
        chan_map = chan_map + chan_offset; % make channel numbers unique within array

        if isfile(fullfile(spike_times_folders{ks_batch}, 'cluster_info.tsv'))
            cluster_info = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_info.tsv'));
            cluster_id = cluster_info.cluster_id + cluster_offset;
            group = cluster_info.group;
            n_spikes = cluster_info.n_spikes;
            chan_best = cluster_info.ch + double(chan_offset);
        else
            cluster_group = tdfread(fullfile(spike_times_folders{ks_batch}, 'cluster_group.tsv'));
            cluster_group.cluster_id = cluster_group.cluster_id + cluster_offset;

            % account for blank units which may not be in cluster_group
            blank_cluster_id = setdiff(unique(spk_clusters), cluster_group.cluster_id);
            temp_cluster_id = [cluster_group.cluster_id; blank_cluster_id];
            temp_group = cluster_group.group;
            temp_group(end+1:end+length(blank_cluster_id),:) = ' ';
            [temp_cluster_id_sorted, I] = sort(temp_cluster_id);
            temp_group_sorted = temp_group(I,:);
            cluster_id = temp_cluster_id_sorted;
            group = temp_group_sorted;

            % find best channel of each cluster
            templates = readNPY(fullfile(spike_times_folders{ks_batch}, 'templates.npy'));
            % n_spikes = accumarray(spk_clusters+1, spk_clusters, [], @numel);
            % n_spikes = n_spikes(n_spikes>0);

            n_spikes = accumarray(spk_clusters+1, 1, [], @sum);
            n_spikes = n_spikes(cluster_id + 1);

            [~,I]= max(sum(templates.^2,2),[],3);
            chan_best = chan_map(I); % best channel for each unique cluster
        end

        % Find indices of units labeled "good", "mua", or ""
        isGood = cellfun(@(x) strcmpi(deblank(x), 'good'), cellstr(group));
        isMua = cellfun(@(x) strcmpi(deblank(x), 'mua'), cellstr(group));
        isBlank = cellfun(@(x) isempty(deblank(x)), cellstr(group));
        hasMinSpikes = n_spikes > minSpikes;

        % get cluster IDs
        SU_clusterIDs{ks_batch} = cluster_id(isGood & hasMinSpikes);
        MU_clusterIDs{ks_batch} = cluster_id((isMua | isBlank) & hasMinSpikes);
        allUnit_clusterIDs{ks_batch} =[SU_clusterIDs{ks_batch}; MU_clusterIDs{ks_batch}];

        % get channels
        SU_chanNums{ks_batch} = chan_best(isGood & hasMinSpikes);
        MU_chanNums{ks_batch} = chan_best((isMua | isBlank) & hasMinSpikes);
        allUnit_chanNums{ks_batch} =  [SU_chanNums{ks_batch}; MU_chanNums{ks_batch}];

        % get ks_batch
        SU_ks_batch{ks_batch} = ks_batch.*ones(size(SU_clusterIDs{ks_batch}));
        MU_ks_batch{ks_batch} = ks_batch.*ones(size(MU_clusterIDs{ks_batch}));
        allUnit_ks_batch{ks_batch} = [SU_ks_batch{ks_batch}; MU_ks_batch{ks_batch}];

        % Get rid of spike times and clusterIDs that correspond to bad units
        spk_times = spk_times(ismember(spk_clusters, allUnit_clusterIDs{ks_batch}));
        spk_clusters = spk_clusters(ismember(spk_clusters, allUnit_clusterIDs{ks_batch}));
        ks_batchForEachSpk= repelem(ks_batch, numel(spk_times))';

        % bin spikes by trial
        [~,~,spk_times_bin{ks_batch}] = histcounts(spk_times, stimIntervals);

        % each cell gives spike times for stim ON and stim OFF periods

        spk_times_cellArray{ks_batch} =accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            spk_times(:), ...
            [nBins + 1, 1], ...
            @(x){x}, ...
            {[]});

        clusterIDForEachSpk_cellArray{ks_batch} = accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            spk_clusters(:), ...
            [nBins + 1, 1], ...
            @(x){x}, ...
            {[]});

        ks_batchForEachSpk_cellArray{ks_batch} = accumarray(...
            spk_times_bin{ks_batch}(:)+1, ...
            ks_batchForEachSpk(:), ...
            [nBins + 1, 1], ...)
            @(x){x}, ...
            {[]});

        clusterIDForEachSpk_stimON_cellArray = clusterIDForEachSpk_cellArray{ks_batch}(2:2:end);
        spk_times_stimON_cellArray =  spk_times_cellArray{ks_batch}(2:2:end);
        %    ks_batchForEachSpk_stimON_cellArray = ks_batchForEachSpk_stimON_cellArray{ks_batch}(2:2:end);

        stimFrameNumForEachSpk_cellArray{ks_batch} = cell(size(trialData,1),1);

        [~,~,stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx)] =...
            cellfun(@(x,y) histcounts(x, 'BinEdges', y), ...
            spk_times_stimON_cellArray(valid_stimFrameBinIdx),...
            stimFrameBins(valid_stimFrameBinIdx), 'UniformOutput',false);
        %spkFrameIdx_cellArray = transpose(spkFrameIdx_cellArray);

        Robs{ks_batch} = single(zeros(numel(allUnit_clusterIDs{ks_batch}),  sum(numFrames(isTrialOfInterest))));

        for unit = 1:numel(allUnit_clusterIDs{ks_batch})
            unitID = allUnit_clusterIDs{ks_batch}(unit);

            % each cell is a trial, and contains frame indicies where spikes
            % occurred (repeated frames == multiple spikes on that frame)

            stimFrameNumForEachSpk_thisUnit_cellArray = cell(size(trialData,1),1);
            stimFrameNumForEachSpk_thisUnit_cellArray(valid_stimFrameBinIdx) = ...
                cellfun(@(x,y) x(y==unitID), ...
                stimFrameNumForEachSpk_cellArray{ks_batch}(valid_stimFrameBinIdx), ...
                clusterIDForEachSpk_stimON_cellArray(valid_stimFrameBinIdx),...
                'UniformOutput', false);

            spksPerFrame_cellArray =  cell(size(trialData,1),1);

            spksPerFrame_cellArray(isTrialOfInterest) = cellfun(@(x,y) histcounts(x, 'BinEdges',0.5:(y+0.5)),...
                stimFrameNumForEachSpk_thisUnit_cellArray(isTrialOfInterest), num2cell(numFrames(isTrialOfInterest)) ,'UniformOutput',false);

            Robs{ks_batch}(unit,:) = [spksPerFrame_cellArray{:}];
        end

    catch ME
        warning("Could not process %s", ks_folders{ks_batch})
        disp(ME.message)
    end

    nSU{ks_batch} = length(SU_clusterIDs{ks_batch});
    nMU{ks_batch} = length(MU_clusterIDs{ks_batch});
    RobsSU{ks_batch} = Robs{ks_batch}(1:nSU{ks_batch},:);
    RobsMU{ks_batch} = Robs{ks_batch}(nSU{ks_batch}+1:nSU{ks_batch}+nMU{ks_batch},:);

    % update chan_offset
    if ks_batch < num_ks_batch % if we still got a batch ahead
        next_array_label = array_labels{ks_batch+1};
        if strcmpi(next_array_label,  this_array_label)
            chan_offset = max(chan_map)+1; % maximum channel number;
        else
            chan_offset = 0; % new array, reset offset to 0
        end
        %   chan_offsets(ks_batch+1) = chan_offset;
    end
    % update cluster offset, and keep track of them
    cluster_offset = max(cluster_id)+1;
    cluster_offsets(ks_batch+1) = cluster_offset;
    chan_offsets(ks_batch+1) = chan_offset;

end

nSU = sum(cellfun(@length, SU_clusterIDs));
nMU = sum(cellfun(@length, MU_clusterIDs));

toc;

%% online sorting
if online_sorting == 1
    [tscounts, wfcounts, evcounts, contcounts] = plx_info(plexon_fname, true);

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
    if compute_stas
        S = single(stimulus_matrix)';
        tic;
        fprintf('Computing STAs\n');

        for lag = 0:nLags-1
            online_tempSTA(:,:,:,:,lag+1) = online_Robs(:,lag+1:end) * S(1:end-lag,:);
        end

        online_STA = reshape(online_tempSTA, size(online_tempSTA,1), 60, 60, 3, nLags);

        toc;
    else
    end

    lags = 3:6;

    chrom_chan_labels = {'Lum.', 'L-M', 'S'};
    if compute_stas && plotting
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
    end
end

%% Process LFPs
if ~skipLFP
    disp('Processing LFPs')
    tic;
    if (exptdate/10000) < 24  % year less than 2024
        LFPchans{1} = 1:24; % laminar
        LFPchans{2} = [33,40,46,47,52,53,54,59,65,67,71,81,83,89,90,95,98,102,103,109,112,131,138,139,145,146,152,158]; % Nform channels that worked
        LFPchans{3} = 161:256; % utah
    else
        LFPchans{1} = length(pl2.SpikeChannels); % one array so far
    end

    [LFP_adfreq, LFP_n, LFP_ts, ~, ~] = plx_ad_v(plexon_fname, ['FP' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
    for ii=1:length(LFPchans)
        nchans = length(LFPchans{ii});
        %LFP_ad = zeros(nchans, LFP_n);
        LFPs{ii} = zeros(nchans, LFP_n);
        for ch = 1:nchans
            %[~,~,~,~, LFP_ad(i,:)] = plx_ad_v(plexon_fname, ['FP' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
            [~,~,~,~, LFPs{ii}(ch,:)] = plx_ad_v(plexon_fname, ['FP' num2str(LFPchans{ii}(ch), ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);
        end
    end
    LFP_times=(0:LFP_n-1)/LFP_adfreq;
    toc;
end

%% %%%%%%%%%%%% Format output like PackageCloud %%%%%%%%%%%%
%allRobs = vertcat(Robs{:});
allRobsSU = vertcat(RobsSU{:});
allRobsMU = vertcat(RobsMU{:});
allRobs = vertcat(allRobsSU, allRobsMU);

%ETgains
ETgains = [Kofiko_GainX(end), Kofiko_GainY(end)];

%ETstim_location
ETstim_location =  [trialData{find(isTrialOfInterest,1, 'last')}.secondarystim_bar_rect;...
    trialData{find(isTrialOfInterest,1, 'last')}.tertiarystim_bar_rect];

% ETtrace
ETtrace = [[Kofiko_Xpix_frameRate_cellArray{:}]; [Kofiko_Ypix_frameRate_cellArray{:}]];

%ETtrace_raw
ETtrace_raw = transpose([vertcat(rightEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeX_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyeY_plexon_calib_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(rightEyePupil_plexon_cellArray{2*find(isTrialOfInterest)}),...
    vertcat(leftEyePupil_plexon_cellArray{2*find(isTrialOfInterest)})]);

%Robs
RobsSU = allRobs(1:nSU,:);

%RobsMU
RobsMU = allRobs(nSU+1:end,:);

%RobsMU_probe_ID
RobsMU_probe_ID = vertcat(MU_chanNums{:});

%RobsMU_rating
RobsMU_rating = [];
%blockID

% Robs_probe_ID
Robs_probe_ID = vertcat(SU_chanNums{:});

%Robs rating
Robs_rating =[];
blockID = [BlockIDPerFrame{isTrialOfInterest}];

%block_inds
block_inds = [cumsum(numFrames(isTrialOfInterest))' - (numFrames(find(isTrialOfInterest,1,'first'))-1); cumsum(numFrames(isTrialOfInterest))'];

%cloud_area
cloud_area = [StimulusAreaPerFrame{isTrialOfInterest}];

%cloud_binary
cloud_binary = [useBinaryPerFrame{isTrialOfInterest}];

%cloud_scale
cloud_scale = transpose(repelem(spatialscale(isTrialOfInterest), numFrames(isTrialOfInterest)));
% spatialScalePerFrame(isTrialOfInterest);

%datafilts
datafilts = ones(size(RobsSU));

%datafiltsMU
datafiltsMU = ones(size(RobsMU));

%dt
dt = 1/60.0; %0.0160; % why not 0.0167?

%electrode_info
electrode_info =[];

%exptdate and exptname
exptdate = filenameP(1:6);
exptname = filenameP;

%fix_location
fix_location = unique(vertcat(vars.m_pt2iFixationSpot{isTrialOfInterest}), 'rows');

%fix_size
fix_size = trialData{find(isTrialOfInterest, 1, 'last')}.m_fFixationSizePix -1;

%pixel_size
pixel_size = 1;

%sacc_inds
sacc_inds=[];

%ss
clusterIDs= cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), clusterIDForEachSpk_cellArray, 'UniformOutput', false);
clusterIDs = transpose(vertcat(clusterIDs{:}));
clusterIDs = vertcat(clusterIDs);

assert(numel(unique(clusterIDs)) == nSU + nMU);

%spike_ts
spike_ts_raw = cellfun(@(x) vertcat(x{2*find(isTrialOfInterest)}), spk_times_cellArray, 'UniformOutput', false);
spike_ts_raw = transpose(vertcat(spike_ts_raw{:}));
spike_ts_raw = vertcat(spike_ts_raw);

% this will get you spike times relative to trial start, ie in range 0 to 4:
trlsecs = unique([vars.m_fStimulusON_MS{isTrialOfInterest}])/1e3;
trialStart = stimStartTimes(isTrialOfInterest);
trialStop  = stimStopTimes(isTrialOfInterest);
spike_ts_raw = spike_ts_raw(:);
trialStart   = trialStart(:);
trialStop    = trialStop(:);

trialIdx = discretize(spike_ts_raw, [trialStart; inf]);

valid = ~isnan(trialIdx) & trialIdx >= 1 & trialIdx <= numel(trialStart);

valid2 = false(size(valid));
valid2(valid) = spike_ts_raw(valid) >= trialStart(trialIdx(valid)) & ...
    spike_ts_raw(valid) <= trialStop(trialIdx(valid));

valid = valid & valid2;

spike_ts_raw = spike_ts_raw(valid);
trialIdx     = trialIdx(valid);
spike_ts     = spike_ts_raw - trialStart(trialIdx) + (trialIdx - 1) * trlsecs;

%stim
stim = reshape(stimulus_matrix, 60,60,3,[]);

%stimET
if stimET
    stimET = reshape(stimulusET_matrix, 60,60,3,[]);
end

%stim_area
stim_area = [vars.m_aiStimulusArea{isTrialOfInterest}]';

%stim_location
modal_stim_area = mode([vars.m_aiStimulusArea{isTrialOfInterest}]');
if computerLocation < 10
    % this line does not work in Dan lab
    medianTiledStimulusRect = median(cat(3, vars.m_aiTiledStimulusRect{:}),3, 'omitmissing');
else
    % previous code that worked
    if isfield(trialData{1}, 'm_aiTiledStimulusRect')
        TiledStimulusRect = cellfun(@(x) x.m_aiTiledStimulusRect, trialData, 'UniformOutput', false);
    else
        StimulusRect = cellfun(@(x) x.m_aiStimulusRect, trialData, 'UniformOutput', false);
        TiledStimulusRect = StimulusRect; %%%%% Only for old datasets before tiling (e.g., 7/07/22)
    end
    medianTiledStimulusRect = median(cat(3, TiledStimulusRect{:}),3);
end

stim_location = [medianTiledStimulusRect(:,1:2) medianTiledStimulusRect(:,1:2)+60];

%stim_location_deltas: first two elements of first row of TiledStimulusRect
%minus median stim_location x and y

if (computerLocation < 10) || isfield(trialData{1}, 'm_aiTiledStimulusRect')
    stim_location_deltas = cellfun(@(x) x(1,1:2) - stim_location(1,1:2), vars.m_aiTiledStimulusRect(isTrialOfInterest), 'UniformOutput',false);
else
    stim_location_deltas = cellfun(@(x) x(1,1:2) - stim_location(1,1:2), TiledStimulusRect(isTrialOfInterest), 'UniformOutput',false);
end
stim_location_deltas = vertcat(stim_location_deltas{:});

%stimscale
stimscale = (stim_location(3)-stim_location(1))/60;

%stimtype
stimtype = [DualstimPrimaryuseRGBCloudPerFrame{isTrialOfInterest}];

%trialID
trialID = [TrialIDPerFrame{isTrialOfInterest}];

%trial_start_ts
trial_start_ts = [stimStartTimes(isTrialOfInterest)]';

%useLeye
useLeye = [UseLeyePerFrame{isTrialOfInterest}];

%useReye
useReye = [UseReyePerFrame{isTrialOfInterest}];

%valid_data
totalFrames = sum(numFrames(isTrialOfInterest));
tvec=1:totalFrames;
Block_offsetinds = block_inds(2,:);

bad_inds_block=sort([block_inds(1,:), block_inds(1,:)+1, block_inds(1,:)+2,block_inds(1,:)+3,block_inds(1,:)+4,block_inds(1,:)+5,block_inds(1,:)+6],1);

bad_inds_fix = unique([find(abs(ETtrace(1,:))>ETdist_thresh),find(abs(ETtrace(2,:))>ETdist_thresh)]);
bad_inds_all = unique([bad_inds_block,bad_inds_fix,bad_inds_fix-1, bad_inds_fix-2, bad_inds_fix+1, bad_inds_fix+2]); % remove indices immediately preceding and following eye movement artifacts
use_inds_fix=setdiff(tvec,bad_inds_all);

% find sequences less than 10 due to eye movement removal, and exclude them
% to avoid clogging up the modeling pipeline with tiny snippets
diffs=diff([1,use_inds_fix]);
[~,X]=find(diff(diffs)<10);
for k= X, use_inds_fix(diffs(k):diffs(k+1)-1)=0; end
use_inds_fix(use_inds_fix==0)=[];

valid_data = use_inds_fix;

% added by me: ks batch
spikeSortingBatch = vertcat(SU_ks_batch{:});
spikeSortingBatchMU = vertcat(MU_ks_batch{:});

cluster = vertcat(SU_clusterIDs{:});
clusterMU = vertcat(MU_clusterIDs{:});

%% remap cluster ids from 1 to number of clusters
% [clusterIDs_sorted,sortByClusterID] = sort(clusterIDs);
% spike_ts_sorted = spike_ts(sortByClusterID);
%
% spikeIDs = clusterIDs_sorted;

%uniqueClusterIDs = vertcat(allUnit_clusterIDs{:});
uniqueClusterIDs = vertcat(vertcat(SU_clusterIDs{:}), vertcat(MU_clusterIDs{:}));

spikeIDs_unsorted = accumarray(transpose(1:numel(clusterIDs)), clusterIDs', [], @(x) find(uniqueClusterIDs==x));
[spikeIDs, sortBySpikeID] = sort(spikeIDs_unsorted);
spike_ts_sorted = spike_ts(sortBySpikeID);

% Process reward times relative to trial onsets -- oops this is kofiko time
%REWon_ts = g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == 32765);
%REWoff_ts = g_strctDAQParams.LastStrobe.TimeStamp(g_strctDAQParams.LastStrobe.Buffer == 32764);

%% Make list of reward times on each trial
reward_on_ts = [];
reward_off_ts = [];

if computerLocation < 10
    [~, ~, rew_off_bin] = histcounts(REW_OFF_ts, stimIntervals);
    REW_OFF_cellArray = accumarray(rew_off_bin(:)+1,...
        REW_OFF_ts(:),...
        [nBins + 1,1],...
        @(x){x},...
        {[]});

    [~, ~, rew_on_bin] = histcounts(REW_ON_ts, stimIntervals);
    REW_ON_cellArray = accumarray(rew_on_bin(:)+1,...
        REW_ON_ts(:),...
        [nBins + 1,1],...
        @(x){x},...
        {[]});


    REW_ON_cellArray = cellfun(@(rewardTime, trialStartTime) rewardTime - trialStartTime, REW_ON_cellArray, REW_OFF_cellArray)) 


    reward_on_ts = vertcat(REW_ON_cellArray{2*find(isTrialOfInterest)});
    reward_off_ts = vertcat(REW_OFF_cellArray{2*find(isTrialOfInterest)});
else
    for tr = 1:length(trial_start_ts)
        t = trial_start_ts(tr);
        % Identify reward times for the current trial
        rewardON = REW_ON_ts(REW_ON_ts >= t & REW_ON_ts < t+4.0);
        rewardOFF = REW_OFF_ts(REW_OFF_ts >= t & REW_OFF_ts < t+4.0);
        for ii=1:length(rewardON)
            reward_on_ts(end+1) = rewardON(ii)-t + 4*(tr-1);
        end
        for ii=1:length(rewardOFF)
            reward_off_ts(end+1) = rewardOFF(ii)-t + 4*(tr-1);
        end
    end
end

%% Add fields to data struct
data.ETgains = ETgains;
data.ETstim_location = ETstim_location;
data.ETtrace = ETtrace;
data.ETtrace_raw = ETtrace_raw;
data.Robs = RobsSU;
data.RobsMU = RobsMU;
data.RobsMU_probe_ID = RobsMU_probe_ID;
data.RobsMU_rating = RobsMU_rating;
data.Robs_probe_ID = Robs_probe_ID;
data.Robs_rating = [];
data.blockID = blockID;
data.block_inds = block_inds;
data.cloud_area = cloud_area;
data.cloud_binary = cloud_binary;
data.cloud_scale = cloud_scale;
data.datafilts = datafilts;
data.datafiltsMU = datafiltsMU;
data.dt = dt;
data.electrode_info = electrode_info;
data.exptdate = exptdate;
data.exptname = exptname;
data.fix_location = fix_location;
data.fix_size = fix_size;
data.pixel_size = pixel_size;
data.sacc_inds = sacc_inds;
data.spikeIDs = spikeIDs;
data.spike_ts = spike_ts_sorted;
data.stim = stim;
data.stimET = stimET;
data.stim_area = stim_area;
data.stim_location = stim_location;
data.stim_location_deltas = stim_location_deltas;
data.stimscale = stimscale;
data.stimtype = stimtype;
data.trialID = trialID;
data.trial_start_ts = trial_start_ts;
data.useLeye = useLeye;
data.useReye = useReye;
data.reward_on_ts = reward_on_ts;
data.reward_off_ts = reward_off_ts;

data.valid_data = valid_data;

data.spikeSortingBatch = spikeSortingBatch;
data.spikeSortingBatchMU = spikeSortingBatchMU;

data.ks_folders = ks_folders;
data.chan_offsets = chan_offsets;
data.cluster_offsets = cluster_offsets;

%% modify array_labels so its python readable
unique_array_labels = unique(array_labels);

[C, IA, IC] = unique(array_labels);
data.arrayPerSU = IC(data.spikeSortingBatch);
data.arrayPerMU = IC(data.spikeSortingBatchMU);

%array_labels = cellfun(@(num,lab) [num2str(num) lab], num2cell(1:num_ks_batch), array_labels, 'UniformOutput', false);
%data.array_labels = horzcat(array_labels{:});
%data.array_labels = array_labels;
unique_labels = horzcat(cellfun(@(num,lab) ['_' num2str(num) lab], num2cell(1:length(C)), unique_array_labels, 'UniformOutput', false));
data.array_labels = horzcat(unique_labels{:});

data.cluster = cluster;
data.clusterMU = clusterMU;


%% Saving
if saving
    disp('Saving')
    if ~isdir(savepath)
        mkdir(savepath);
    end

    switch targ_ETstimtype
        case 0; curETstimtype = 'NA';
        case 1; curETstimtype='1D';
        case 7; curETstimtype='CC';
    end

    switch_stimtype = unique(vertcat(vars.DualstimPrimaryuseRGBCloud{isTrialOfInterest}));
    switch switch_stimtype
        case 0; curstimstype='GT';
        case 3; curstimstype='HL';
        case 6; curstimstype='HC';
        case 8; curstimstype='CC';
    end


    array_label_filepart = [cellfun(@(x) [x '_'], unique_array_labels(1:end-1), 'UniformOutput', false) unique_array_labels(end)];
    array_label_filepart = horzcat(array_label_filepart{:});

    useofflinesorting = 1;  % why is this flag all the way down here?
    %FullExpt_ET_filename = sprintf( '%s_FullExpt_ET.mat', filenameP );
    FullExpt_ET_filename = sprintf( 'K%s_FullExpt_ET.mat', exptdate );

    if computerLocation < 10
        data_filename=[monkey_name '_' exptname(1:6) '_' array_label_filepart '_' curstimstype '_ET' curETstimtype '_v10.mat'];
        fixinfo_filename=[filenameP '_fixinfo.mat'];
    else
        data_filename = ['K' exptdate '_' curstimstype '_ET' curETstimtype '_v10.mat'];
        fixinfo_filename = sprintf( 'K%s_fixinfo.mat', exptdate );
    end
    LFPfilename = sprintf( 'K%s_LFPs.mat', exptdate );


    %data_filename=[ 'K' exptdate '_' curstimstype '_ET' curETstimtype '_v10.mat'];
    %fixinfo_filename=['K' exptdate '_fixinfo.mat'];
    save(fullfile(savepath, data_filename),  '-struct', 'data', '-v7.3'); % save packaged cloud data
    save(fullfile(savepath, fixinfo_filename), '-struct', 'ETdata', '-v7.3') % save fixinfo
    save(fullfile(savepath, FullExpt_ET_filename), 'PlexET_ad_calib', 'PlexET_times', '-v7.3'); % save FullExpt_ET

    if ~skipLFP
        trial_start_inds = floor(trial_start_ts*1000);
        LFP_ad = LFPs{1};
        if length(LFPs) == 1  % then this is the only array
            save(fullfile(savepath, LFPfilename), 'LFP_ad', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
        elseif length(LFPs) == 3  % then old-school (2022 expt date)
            LFPa2 = LFPs{2};
            LFPa3 = LFPs{3};
            %save(fullfile(savepath, [filenameP '_LFP.mat']), 'LFP_ad', 'LFPa2', 'LFPa3', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
            save(fullfile(savepath, LFPfilename), 'LFP_ad', 'LFPa2', 'LFPa3', 'trial_start_ts', 'trial_start_inds', '-v7.3' )
        else
            disp('Have not programmed in this LFP-array config yet')
        end
    end

end


% trial_start_inds = round(trialstart_plx'.*1000);