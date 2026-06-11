mainKofiko_fname = input('Main kofiko file: ');
kofiko_subfolder = input('Kofiko subfolder: ');

[~, filenameP] = fileparts(mainKofiko_fname);

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


ScreenSizeX_pix = g_strctStimulusServer.m_aiScreenSize(3);
ScreenSizeY_pix = g_strctStimulusServer.m_aiScreenSize(4);

trialTypeOfInterestIdx = strcmpi( vars.m_strTrialType, 'Dual Stim');

cloudIdx= trialTypeOfInterestIdx & vertcat(vars.DualstimPrimaryuseRGBCloud{:})==8;

cloud_rect = mode(vertcat(vars.m_aiStimulusRect{cloudIdx}));

cloud_center_xy = [cloud_rect(1) + (cloud_rect(3)-cloud_rect(1))/2, ...
                cloud_rect(4) - (cloud_rect(4)-cloud_rect(2))/2];

ecc_x = cloud_center_xy(1) - ScreenSizeX_pix/2;
ecc_y = cloud_center_xy(2) - ScreenSizeY_pix/2;

sqrt(ecc_x^2 + ecc_y^2)/60;

fprintf('ecc_x: %d, ecc_y: %d \n', ecc_x, ecc_y)