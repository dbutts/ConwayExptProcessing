function trial = loadKofikoTrialData(kofiko_subfolder,mainKofiko_fname, filenameP)

% formerly returend vars dont think i need it
%% Load kofiko data

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

% Concatenate the g_strctLocalExperimentRecorindg structs and remove empty cells
trial = vertcat(g_strctLocalExperimentRecording{:});
trial(cellfun(@isempty, trial)) = [];

% Unify field names across g_strctLocalExperimentRecording to facilitate using cellfun
fieldNames = cellfun(@(x) fieldnames(x), trial, 'UniformOutput', false);
uniqueFieldNames = unique(vertcat(fieldNames{:}));
missingFieldNames = cellfun(@(x) setdiff(uniqueFieldNames, fieldnames(x)),...
    trial, 'UniformOutput',false);

for i = 1:numel(trial)
    for  j = 1:numel(missingFieldNames{i})
        trial{i}.(missingFieldNames{i}{j}) = [];
    end
end

% Find unique trials by finding unique Flip ON timestamps
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trial);
tempImageFlipON_TS_Kofiko_unique = unique(tempImageFlipON_TS_Kofiko);
[~, uniqueTrialIdx, ~] = intersect(tempImageFlipON_TS_Kofiko, tempImageFlipON_TS_Kofiko_unique);
trial= trial(uniqueTrialIdx);
% Order trials by flip ON
tempImageFlipON_TS_Kofiko = cellfun(@(x) x.m_fImageFlipON_TS_Kofiko, trial);
[~, trial_idx_for_sorting] = sort(tempImageFlipON_TS_Kofiko , 'ascend');
trial = trial(trial_idx_for_sorting);
trial = [trial{:}];

% figure out maximum size of variable in each field
allFieldNames = fieldnames(trial);
for f = 1:numel(allFieldNames)
    try
        temp = cellfun(@size, {trial.(allFieldNames{f})}, 'UniformOutput', false);
        maxSize{f} = max(vertcat(temp{:}),[],1);
    catch
        maxSize{f}= nan;
    end
end

% now for each trial and each field replace empties with nans of
% appropriate size

for t = 1:numel(trial)
    for  f = 1:numel(allFieldNames)
        if isempty(trial(t).(allFieldNames{f}))
            try
            trial(t).(allFieldNames{f}) = nan(maxSize{f});
            catch
                trial(t).(allFieldNames{f}) = nan;
            end
        end
    end
end


end