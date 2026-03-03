function setupStrct = preprocessingSetup

% 3/2/26  -- MJG: wrote it
setupFields = {'userName', 'monkey_name', 'exptDate', 'filenameP', 'dataSource'};

setupDlgArray = {'inputdlg','User name', 'Your name:', {'David Hubel'},...
    'inputdlg', 'Monkey name', 'Enter monkey name:', {'Sprout'},...
    'inputdlg', 'Experiment date', 'Enter experiment date:', {'260101'},...
    'inputdlg', 'File name', 'Enter file name:', {'260101_123456_Sprout'},...
    'questdlg', 'Date source', 'Read data from local machine or isilon?', {'Local'}, {'Local', 'Isilon', 'Cancel'}};

setupStrct = struct;
fieldCounter = 1;
for i = 1:numel(setupDlgArray)

    if strcmpi(setupDlgArray{i}, 'inputdlg')
        dlgType = setupDlgArray{i};
        title = setupDlgArray{i+1};
        prompt = setupDlgArray{i+2};
        default = setupDlgArray{i+3};

        answer =  dlgFun(dlgType,title, prompt, default);

        try
            answer = answer{:};
        catch
        end

        setupStrct.(setupFields{fieldCounter}) = answer;
        fieldCounter = fieldCounter + 1;

    elseif strcmpi(setupDlgArray{i}, 'questdlg')
        dlgType = setupDlgArray{i};
        title = setupDlgArray{i+1};
        prompt = setupDlgArray{i+2};
        default = setupDlgArray{i+3};
        btns = setupDlgArray{i+4};

        answer = dlgFun(dlgType, title, prompt, default, btns);

        try
            answer = answer{:};
        catch
        end

        setupStrct.(setupFields{fieldCounter}) = answer;
        fieldCounter = fieldCounter+1;
    end
end

disp('Setup starting')

% set isilon path
if isunix && ~ismac
    mountBase = '/mnt';
elseif ismac
    mountBase = [filesep 'Volumes'];
elseif ispc
    isilonDriveLetter = inputdlg('Isilon Drive Letter:');
    mountBase = [isilonDriveLetter{1} ':\'];
end
isilonPath = fullfile(mountBase, 'isilon');

if strcmpi(setupStrct.(setupFields{strcmpi(setupFields, 'dataSource')}), 'Local')
    homeDir = char(java.lang.System.getProperty('user.home'));
    processingPath = fullfile(homeDir, 'Git', 'ConwayExptProcessing');

    if ~isdir(processingPath)
        oldProcessingPath = processingPath;
        processingPath = fullfile(isilonPath, 'code', 'ConwayExptProcessing');
        warning('No such directory as %s \nSetting processingPath to %s', oldProcessingPath, processingPath);
    end

    dirpath = fullfile(homeDir, 'Data', 'V1_Fovea', setupStrct.monkey_name, setupStrct.exptDate);

    if ~isdir(dirpath)
        warning('No such directory as %s \n',  dirpath);
    end

    stimpath = fullfile(homeDir, 'Processing', 'Cloudstims_calib_04_2024/');

    if ~isdir(stimpath)
        oldStimpath = stimpath;
        stimpath = fullfile(isilonPath, 'PROJECTS', 'V1_Fovea', 'stimuli', 'Cloudstims_calib_04_2024');
        warning('No such directory as %s \nSetting dirPath and pl2path to %s', oldStimpath, stimpath);
    end
else
    processingPath = fullfile(isilonPath, 'code', 'ConwayExptProcessing');
    dirpath = fullfile(isilonPath, 'DATA', 'monkey_ephys');
    stimpath = fullfile(isilonPath, 'PROJECTS', 'V1_Fovea', 'stimuli', 'Cloudstims_calib_04_2024');
end

addpath(processingPath); % add necessary dependencies
addpath(genpath(fullfile(processingPath, 'Dependencies')));
addpath(fullfile(processingPath, 'Tools'));
try
    pl2dir = dir(fullfile(dirpath, '**/*.pl2'));
    pl2folders = {pl2dir.folder};
    pl2files = {pl2dir.name};

    pl2path = pl2folders{cellfun(@(x) contains(x, setupStrct.filenameP), pl2files)};
    dirpath = pl2path;
    plexon_fname = fullfile(pl2path, [setupStrct.filenameP '.pl2']);
catch
    pl2path ='';
    plexon_fname = '';
end
outputdir = fullfile(homeDir, 'Data', 'V1_Fovea');
disp('setup complete')


%% multiple arrays

% dialog to set up arrays

arrayLabels = dlgFun('inputdlg', 'Array labels',  {'Array lablel(s):'}, {'UT1, UT2, lam'});
arrayLabels = cellfun(@strip, split(arrayLabels, ","), 'UniformOutput', false); % split comma separated values and remove whitespace
arrayLabels = transpose(arrayLabels);


nChans = dlgFun('inputdlg', 'Array channels', cellfun(@(x) [x ' channels:'], arrayLabels, 'UniformOutput', false));

try
    nChans = cellfun(@str2num, nChans);
    nChans = transpose(nChans);


    if numel(nChans) > 1
        chnOffsets = cumsum([0 nChans(1:end-1)]);
    else
        chnOffsets = 0;
    end
catch
    warning('Issue with array channel numbers')
    chnOffsets = [];
end


arraySpacing = dlgFun('inputdlg', 'Array spacing',cellfun(@(x) [x ' spacing in um (x,y):'], arrayLabels, 'UniformOutput', false));
arraySpacing = regexp(arraySpacing, '\d*', 'match');
arraySpacing = transpose(cellfun(@str2num, vertcat(arraySpacing{:})));


curChannels = dlgFun('inputdlg', 'Current channels', cellfun(@(x) [x ' Current channels (e.g. 1:64):'], arrayLabels, 'UniformOutput', false));
curChannels = regexp(curChannels, '\d*:\d*', 'match');

preconverted = dlgFun('questdlg', 'Preconversion', 'Have the data been preconverted (raw to binary/.dat)?', {'No'}, {'Yes', 'No', 'Cancel'});

setupStrctFilename = fullfile(dirpath, ['setupStrct' datestr(now, 'mm_dd_yyyy_HH_MM_SS') '.mat']);

setupStrct.processingPath = processingPath;
setupStrct.dirpath = dirpath;
setupStrct.stimpath = stimpath;
setupStrct.pl2path = pl2path;
setupStrct.plexon_fname = plexon_fname;
setupStrct.outputdir = outputdir;
setupStrct.arrayLabels = arrayLabels;
setupStrct.nChans = nChans;
setupStrct.chnOffsets = chnOffsets;
setupStrct.arraySpacing = arraySpacing;
setupStrct.curChannels = curChannels;
setupStrct.preconverted = preconverted;
setupStrct.filename = setupStrctFilename;

try
    save(setupStrctFilename,'-struct', setupStrct);

catch
    setupStrctFilename = fullfile(['setupStrct_' datestr(now, 'mm_dd_yyyy_HH_MM_SS') '.mat']);
    save(setupStrctFilename,'-struct', 'setupStrct');
    warning('setupStrct saved in present working directory')
end

function answer = dlgFun(dlgType,title, prompt, default, btns)
attempts = 0;
answer = {};

if nargin < 4
    default = cellfun(@num2str, cell(size(prompt)), 'UniformOutput', false);
end

switch dlgType
    case 'inputdlg'
        while attempts < 3 & isempty(answer)
            answer = inputdlg(prompt, title, 1, default);
            attempts = attempts + 1;
        end
    case 'questdlg'
        while attempts < 3 & isempty(answer)
            answer = questdlg(prompt, title, btns{1}, btns{2}, btns{3}, default);
            attempts = attempts+1;
        end
end

if attempts == 3
  error('User quit')  

end