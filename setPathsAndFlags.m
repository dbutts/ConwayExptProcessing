%function setPathsAndFlags

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


%end