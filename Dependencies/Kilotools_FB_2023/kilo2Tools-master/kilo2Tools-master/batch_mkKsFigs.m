%% batch_mkKsFigs
%
% make a bunch of figures for each kiloSort direcory (ksDir)
% should be done AFTER manula curation stage in Phy

%% add necessary paths to toolboxes:
cd('D:\Code\Toolboxes\kilo2Tools')
paths = addPathsForSpikeSorting;

%% list of paths to raw ephys files
ksDirList = {...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150717a_msc_pre'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150717b_msc_post'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150719a_msc_pre'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150719b_msc_post'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150721a_msc_pre'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150721b_msc_post'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150724a_msc_pre'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150724b_msc_post'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150726a_msc_pre'...
    %     'Y:\LAB PROJECTS\fefAttention\sna20150726b_msc_post'...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171127a-msc-mergedPreAndPost',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20170921a-msc-mergedPreAndPost',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171021a-msc-mergedPreAndPost',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171214a-msc-mergedPreAndPost',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171129a-msc-mergedPreAndPost',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171202a-msc-mergedPreAndPost',...
%     'Y:\LAB PROJECTS\scVprobeAfterPandemic\data\pop20200902',...
        'Y:\LAB PROJECTS\scVprobeAfterPandemic\data\pop20200914',...
    % 'Y:\LAB PROJECTS\fstAttention_mergedPreAndPost\s20171208a-msc-mergedPreAndPost',...
    };

nFiles = numel(ksDirList);

%% (1) load sp, (2) mkfigs, (3) convert to su, (4) mkfigs

for iF = 1:nFiles
    tic
    ksDir = fullfile(ksDirList{iF}, 'kiloSorted2');
    
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~~~~  ' ksDir])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % create sp:
    sp = getSp(ksDir, 'waves', true);
    
    % make sure there's a figures folder:
    if ~exist(fullfile(ksDir, 'figures'), 'dir')
        mkdir(fullfile(ksDir, 'figures'));
    end
    
    % set options:
    opts.saveFigs = 1;
    opts.dirFigs = ksDir;
    
    % make dem figures and save:
    mkfig.waveform_overChannels(sp, opts);
    mkfig.waveformOverChannels_perCluster(sp, opts);
    mkfig.waveformAndSpikeCount_overChannels(sp, opts)
    su = sp2su(sp, ksDir);
    mkfig.unitSummary(su, opts);
    
    close all
    timePassed(iF) = toc;
end


%%


