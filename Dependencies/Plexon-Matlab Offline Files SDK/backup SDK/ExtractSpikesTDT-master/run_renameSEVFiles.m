% require format: RZZZZ_YYYYMMDDX-N_data_ch??.sev
dataDir = '/Users/mattgaidica/Documents/Data/ChoiceTask/R0036/R0036-rawdata/R0036_20150225a/R0036_20150225a'; %CHANGE THIS 1/2
baseName = 'R0036_20150225a_data_'; %CHANGE THIS 2/2
dirFiles = dir(fullfile(dataDir,'*.sev'));

for i=1:length(dirFiles)
    C = strsplit(dirFiles(i).name,'_');
    C = strsplit(C{end},'.'); %C{1} = chXX
    disp(fullfile(dataDir,[baseName C{1} '.sev']))
    
    movefile(fullfile(dataDir,dirFiles(i).name),...
        fullfile(dataDir,[baseName C{1} '.sev']));
end