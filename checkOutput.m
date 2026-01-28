function [] = checkOutput(outputdir,filenameP)
%%% check to see if all the analysis files are saved


etna = exist([outputdir filenameP '_laminar_CC_ETNA_v09.mat'], 'file'); % need to edit how the file is saved so that it includes time
fullet = exist([outputdir filenameP filesep 'Analysis' filesep filenameP '_FullExpt_ET.mat'], 'file');
fixInfo = exist([outputdir filenameP filesep 'Analysis' filesep filenameP '_fixinfo.mat'], 'file');

fprintf('\n\n**CheckOutput()**\n')
disp('All all of the analysis files saved?');
disp([outputdir filenameP '_laminar_CC_ETNA_v09.mat:'])
disp(string(etna))
disp([outputdir filenameP filesep filenameP '_FullExpt_ET.mat:'])
disp(string(fullet))
disp([outputdir filenameP filesep 'Analysis' filesep filenameP '_fixinfo.mat:'])
disp(string(fixInfo))

fprintf('\nDoes fixinfo have the necessary information?\n')
 % Load metadata from the .mat file
info = whos('-file', [outputdir filenameP filesep 'Analysis' filesep filenameP '_fixinfo.mat']);

% If no variables exist, it's clearly empty
if isempty(info)
    disp('It is empty')
end

% Check each variable's size
for i = 1:length(info)
    % Check if the variable has non-zero size
    if any(info(i).size > 0)
       disp([info(i).name 'is not empty'])
    end
end





