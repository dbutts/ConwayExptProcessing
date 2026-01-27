function exportSessionConf(sessionName,saveDir)
% This function exports a configuration file for a session so
% that processing can be offloaded to a non-networked machine.
% INPUTS:
%   sessionName : ex. R0036_20150225a
%   saveDir: where you want this config file save, right now the
%   extractSpikesTDT script prompts for the location of this file

[~,ratID] = sql_getSubjectFromSession(sessionName);
chMap = sql_getChannelMap(ratID);
% you can pass tetrodes into extractSpikesTDT, assume you can't
% do it here (the DB needs to be populated)
tetrodeList = chMap.tetNames;
validMasks = sql_getAllTetChannels(sessionName);

filename = ['session_conf_',sessionName,'.mat'];
save(fullfile(saveDir,filename),...
    'chMap','ratID','tetrodeList','validMasks');