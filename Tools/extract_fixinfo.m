function ETdata = extract_fixinfo( exptdata, metadata_struct, targ_trialtype, output_directory)
% usage: ETdata = extract_fixinfo( ExptTrials, ExptInfo, 'Dotgrid', outputdir)

% extract fixation dot and relevant eyetracking information for a specific 
% trial type (e.g. Fivedot, Dotgrid) to use for calibration

ETdata.fixloc   = [];
ETdata.plxonset = [];
%ETdata.fixsize  = [];

ii=1;
for tt=1:length(exptdata)
    if any(strcmp(exptdata{tt,1}.m_strTrialType, targ_trialtype))
        ETdata.fixloc(ii,:) = exptdata{tt,1}.m_pt2iFixationSpot;
        ETdata.plxonset(ii,:) = exptdata{tt,1}.PlexonOnsetTime;
        %ETdata.fixsize(ii,:) = exptdata{tt,1}.m_fFixationSizePix;  
        ii=ii+1;
    end
end

cur_filename=[metadata_struct.exptname '_fixinfo.mat'];
cd(output_directory)
save(cur_filename, '-struct', 'ETdata', '-v7.3')

disp('Done saving ET info for calibration.') 

return