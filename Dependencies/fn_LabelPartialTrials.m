function [ExptTrials] = fn_LabelPartialTrials(ExptTrials,opts)
%fn_LabelPartialTrials Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(opts, 'trl_fix_thresh')
    opts.trl_fix_thresh = 1; % default) indicates fraction of a trial the monkey has to fixate in order for the trial to be included in analysis
end

    numTrials = size(ExptTrials,1);

    for iTrials=1:numTrials
        % try and check for partial trials
        xtrace = (ExptTrials{iTrials,1}.m_afEyeXPositionScreenCoordinates-960);
        ytrace = (ExptTrials{iTrials,1}.m_afEyeYPositionScreenCoordinates-540);  
    %        bad_inds = find(abs(xtrace)>45 | abs(ytrace)>45);
    %        use_inds_trial = ones(1,length(xtrace)); use_inds_trial(bad_inds)=0;
        xfix = length(find(abs(xtrace)<45))./length(xtrace);
        yfix = length(find(abs(ytrace)<45))./length(ytrace);
    
        if xfix>opts.trl_fix_thresh && yfix>opts.trl_fix_thresh;
            ExptTrials{iTrials,1}.m_bMonkeyFixatedOverride=1;
        else
            ExptTrials{iTrials,1}.m_bMonkeyFixatedOverride=0;
        end
    end

end