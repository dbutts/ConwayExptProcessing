%% get indices
% correct and incorrect
Corr_inds=[]; Incorr_inds=[]; Ttasktype=[]; Toutcome=[]; Tcatchtrial=[];
for ii=1:length(ExptTrials)
    Toutcome{ii}=ExptTrials{ii, 1}.m_strctTrialOutcome.m_strResult; % self-explanatory
    if strcmp(ExptTrials{ii, 1}.m_strctTrialOutcome.m_strResult, 'Incorrect')
        Incorr_inds = [Incorr_inds,ii];
    elseif strcmp(ExptTrials{ii, 1}.m_strctTrialOutcome.m_strResult, 'Correct')
        Corr_inds = [Corr_inds,ii];
    end
    Ttasktype(ii) =     ExptTrials{ii, 1}.m_iMkTurkTaskType; % 1 means a color trial, 2 means a shape trial
    Tcatchtrial(ii) =   ExptTrials{ii, 1}.isCatchTrial; % 1 means the task type has just switched
    TcueID(ii) =        ExptTrials{ii, 1}.m_iMkTurkTargetID; % the correct color/shape this trial
    TchoiceIDs(ii,:) =  ExptTrials{ii, 1}.m_iMkTurkChoiceIDs; % the possible color/shape IDs this trial (Inote that the correct one is always in the first index)
end

%% RGB Color info from fnInitializeMkTurkPlusTextures.m
cRED = [246, 142, 209, 121, 142, 81, 83, 48, 162, 93, 236, 136, 191, 110];
cGREEN = [134, 77, 165, 96, 187, 108, 190, 109, 163, 94, 125, 73, 165, 94];
cBLUE = [175, 97, 95, 52, 129, 71, 218, 120, 255, 153, 255, 142, 198, 109];
cAll=[cRED;cGREEN;cBLUE];
%%
[ld,rg,yv] = rgb2ldrgyv(cAll./255);
cDKLAll = [ld;rg;yv];
