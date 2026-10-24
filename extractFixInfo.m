function ETdata = extractFixInfo(trial, stimTiming)
calibrationTasks = {'Fivedot','FiveDot', 'Dotgrid'};
isCalibrationTrial = cellfun(@(x) any(strcmpi(x, calibrationTasks)), {trial.m_strTrialType});

ETdata.fixloc = vertcat(trial(isCalibrationTrial).m_pt2iFixationSpot);
ETdata.plxonset = stimTiming.stimStartTimes(isCalibrationTrial);

end