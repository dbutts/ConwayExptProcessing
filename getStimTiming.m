function stimTiming = getStimTiming(trial, B)

if nargin < 2
    B = ones(2,1);
end

stimStartTimes = [ones(size([trial.m_fImageFlipON_TS_Kofiko]')),...
    [trial.m_fImageFlipON_TS_Kofiko]']*B;

stimStopTimes = stimStartTimes + [trial.m_fStimulusON_MS]'/1e3;

stimIntervals = [stimStartTimes stimStopTimes]';
stimIntervals = stimIntervals(:);

% Number of frames per trial
numFrames =  transpose(min([trial.numFrames], [trial.numFrames] ./ [trial.repframes]));

stimTiming.stimStartTimes = stimStartTimes;
stimTiming.stimStopTimes = stimStopTimes;
stimTiming.stimIntervals = stimIntervals;
stimTiming.numFrames = numFrames;


end

