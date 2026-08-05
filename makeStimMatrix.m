function [stim1_cellArray, stim2_cellArray, stim3_cellArray] = makeStimMatrix(stimpath, trial, LumScale)

if nargin < 3
    LumScale = 1;
end

cache = containers.Map('KeyType','char','ValueType','any');

nTrials = numel(trial);

primaryStimPresent =  [trial.m_aiStimulusArea] > 2;
secondaryStimPresent = [trial.m_aiSecondaryStimulusArea] > 2;

stim1_cellArray = cell(nTrials,1);
stim2_cellArray = cell(nTrials,1);
stim3_cellArray = cell(nTrials,1);

stimseq = cellfun(@(x,y) x(1:y:end), {trial.stimseq}, num2cell(cellfun(@(x) max(0, x), {trial.repframes})), 'UniformOutput', false);
stimseq_ET_Cclouds = cellfun(@(x,y) x(1:y:end), {trial.stimseq_ET_Cclouds}, num2cell(cellfun(@(x) max(0, x), {trial.repframes})), 'UniformOutput', false);
stimseq_ET_Aclouds = cellfun(@(x,y) x(1:y:end), {trial.stimseq_ET_Aclouds}, num2cell(cellfun(@(x) max(0, x), {trial.repframes})), 'UniformOutput', false);
stimseq_ET_bars = cellfun(@(x,y) x(1:y:end), {trial.stimseq_ET_bars}, num2cell(cellfun(@(x) max(0, x), {trial.repframes})), 'UniformOutput', false);

% for each trial, if there is a primary stimulus, figure out its type
% (ccloud, harlety, etc), and parameters (spatial scale, etc.), load the
% images and put them into stimulus_cellArray

% if there is a secondary stimulus, load all possible ones (cclouds,
% aclouds, etc) and fill out corresponding cell arrays

% Load primary stimuli

for t = 1:nTrials
    %% DUAL STIM
    if strcmpi(trial(t).m_strTrialType, 'Dual Stim')
        if primaryStimPresent(t)
            % imScalar =  fix(trial(t).m_aiStimulusArea/ 60);
            % imScalar = max([imScalar,1]);

            switch trial(t).DualstimPrimaryuseRGBCloud
                case 0 % ground truth
                case 1 % bars (Felix says to ignore)
                case 2 % bars (Felix says to ignore)
                case 3 % achromatic hartleys
                case 4 % L-M hartleys
                case 5 % S hartleys
                case 6 % Full hartleys

                    [stim1_cellArray(t), cache] = loadHartleys_helper(stimpath, stimseq(t), cache);

                case 7 % Achromatic clouds
                case 8 % Chromatic clouds
                    [stim1_cellArray(t), cache] = loadCclouds_helper( ...
                        trial(t), stimpath, stimseq(t), cache, LumScale);
            end

            if secondaryStimPresent(t)
                % imScalar =  fix(trial(t).m_aiSecondaryStimulusArea / 60);
                % imScalar = max([imScalar,1]);

                try
                    % note: spatialscale is that of color clouds, not achrom clouds it
                    % seems
                    switch trial(t).DualstimSecondaryUseCloud
                        case 0 % secondary stim is vertical bars, tertiary horiztonal
                        case 1  % both stims alternatve between vertical adn horizontal bars
                        case 2  % secondary stim alternates between vertical and horizontal bars, tertiary stim is achrom cloud
                        case 3  % secondary stim is achrom cloud, tertiary vertical/horizontal bars
                        case 4  % both stims achrom clouds
                        case 5  % secondary stim vertical/horizontal bars, tertiary stim ccloud

                            % SECONDARY STIMULUS

                            % TERTIARY STIMULUS

                            [stim3_cellArray(t), cache] = loadCclouds_helper( ...
                                trial(t), stimpath, stimseq_ET_Cclouds(t), cache, LumScale);

                        case 6  % secondary stim ccloud, tertiary stim vertical/horizontal bars

                            % SECONDARY STIMULUS

                            [stim2_cellArray(t), cache] = loadCclouds_helper( ...
                                trial(t), stimpath, stimseq_ET_Cclouds(t), cache, LumScale);
                            % TERTIARY STIMULUS

                        case 7  % both stim cclouds

                            % SECONDARY STIMULUS
                            [stim2_cellArray(t), cache] = loadCclouds_helper( ...
                                trial(t), stimpath, stimseq_ET_Cclouds(t), cache, LumScale);

                            % TERTIARY STIMULUS

                            stim3_cellArray(t) = stim2_cellArray(t);

                        case 8  % secondary stim achrom cloud, tertiary stim ccloud

                            % SECONDARY STIMULUS

                            % TERTIARY STIMULUS
                            [stim3_cellArray(t), cache] = loadCclouds_helper( ...
                                trial(t), stimpath, stimseq_ET_Cclouds(t), cache, LumScale);
                    end
                catch ME
                    disp(ME)
                end
            end
        end
    elseif strcmpi(trial(t).m_strTrialType, 'Disc Probe')
      stim1_cellArray{t} = transpose(trial(t).DiscprobeColor);
    end
    % reshape
    % stim1_cellArray = cellfun(@(x) permute(x, [3 1 2 4]), stim1_cellArray, 'UniformOutput', false);
    % stim1_cellArray = cellfun(@(x) reshape(x, size(x,1), prod(size(x, 2:4))), stim1_cellArray, 'UniformOutput', false);
    % stim1_cellArray = cellfun(@transpose, stim1_cellArray, 'UniformOutput', false);
    % %stim1_matrix = horzcat(stim1_cellArray{isTrialOfInterest});
    %
    % stim2_cellArray = cellfun(@(x) permute(x, [3 1 2 4]), stim2_cellArray, 'UniformOutput', false);
    % stim2_cellArray = cellfun(@(x) reshape(x, size(x,1), prod(size(x, 2:4))), stim2_cellArray, 'UniformOutput', false);
    % stim2_cellArray = cellfun(@transpose, stim2_cellArray, 'UniformOutput', false);
    % %stim2_matrix = horzcat(stim2_cellArray{isTrialOfInterest});
    %
    % stim3_cellArray = cellfun(@(x) permute(x, [3 1 2 4]), stim3_cellArray, 'UniformOutput', false);
    % stim3_cellArray = cellfun(@(x) reshape(x, size(x,1), prod(size(x, 2:4))), stim3_cellArray, 'UniformOutput', false);
    % stim3_cellArray = cellfun(@transpose, stim3_cellArray, 'UniformOutput', false);

    %stim3_matrix = horzcat(stim3_cellArray{isTrialOfInterest});
end
end