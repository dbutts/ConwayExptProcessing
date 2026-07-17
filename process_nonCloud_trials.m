% disc probe analysis

%isDiscProbeTrial = cellfun(@(x) strcmpi(x, 'Disc Probe'), vars.m_strTrialType);
dpcolors_RGB = vertcat(vars.DiscprobeColor{:});


unique_dpcolors_RGB = unique(dpcolors_RGB(all(isfinite(dpcolors_RGB), 2),:), 'rows');
si
% for each color, get trial indices where it was shown 

% psth matrix: columns are color, rows are unit, elements give nu

clusterIDs = vertcat(allUnit_clusterIDs{:});

% make cell array where each row is trial or non trial interval, each
% element contains cluster IDs corresponding to spikes in that trial
C = horzcat(clusterIDForEachSpk_cellArray{:});
C = cellfun(@(x) vertcat(x{:}), num2cell(C,2), 'UniformOutput', false);

S = horzcat(spk_times_cellArray{:});
S = cellfun(@(x) vertcat(x{:}), num2cell(S,2), 'UniformOutput', false);

% stim on times
%S = S(2:2:end);
%C = C(2:2:end);

temp = [0;repelem(stimStartTimes,2)];
temp = temp(1:end-1);


S = cellfun(@(x,y) x - y, S, num2cell(temp), 'UniformOutput', false);

for color = 1:size(unique_dpcolors_RGB,1)
    trialIndices = 2*find(all(dpcolors_RGB == unique_dpcolors_RGB(color,:),2));
    trialIndices = sort([trialIndices; trialIndices+1]);
    trialIndices(trialIndices > size(S,1)) = [];

    temp = S(trialIndices);
    for cluster = 1:numel(clusterIDs)
        psth{cluster, color} = cellfun(@(X,y) y(X == clusterIDs(cluster)), C(trialIndices), temp, 'UniformOutput', false); 
    end
end

figure; hold on
temp = psth{45,144};
for i = 1:size(temp,1)/2
     scatter(1000*vertcat(temp{2*i-1:2*i}), i, 'ko', 'filled')
end
xlim([0 300])