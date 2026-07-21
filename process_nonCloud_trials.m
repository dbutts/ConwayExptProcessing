% make cell array where each row is trial or non trial interval, each
% element contains cluster IDs corresponding to spikes in that trial

v_bkg = 0.5.*ones(3,1);
C = horzcat(clusterIDForEachSpk_cellArray{:});
C = cellfun(@(x) vertcat(x{:}), num2cell(C,2), 'UniformOutput', false);

S = horzcat(spk_times_cellArray{:});
S = cellfun(@(x) vertcat(x{:}), num2cell(S,2), 'UniformOutput', false);

% get spike times relative to trial start time 
temp = [0;repelem(stimStartTimes,2)];
temp = temp(1:end-1);
S = cellfun(@(x,y) x - y, S, num2cell(temp), 'UniformOutput', false);

%
%rows: clusters; columns: trials
clusterIDs = vertcat(allUnit_clusterIDs{:});
R = cell(numel(clusterIDs), numel(S)/2);
for i = 1:numel(S)/2
    c = C{2*i : 2*i+1};
    s = S{2*i : 2*i+1};
    for k = 1:numel(clusterIDs)
        R{k,i} = s(c == clusterIDs(k));
    end
end


%% disc probe analysis
firingRateMatrix = cellfun(@(x) numel(x(x >= 0.05 & x <= 0.15))/0.1, R);

% for each unit
 
% for each color

% get mean firing rate in 50-150 ms window after stim onset 


% vertically concatenate disc probe RGB values across trials
dpcolors_RGB = vertcat(vars.DiscprobeColor{:});
unique_dpcolors_RGB = unique(dpcolors_RGB(all(isfinite(dpcolors_RGB), 2),:), 'rows');

dpcolors_dkl =(T_DKL2RGB \ (((unique_dpcolors_RGB'./255) - v_bkg)./v_bkg));
sat = sqrt(dpcolors_dkl(2,:).^2 + dpcolors_dkl(3,:).^2);
hueAngles = rad2deg(atan2(dpcolors_dkl(3,:),dpcolors_dkl(2,:)));

for color = 1:size(unique_dpcolors_RGB,1)
    trialIndices = all(dpcolors_RGB == unique_dpcolors_RGB(color,:),2);
    dpTuningMatrix(:,color) = mean(firingRateMatrix(:,trialIndices),2);
end


%figure, histogram(1e3.*vertcat(R{117,[2*find(isDiscProbeTrial); 2*find(isDiscProbeTrial)+1]}), 'BinEdges', binEdges)


% for color = 1:size(unique_dpcolors_RGB,1)
%     trialIndices = 2*find(all(dpcolors_RGB == unique_dpcolors_RGB(color,:),2));
%     trialIndices = sort([trialIndices; trialIndices+1]); % stim on and subsequent stim off period
%     trialIndices(trialIndices > size(S,1)) = [];
% 
%     temp = S(trialIndices);
%     for cluster = 1:numel(clusterIDs)
%         %raster{cluster, color} = cellfun(@(X,y) y(X == clusterIDs(cluster)), C(trialIndices), temp, 'UniformOutput', false); 
%         raster{cluster, color} = 
%     end
% end
% 
% % plot psths across colors
% for i = 1:size(raster, 1)
%     temp = vertcat(raster{i,:});
%     temp = vertcat(temp{:});
% 
%     binEdges = 0:10:400;
%     figure, histogram(1000.*temp, 'BinEdges', binEdges, 'Normalization', 'countdensity')
% end