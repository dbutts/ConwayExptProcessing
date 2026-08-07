function cellArray = binByStimIntervals(t,x, stimIntervals)

[~,~,t_bin] = histcounts(t, stimIntervals);
nBins = numel(stimIntervals) - 1;

cellArray = accumarray( ...
   t_bin(:) + 1, ...
   x(:), ...
   [nBins + 1, 1], ...   % force size (extra 1 for bin 0)
   @(y){y}, ...
   {[]} ...              % fill empty bins with empty cells
   );

end